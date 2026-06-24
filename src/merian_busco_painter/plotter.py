"""Plot BUSCO locations colored by Merian element assignments."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager

BAR_HEIGHT = 0.62
LABEL_OFFSET_FACTOR = 0.02
LABEL_PADDING_FACTOR = 0.06
DEFAULT_PANEL_SIZE = 20
PLOT_BODY_WIDTH_CM = 22
ROW_HEIGHT_CM = 1.10
MIN_PLOT_HEIGHT_CM = 18
MAX_PANEL_COLUMNS = 3
COMPACT_LABEL_THRESHOLD = 80
ASSEMBLY_MODES = ("auto", "final", "draft")
DEFAULT_LABEL_WRAP = 4


def resolve_open_sans_font(env_var: str = "GENOMENOTES_FONT") -> str | None:
    """Locate a regular Open Sans font file."""
    import os

    explicit = os.environ.get(env_var)
    if explicit and Path(explicit).is_file():
        return explicit

    def pick_upright(paths: list[Path]) -> str | None:
        regular = [path for path in paths if "Regular" in str(path)]
        if regular:
            return str(sorted(regular)[0])
        upright = [path for path in paths if "italic" not in path.name.lower()]
        if upright:
            return str(sorted(upright)[0])
        return str(sorted(paths)[0]) if paths else None

    user_fonts = Path.home() / "Library" / "Fonts"
    chosen = pick_upright(list(user_fonts.glob("OpenSans*.ttf")))
    if chosen:
        return chosen

    package_root = Path(__file__).resolve().parent.parent
    font_dir = package_root / "assets" / "fonts"
    if font_dir.is_dir():
        hits: list[Path] = []
        for pattern in ("OpenSans-Regular.ttf", "OpenSans*.ttf", "open-sans*.ttf"):
            hits.extend(font_dir.glob(pattern))
        chosen = pick_upright(hits)
        if chosen:
            return chosen

    return None


def setup_font() -> None:
    """Configure Open Sans when available, otherwise use matplotlib defaults."""
    try:
        font_path = resolve_open_sans_font()
        if font_path:
            font_manager.fontManager.addfont(font_path)
            plt.rcParams["font.family"] = "Open Sans"
            plt.rcParams["font.style"] = "normal"
            plt.rcParams["font.weight"] = "normal"
            print(f"[INFO] Using Open Sans font: {font_path}")
        else:
            print("[WARN] Open Sans not found, using default font")
    except Exception as exc:
        print(f"[WARN] Could not load Open Sans: {exc}")


def detect_lengths_format(lengths_file: Path) -> str:
    """Return 'final' for Chrom/Length_Mb tables, or 'draft' for .fai files."""
    with open(lengths_file) as fh:
        first_line = fh.readline().strip().split("\t")

    if first_line[:2] == ["Chrom", "Length_Mb"]:
        return "final"
    if len(first_line) >= 2:
        return "draft"
    raise ValueError(f"Could not detect lengths file format: {lengths_file}")


def load_lengths(lengths_file: Path, assembly_mode: str = "auto") -> pd.DataFrame:
    """Load either final chrom_lengths.tsv or draft .fai lengths into bp."""
    if assembly_mode not in ASSEMBLY_MODES:
        raise ValueError(f"Unknown assembly mode: {assembly_mode}")

    detected = detect_lengths_format(lengths_file)
    if assembly_mode != "auto" and assembly_mode != detected:
        expected = "Chrom/Length_Mb TSV" if assembly_mode == "final" else ".fai"
        raise ValueError(
            f"{lengths_file} looks like {detected!r} lengths, but "
            f"--assembly-mode {assembly_mode} expects {expected}"
        )

    if detected == "final":
        chrom_lengths = pd.read_csv(lengths_file, sep="\t")
        chrom_lengths["length"] = chrom_lengths["Length_Mb"] * 1e6
        return chrom_lengths.rename(columns={"Chrom": "query_chr"})[
            ["query_chr", "length"]
        ]

    return pd.read_csv(
        lengths_file,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["query_chr", "length"],
    )


def load_data(
    location_file: Path, lengths_file: Path | None = None, assembly_mode: str = "auto"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load BUSCO locations and chromosome/scaffold lengths."""
    locations = pd.read_csv(location_file, sep="\t", keep_default_na=False)
    locations["query_chr"] = locations["query_chr"].str.replace(":.*", "", regex=True)
    locations["position"] = pd.to_numeric(locations["position"], errors="coerce")

    if lengths_file:
        chrom_lengths = load_lengths(lengths_file, assembly_mode=assembly_mode)
    else:
        if assembly_mode in {"final", "draft"}:
            print(
                "[WARN] No lengths file supplied; estimating lengths from BUSCO "
                f"positions despite --assembly-mode {assembly_mode}"
            )
        chrom_lengths = locations.groupby("query_chr")["position"].max().reset_index()
        chrom_lengths["length"] = chrom_lengths["position"] * 1.05

    return locations, chrom_lengths


def filter_chromosomes(
    locations: pd.DataFrame, minimum_buscos: int = 3
) -> pd.DataFrame:
    """Filter chromosomes by minimum BUSCO count, keeping accession placeholders."""
    is_placeholder = locations["buscoID"] == "NA"
    busco_counts = (
        locations[~is_placeholder]
        .groupby("query_chr")
        .size()
        .reset_index(name="n_busco")
    )
    valid_chroms = busco_counts[busco_counts["n_busco"] >= minimum_buscos][
        "query_chr"
    ]
    placeholder_chroms = locations.loc[is_placeholder, "query_chr"]
    keep_chroms = set(valid_chroms).union(placeholder_chroms)
    return locations[locations["query_chr"].isin(keep_chroms)].copy()


def format_merian_label(merians: list[str], wrap: int = DEFAULT_LABEL_WRAP) -> str:
    """Return a Merian label, optionally wrapped by element count."""
    if wrap <= 0:
        return "; ".join(merians)

    wrapped_lines = [
        "; ".join(merians[index : index + wrap])
        for index in range(0, len(merians), wrap)
    ]
    return "\n".join(wrapped_lines)


def calculate_merian_labels(
    locations: pd.DataFrame, threshold: int = 5, wrap: int = DEFAULT_LABEL_WRAP
) -> dict[str, str]:
    """Return Merian labels ordered by median BUSCO position per chromosome."""
    counts = (
        locations.groupby(["query_chr", "assigned_chr"])
        .agg(n=("assigned_chr", "size"), median_position=("position", "median"))
        .reset_index()
    )
    counts = counts[counts["n"] >= threshold]

    def merian_sort_key(merian: str) -> tuple[bool, int]:
        return (merian != "MZ", int(merian[1:]) if merian != "MZ" else 0)

    labels: dict[str, str] = {}
    for chrom in counts["query_chr"].unique():
        chrom_counts = counts[counts["query_chr"] == chrom].copy()
        chrom_counts["merian_sort"] = chrom_counts["assigned_chr"].map(
            merian_sort_key
        )
        chrom_counts = chrom_counts.sort_values(
            ["median_position", "merian_sort"], kind="stable"
        )
        labels[chrom] = format_merian_label(
            chrom_counts["assigned_chr"].tolist(), wrap=wrap
        )

    return labels


def get_merian_colors() -> dict[str, object]:
    """Colorblind-friendly palette for Merian elements."""
    return {
        "M1": (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
        "M2": (0.6823529411764706, 0.7803921568627451, 0.9098039215686274),
        "M3": (1.0, 0.4980392156862745, 0.054901960784313725),
        "M4": (1.0, 0.7333333333333333, 0.47058823529411764),
        "M5": (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
        "M6": (0.596078431372549, 0.8745098039215686, 0.5411764705882353),
        "M7": (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
        "M8": (1.0, 0.596078431372549, 0.5882352941176471),
        "M9": (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
        "M10": (0.7725490196078432, 0.6901960784313725, 0.8352941176470589),
        "M11": (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
        "M12": (0.7686274509803922, 0.611764705882353, 0.5803921568627451),
        "M13": (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
        "M14": (0.9686274509803922, 0.7137254901960784, 0.8235294117647058),
        "M15": (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
        "M16": (0.7803921568627451, 0.7803921568627451, 0.7803921568627451),
        "M17": (0.0, 0.502, 0.502),
        "M18": (0.5, 0.75, 0.75),
        "M19": (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
        "M20": (0.6196078431372549, 0.8549019607843137, 0.8980392156862745),
        "M21": (0.2235294117647059, 0.23137254901960785, 0.4745098039215686),
        "M22": (0.3215686274509804, 0.32941176470588235, 0.6392156862745098),
        "M23": (0.4196078431372549, 0.43137254901960786, 0.8117647058823529),
        "M24": (0.611764705882353, 0.6196078431372549, 0.8705882352941177),
        "M25": (0.38823529411764707, 0.4745098039215686, 0.2235294117647059),
        "M26": (0.5490196078431373, 0.6352941176470588, 0.3215686274509804),
        "M27": (0.7098039215686275, 0.8117647058823529, 0.4196078431372549),
        "M28": (0.807843137254902, 0.8588235294117647, 0.611764705882353),
        "M29": (0.5490196078431373, 0.42745098039215684, 0.19215686274509805),
        "M30": (0.7411764705882353, 0.6196078431372549, 0.2235294117647059),
        "M31": (0.9058823529411765, 0.7294117647058823, 0.3215686274509804),
        "MZ": (0.25, 0.25, 0.25, 1.0),
    }


def get_merian_colors_spectrum() -> dict[str, object]:
    turbo = plt.cm.turbo
    colors: dict[str, object] = {}
    for index, merian in enumerate(f"M{i}" for i in range(1, 32)):
        colors[merian] = turbo(0.05 + (index / 30) * 0.90)
    colors["MZ"] = (0.25, 0.25, 0.25)
    return colors


def get_merian_colors_merianbow() -> dict[str, str]:
    hex_colors = [
        "#666666",
        "#D50062",
        "#F10059",
        "#FF104E",
        "#FF4F44",
        "#FF7839",
        "#FF9E2F",
        "#FF9700",
        "#CF8F00",
        "#9C8600",
        "#6A7A00",
        "#336C00",
        "#0A8500",
        "#009F00",
        "#00B846",
        "#00D27C",
        "#00EBB3",
        "#00D5BF",
        "#00BEC9",
        "#00A6D1",
        "#008FD5",
        "#0079D4",
        "#008BFC",
        "#009BFF",
        "#00A8FF",
        "#00B3FF",
        "#87BCFF",
        "#B198FF",
        "#C971FF",
        "#D546E1",
        "#D500B0",
        "#CC007E",
    ]
    colors = {"MZ": hex_colors[0]}
    for index in range(1, 32):
        colors[f"M{index}"] = hex_colors[index]
    return colors


def get_merian_colors_merianbow4() -> dict[str, str]:
    hex_colors = [
        "#666666",
        "#710093",
        "#BC007B",
        "#EA005B",
        "#fc1d1d",
        "#FD9514",
        "#E9CB19",
        "#87BF13",
        "#00AB3E",
        "#00f2a1",
        "#005c66",
        "#00589E",
        "#006DDB",
        "#0080ff",
        "#A676FF",
        "#ff1aef",
        "#FF82CD",
        "#FF6B70",
        "#EE6A15",
        "#A16C00",
        "#4E6400",
        "#005200",
        "#00e251",
        "#009286",
        "#00B0E0",
        "#00C8FF",
        "#A1D0FF",
        "#ABA9E5",
        "#AE83B6",
        "#A56183",
        "#8E454F",
        "#6B3122",
    ]
    colors = {"MZ": hex_colors[0]}
    for index in range(1, 32):
        colors[f"M{index}"] = hex_colors[index]
    return colors


def get_palette(palette: str) -> dict[str, object]:
    if palette == "spectrum":
        print("[INFO] Using colorblind-friendly spectrum palette")
        return get_merian_colors_spectrum()
    if palette == "merianbow":
        print("[INFO] Using MerianBow palette (original)")
        return get_merian_colors_merianbow()
    if palette == "merianbow4":
        print("[INFO] Using MerianBow4 palette (CVD-optimised)")
        return get_merian_colors_merianbow4()
    print("[INFO] Using categorical palette")
    return get_merian_colors()


def plot_merian_chromosomes(
    locations: pd.DataFrame,
    chrom_lengths: pd.DataFrame,
    output_prefix: str,
    minimum_buscos: int = 3,
    palette: str = "categorical",
    label_threshold: int = 5,
    panel_size: int = DEFAULT_PANEL_SIZE,
    max_columns: int = MAX_PANEL_COLUMNS,
    label_wrap: int = DEFAULT_LABEL_WRAP,
) -> None:
    """Create the main Merian plot with chromosome labels."""
    setup_font()

    locations = filter_chromosomes(locations, minimum_buscos)
    is_placeholder = locations["buscoID"] == "NA"
    valid_merians = ["MZ"] + [f"M{i}" for i in range(1, 32)]
    locations = locations[
        is_placeholder | locations["assigned_chr"].str.upper().isin(valid_merians)
    ].copy()
    locations["assigned_chr"] = locations["assigned_chr"].str.upper()

    if locations.empty:
        raise ValueError(
            f"No chromosomes/scaffolds have at least {minimum_buscos} BUSCOs"
        )

    plotted_chroms = set(locations["query_chr"].dropna().unique())
    chrom_lengths = chrom_lengths[chrom_lengths["query_chr"].isin(plotted_chroms)].copy()
    if chrom_lengths.empty:
        raise ValueError("No plotted chromosomes/scaffolds have matching lengths")

    label_locations = locations[locations["buscoID"] != "NA"].copy()
    merian_labels = calculate_merian_labels(
        label_locations, threshold=label_threshold, wrap=label_wrap
    )
    chrom_order = chrom_lengths.sort_values("length", ascending=False)[
        "query_chr"
    ].tolist()

    n_chroms = len(chrom_order)
    print(
        f"[INFO] Plotting {len(locations)} BUSCOs across "
        f"{n_chroms} chromosomes/scaffolds after filtering..."
    )

    merian_colors = get_palette(palette)
    panel_size = max(1, int(panel_size))
    max_columns = max(1, int(max_columns))
    ncols = min(max_columns, max(1, math.ceil(n_chroms / panel_size)))
    chroms_per_panel = max(1, math.ceil(n_chroms / ncols))
    print(
        f"[INFO] Layout: {ncols} column(s), up to "
        f"{chroms_per_panel} chromosomes/scaffolds per column"
    )

    compact_layout = n_chroms > COMPACT_LABEL_THRESHOLD
    label_fontsize = 8 if compact_layout else 10
    panel_height = max(
        MIN_PLOT_HEIGHT_CM, ROW_HEIGHT_CM * min(chroms_per_panel, max(1, n_chroms))
    )

    panel_chroms_list: list[list[str]] = []
    panel_limits: list[float] = []
    for col_idx in range(ncols):
        start_idx = col_idx * chroms_per_panel
        end_idx = min((col_idx + 1) * chroms_per_panel, n_chroms)
        panel_chroms = chrom_order[start_idx:end_idx]
        panel_chroms_list.append(panel_chroms)
        if panel_chroms:
            panel_max_length = chrom_lengths[
                chrom_lengths["query_chr"].isin(panel_chroms)
            ]["length"].max()
            panel_limits.append(panel_max_length * (1 + LABEL_PADDING_FACTOR))
        else:
            panel_limits.append(1.0)

    fig_width = PLOT_BODY_WIDTH_CM / 2.54
    panel_widths = [1] * ncols if compact_layout else panel_limits
    fig, axes = plt.subplots(
        nrows=1,
        ncols=ncols,
        figsize=(fig_width, panel_height / 2.54),
        squeeze=False,
        gridspec_kw={"width_ratios": panel_widths},
    )
    axes = axes[0]

    for col_idx in range(ncols):
        ax = axes[col_idx]
        panel_chroms = panel_chroms_list[col_idx]
        if not panel_chroms:
            ax.axis("off")
            continue

        y_positions = {chrom: i for i, chrom in enumerate(reversed(panel_chroms))}
        panel_limit = panel_limits[col_idx]

        for chrom in panel_chroms:
            y = y_positions[chrom]
            length = chrom_lengths[chrom_lengths["query_chr"] == chrom][
                "length"
            ].values[0]
            bar_bottom = y - BAR_HEIGHT / 2

            rect = patches.Rectangle(
                (0, bar_bottom),
                length,
                BAR_HEIGHT,
                facecolor="white",
                edgecolor="black",
                linewidth=0.5,
            )
            ax.add_patch(rect)

            chrom_buscos = locations[locations["query_chr"] == chrom]
            for _, busco in chrom_buscos.iterrows():
                if pd.notna(busco["position"]):
                    color = merian_colors.get(busco["assigned_chr"], (0.85, 0.85, 0.85))
                    tile = patches.Rectangle(
                        (busco["position"] - 25000, bar_bottom),
                        50000,
                        BAR_HEIGHT,
                        facecolor=color,
                        edgecolor="none",
                    )
                    ax.add_patch(tile)

            if chrom in merian_labels:
                ax.text(
                    length * (1 + LABEL_OFFSET_FACTOR),
                    y,
                    merian_labels[chrom],
                    va="center",
                    ha="left",
                    fontsize=label_fontsize,
                    color="#333333",
                    linespacing=0.85,
                )

        ax.set_xlim(0, panel_limit)
        ax.set_ylim(-0.6, len(panel_chroms) - 0.4)
        ax.set_xlabel("")
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1e6:.0f}"))
        ax.set_yticks([y_positions[chrom] for chrom in panel_chroms])
        ax.set_yticklabels(panel_chroms, fontsize=label_fontsize)
        ax.set_ylabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    legend_elements = [
        patches.Patch(facecolor=merian_colors[merian], label=merian)
        for merian in ["MZ"] + [f"M{i}" for i in range(1, 32)]
    ]
    fig.supxlabel("Position (Mb)", fontsize=11)
    plt.tight_layout()
    fig.legend(
        handles=legend_elements,
        title="Merian elements",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        fontsize=10,
        title_fontsize=11,
        ncol=2,
        columnspacing=1.2,
    )

    for ext in ("png", "svg"):
        output_file = f"{output_prefix}.{ext}"
        dpi = 300 if ext == "png" else None
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
        print(f"[INFO] Saved: {output_file}")

    plt.close()


def plot_locations(
    location_file: Path,
    output_prefix: str,
    lengths_file: Path | None = None,
    assembly_mode: str = "auto",
    minimum_buscos: int = 3,
    palette: str = "categorical",
    label_threshold: int = 5,
    panel_size: int = DEFAULT_PANEL_SIZE,
    max_columns: int = MAX_PANEL_COLUMNS,
    label_wrap: int = DEFAULT_LABEL_WRAP,
) -> None:
    print("[INFO] Loading data...")
    locations, chrom_lengths = load_data(
        location_file, lengths_file, assembly_mode=assembly_mode
    )
    print(
        f"[INFO] Loaded {len(locations)} BUSCOs and "
        f"{len(chrom_lengths)} chromosome/scaffold lengths..."
    )
    plot_merian_chromosomes(
        locations,
        chrom_lengths,
        output_prefix,
        minimum_buscos=minimum_buscos,
        palette=palette,
        label_threshold=label_threshold,
        panel_size=panel_size,
        max_columns=max_columns,
        label_wrap=label_wrap,
    )
    print("[INFO] Done.")
