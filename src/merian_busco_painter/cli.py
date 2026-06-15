"""Command line interface for merian-busco-painter."""

from __future__ import annotations

import argparse
from importlib.resources import files
from pathlib import Path

from merian_busco_painter import __version__
from merian_busco_painter.painter import paint_buscos
from merian_busco_painter.plotter import (
    ASSEMBLY_MODES,
    DEFAULT_PANEL_SIZE,
    MAX_PANEL_COLUMNS,
    plot_locations,
)


def default_reference_table() -> Path:
    ref = files("merian_busco_painter").joinpath(
        "data/Merian_elements_full_table.tsv"
    )
    return Path(str(ref))


def configure_paint_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--reference-table",
        "--reference_table",
        "-r",
        dest="reference_table",
        type=Path,
        default=default_reference_table(),
        help="Reference BUSCO-to-Merian table (default: bundled odb10 table)",
    )
    parser.add_argument(
        "--query-table",
        "--query_table",
        "-q",
        dest="query_table",
        type=Path,
        required=True,
        help="BUSCO full_table.tsv for the assembly to paint",
    )
    parser.add_argument(
        "--prefix",
        "-p",
        default="buscopainter",
        help="Output directory-like prefix or filename stem",
    )
    parser.add_argument(
        "--accession",
        "-a",
        help="Assembly accession for NCBI chromosome lengths",
    )
    parser.add_argument(
        "--write-summary",
        "--write_summary",
        dest="write_summary",
        action="store_true",
        help="Write per-chromosome BUSCO count summary",
    )
    parser.set_defaults(func=run_paint)


def add_paint_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "paint", help="Map BUSCO full_table rows to Merian elements"
    )
    configure_paint_parser(parser)
    return parser


def configure_plot_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-f",
        "--file",
        type=Path,
        required=True,
        help="Painter all_location.tsv output",
    )
    parser.add_argument(
        "-l",
        "--lengths",
        type=Path,
        default=None,
        help="Final chrom_lengths.tsv or draft .fai lengths file",
    )
    parser.add_argument(
        "--assembly-mode",
        choices=ASSEMBLY_MODES,
        default="auto",
        help=(
            "auto detects the lengths format; final expects Chrom/Length_Mb TSV; "
            "draft expects .fai"
        ),
    )
    parser.add_argument(
        "-p",
        "--prefix",
        default="buscopainter",
        help="Output prefix for PNG/SVG plot files",
    )
    parser.add_argument(
        "-m",
        "--minimum",
        type=int,
        default=3,
        help="Minimum BUSCOs per chromosome/scaffold (default: 3)",
    )
    parser.add_argument(
        "--palette",
        choices=["categorical", "spectrum", "merianbow", "merianbow4"],
        default="categorical",
        help="Color palette for Merian elements",
    )
    parser.add_argument(
        "--label-threshold",
        type=int,
        default=5,
        help="Minimum BUSCOs for a Merian label to appear (default: 5)",
    )
    parser.add_argument(
        "--panel-size",
        type=int,
        default=DEFAULT_PANEL_SIZE,
        help=(
            "Target chromosomes/scaffolds per panel before splitting into "
            f"columns (default: {DEFAULT_PANEL_SIZE})"
        ),
    )
    parser.add_argument(
        "--max-columns",
        type=int,
        default=MAX_PANEL_COLUMNS,
        help=f"Maximum number of plot columns (default: {MAX_PANEL_COLUMNS})",
    )
    parser.set_defaults(func=run_plot)


def add_plot_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "plot", help="Plot Merian assignments for final or draft assemblies"
    )
    configure_plot_parser(parser)
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="merian-busco-painter",
        description="Paint BUSCO full_table.tsv files with Merian assignments.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_paint_parser(subparsers)
    add_plot_parser(subparsers)
    return parser


def run_paint(args: argparse.Namespace) -> None:
    outputs = paint_buscos(
        reference_table=args.reference_table,
        query_table=args.query_table,
        prefix=args.prefix,
        accession=args.accession,
        write_summary=args.write_summary,
    )
    print("[INFO] Outputs written:")
    print(f"[INFO]   {outputs.all_locations}")
    if outputs.wrote_lengths:
        print(f"[INFO]   {outputs.chrom_lengths}")
    if outputs.wrote_summary:
        print(f"[INFO]   {outputs.summary}")


def run_plot(args: argparse.Namespace) -> None:
    plot_locations(
        location_file=args.file,
        lengths_file=args.lengths,
        assembly_mode=args.assembly_mode,
        output_prefix=args.prefix,
        minimum_buscos=args.minimum,
        palette=args.palette,
        label_threshold=args.label_threshold,
        panel_size=args.panel_size,
        max_columns=args.max_columns,
    )


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


def paint_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="BUSCO to Merian mapper")
    configure_paint_parser(parser)
    args = parser.parse_args(argv)
    run_paint(args)


def plot_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Plot BUSCO Merian assignments")
    configure_plot_parser(parser)
    args = parser.parse_args(argv)
    run_plot(args)
