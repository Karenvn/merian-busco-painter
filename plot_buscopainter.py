#!/usr/bin/env python3
"""Plot BUSCO locations coloured by Merian element assignments."""

import argparse
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import font_manager

def resolve_open_sans_font(env_var="GENOMENOTES_FONT"):
    """Locate a regular Open Sans font file.
    Preference order:
      1. $GENOMENOTES_FONT (explicit override)
      2. ~/Library/Fonts/OpenSans-Regular.ttf (Homebrew / Font Book user install)
      3. Any regular OpenSans*.ttf 
    Returns a string path or None if nothing found.
    """
    import os
    from pathlib import Path

    # 1. explicit override
    p = os.environ.get(env_var)
    if p and Path(p).is_file():
        return p

    # helper to pick first non-italic from list of paths
    def pick_upright(paths):
        # First try to find Regular specifically
        regular = [x for x in paths if "Regular" in str(x)]
        if regular:
            return str(sorted(regular)[0])
        # Otherwise exclude italic
        upright = [x for x in paths if "italic" not in x.name.lower()]
        if upright:
            return str(sorted(upright)[0])
        return str(sorted(paths)[0]) if paths else None

    # 2. user fonts (Homebrew install ends up here)
    user_fonts = Path.home() / "Library" / "Fonts"
    hits = list(user_fonts.glob("OpenSans*.ttf"))
    chosen = pick_upright(hits)
    if chosen:
        return chosen

    # 3. vendored copy in repo assets
    script_root = Path(__file__).resolve().parent.parent
    pkg_dir = script_root / "assets" / "fonts"
    if pkg_dir.is_dir():
        hits = []
        for pat in ("OpenSans-Regular.ttf", "OpenSans*.ttf", "open-sans*.ttf"):
            hits.extend(pkg_dir.glob(pat))
        chosen = pick_upright(hits)
        if chosen:
            return chosen

    return None

def setup_font():
    """Configure OpenSans font if available, fallback to system default"""
    try:
        font_path = resolve_open_sans_font()
        if font_path:
            font_manager.fontManager.addfont(font_path)
            plt.rcParams['font.family'] = 'Open Sans'
            plt.rcParams['font.style'] = 'normal'
            plt.rcParams['font.weight'] = 'normal'
            print(f"[INFO] Using OpenSans font: {font_path}")
        else:
            print("[WARN] OpenSans not found, using default font")
    except Exception as e:
        print(f"[WARN] Could not load OpenSans: {e}")

def load_data(location_file, lengths_file=None):
    """Load BUSCO locations and chromosome lengths"""
    # Read locations
    locations = pd.read_csv(location_file, sep='\t')
    
    # Clean chromosome names (remove any :scaffolds)
    locations['query_chr'] = locations['query_chr'].str.replace(':.*', '', regex=True)
    
    # Read or estimate chromosome lengths
    if lengths_file:
        chrom_lengths = pd.read_csv(lengths_file, sep='\t')
        chrom_lengths['length'] = chrom_lengths['Length_Mb'] * 1e6  # Convert to bp
        chrom_lengths = chrom_lengths.rename(columns={'Chrom': 'query_chr'})
    else:
        # Fallback: estimate from max position + 5% padding
        chrom_lengths = locations.groupby('query_chr')['position'].max().reset_index()
        chrom_lengths['length'] = chrom_lengths['position'] * 1.05
    
    return locations, chrom_lengths

def filter_chromosomes(locations, minimum_buscos=3):
    """Filter chromosomes by minimum BUSCO count"""
    busco_counts = locations.groupby('query_chr').size().reset_index(name='n_busco')
    valid_chroms = busco_counts[busco_counts['n_busco'] >= minimum_buscos]['query_chr']
    return locations[locations['query_chr'].isin(valid_chroms)]

def calculate_merian_labels(locations, threshold=5):
    """
    Calculate Merian element labels for each chromosome.
    Returns {chromosome: 'M1; M3; M18'} for elements with >= threshold BUSCOs.
    """
    # Count BUSCOs per chromosome-Merian combination
    counts = (locations.groupby(['query_chr', 'assigned_chr'])
              .size()
              .reset_index(name='n'))
    
    # Filter by threshold
    counts = counts[counts['n'] >= threshold]
    
    # Group by chromosome and join Merian elements
    labels = {}
    for chrom in counts['query_chr'].unique():
        chrom_merians = counts[counts['query_chr'] == chrom]['assigned_chr'].tolist()
        # Sort Merian elements (MZ first, then M1-M31 numerically)
        sorted_merians = sorted(chrom_merians, 
                               key=lambda x: (x != 'MZ', int(x[1:]) if x != 'MZ' else 0))
        labels[chrom] = '; '.join(sorted_merians)
    
    return labels

def get_merian_colors():
    """Colorblind-friendly palette for Merian elements"""
    colors = {
        'M1': (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
        'M2': (0.6823529411764706, 0.7803921568627451, 0.9098039215686274),
        'M3': (1.0, 0.4980392156862745, 0.054901960784313725),
        'M4': (1.0, 0.7333333333333333, 0.47058823529411764),
        'M5': (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
        'M6': (0.596078431372549, 0.8745098039215686, 0.5411764705882353),
        'M7': (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
        'M8': (1.0, 0.596078431372549, 0.5882352941176471),
        'M9': (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
        'M10': (0.7725490196078432, 0.6901960784313725, 0.8352941176470589),
        'M11': (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
        'M12': (0.7686274509803922, 0.611764705882353, 0.5803921568627451),
        'M13': (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
        'M14': (0.9686274509803922, 0.7137254901960784, 0.8235294117647058),
        'M15': (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
        'M16': (0.7803921568627451, 0.7803921568627451, 0.7803921568627451),
        'M17': (0.0, 0.502, 0.502),  # dark teal
        'M18': (0.5, 0.75, 0.75), # light teal
        'M19': (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
        'M20': (0.6196078431372549, 0.8549019607843137, 0.8980392156862745),
        'M21': (0.2235294117647059, 0.23137254901960785, 0.4745098039215686),
        'M22': (0.3215686274509804, 0.32941176470588235, 0.6392156862745098),
        'M23': (0.4196078431372549, 0.43137254901960786, 0.8117647058823529),
        'M24': (0.611764705882353, 0.6196078431372549, 0.8705882352941177),
        'M25': (0.38823529411764707, 0.4745098039215686, 0.2235294117647059),
        'M26': (0.5490196078431373, 0.6352941176470588, 0.3215686274509804),
        'M27': (0.7098039215686275, 0.8117647058823529, 0.4196078431372549),
        'M28': (0.807843137254902, 0.8588235294117647, 0.611764705882353),
        'M29': (0.5490196078431373, 0.42745098039215684, 0.19215686274509805),
        'M30': (0.7411764705882353, 0.6196078431372549, 0.2235294117647059),
        'M31': (0.9058823529411765, 0.7294117647058823, 0.3215686274509804),
        'MZ': (0.25, 0.25, 0.25, 1.0)
    }
    return colors

def get_merian_colors_spectrum():
    """Generate colorblind-friendly spectrum using perceptually uniform colormaps"""
    turbo = plt.cm.turbo
    merian_elements = [f'M{i}' for i in range(1, 32)]
    colors = {}
    
    for i, merian in enumerate(merian_elements):
        position = 0.05 + (i / 30) * 0.90
        colors[merian] = turbo(position)
    
    colors['MZ'] = (0.25, 0.25, 0.25)
    return colors

def get_merian_colors_merianbow():
    """Original MerianBow palette - full rainbow with regular color space steps"""
    hex_colors = [
        "#666666", "#D50062", "#F10059", "#FF104E", "#FF4F44", "#FF7839", "#FF9E2F",
        "#FF9700", "#CF8F00", "#9C8600", "#6A7A00", "#336C00", "#0A8500", "#009F00",
        "#00B846", "#00D27C", "#00EBB3", "#00D5BF", "#00BEC9", "#00A6D1", "#008FD5",
        "#0079D4", "#008BFC", "#009BFF", "#00A8FF", "#00B3FF", "#87BCFF", "#B198FF",
        "#C971FF", "#D546E1", "#D500B0", "#CC007E"
    ]
    
    colors = {'MZ': hex_colors[0]}
    for i in range(1, 32):
        colors[f'M{i}'] = hex_colors[i]
    
    return colors

def get_merian_colors_merianbow4():
    """MerianBow4 palette - optimized for colorblind accessibility"""
    hex_colors = [
        "#666666", "#710093", "#BC007B", "#EA005B", "#fc1d1d", "#FD9514", "#E9CB19", 
        "#87BF13", "#00AB3E", "#00f2a1", "#005c66", "#00589E", "#006DDB", "#0080ff", 
        "#A676FF", "#ff1aef", "#FF82CD", "#FF6B70", "#EE6A15", "#A16C00", "#4E6400",
        "#005200", "#00e251", "#009286", "#00B0E0", "#00C8FF", "#A1D0FF", "#ABA9E5", 
        "#AE83B6", "#A56183", "#8E454F", "#6B3122"
    ]
    
    colors = {'MZ': hex_colors[0]}
    for i in range(1, 32):
        colors[f'M{i}'] = hex_colors[i]
    
    return colors

def plot_merian_chromosomes(locations, chrom_lengths, output_prefix, minimum_buscos=3,
                           palette='categorical', label_threshold=5, panel_size=40):
    """Create the main Merian plot with chromosome labels"""
    setup_font()
    
    # Filter by minimum BUSCOs, but keep track of all chromosomes from chrom_lengths
    # This ensures sex chromosomes (like W) are always plotted even if empty
    locations_filtered = filter_chromosomes(locations, minimum_buscos)
    
    # Keep track of placeholder rows (where buscoID is NA)
    is_placeholder = locations_filtered['buscoID'] == 'NA'

    # Filter only valid Merian assignments for actual BUSCOs
    valid_merians = ['MZ'] + [f'M{i}' for i in range(1, 32)]
    locations_filtered = locations_filtered[
        is_placeholder | locations_filtered['assigned_chr'].str.upper().isin(valid_merians)
    ]
    locations_filtered['assigned_chr'] = locations_filtered['assigned_chr'].str.upper()
    
    # Use filtered locations for calculating labels, but keep all chromosomes for plotting
    locations = locations_filtered
    
    # Calculate Merian labels for each chromosome
    merian_labels = calculate_merian_labels(locations, threshold=label_threshold)
    
    # Order chromosomes by length (descending) - USE ALL CHROMOSOMES FROM chrom_lengths
    chrom_order = chrom_lengths.sort_values('length', ascending=False)['query_chr'].tolist()
    # DO NOT FILTER - keep all chromosomes even if they have no BUSCOs
    
    n_chroms = len(chrom_order)
    
    # Get colors based on palette choice
    if palette == 'spectrum':
        merian_colors = get_merian_colors_spectrum()
        print("[INFO] Using colourblind-friendly spectrum palette")
    elif palette == 'merianbow':
        merian_colors = get_merian_colors_merianbow()
        print("[INFO] Using MerianBow palette (original)")
    elif palette == 'merianbow4':
        merian_colors = get_merian_colors_merianbow4()
        print("[INFO] Using MerianBow4 palette (CVD-optimised)")
    else:
        merian_colors = get_merian_colors()
        print("[INFO] Using categorical palette")

    panel_size = max(1, int(panel_size))
    ncols = max(1, math.ceil(n_chroms / panel_size))
    panel_height = max(15, 0.9 * min(panel_size, max(1, n_chroms)))
    panel_width = 30 / 2.54

    fig, axes = plt.subplots(
        nrows=1,
        ncols=ncols,
        figsize=(panel_width * ncols, panel_height / 2.54),
        squeeze=False
    )
    axes = axes[0]

    for col_idx in range(ncols):
        ax = axes[col_idx]
        start_idx = col_idx * panel_size
        end_idx = min((col_idx + 1) * panel_size, n_chroms)
        panel_chroms = chrom_order[start_idx:end_idx]

        if not panel_chroms:
            ax.axis('off')
            continue

        y_positions = {chrom: i for i, chrom in enumerate(reversed(panel_chroms))}
        panel_max_length = chrom_lengths[
            chrom_lengths['query_chr'].isin(panel_chroms)
        ]['length'].max()

        for chrom in panel_chroms:
            y = y_positions[chrom]
            length = chrom_lengths[chrom_lengths['query_chr'] == chrom]['length'].values[0]

            rect = patches.Rectangle((0, y - 0.4), length, 0.8,
                                     facecolor='white', edgecolor='black', linewidth=0.5)
            ax.add_patch(rect)

            chrom_buscos = locations[locations['query_chr'] == chrom]
            for _, busco in chrom_buscos.iterrows():
                if pd.notna(busco['position']):
                    color = merian_colors.get(busco['assigned_chr'], (0.85, 0.85, 0.85))
                    tile = patches.Rectangle((busco['position'] - 25000, y - 0.4),
                                             50000, 0.8,
                                             facecolor=color, edgecolor='none')
                    ax.add_patch(tile)

            if chrom in merian_labels:
                label_text = merian_labels[chrom]
                ax.text(length * 1.02, y, label_text,
                        va='center', ha='left', fontsize=10,
                        color='#333333')

        ax.set_xlim(0, panel_max_length * 1.25)
        ax.set_ylim(-0.6, len(panel_chroms) - 0.4)
        ax.set_xlabel('Position (Mb)', fontsize=11)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.0f}'))
        ax.set_yticks([y_positions[c] for c in panel_chroms])
        ax.set_yticklabels(panel_chroms, fontsize=10)
        ax.set_ylabel('')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
    
    # Create legend
    legend_elements = []
    for merian in ['MZ'] + [f'M{i}' for i in range(1, 32)]:
        legend_elements.append(patches.Patch(facecolor=merian_colors[merian], 
                                            label=merian))
    
    axes[-1].legend(handles=legend_elements, title='Merian elements',
                    loc='center left', bbox_to_anchor=(1.01, 0.5),
                    frameon=False, fontsize=10, title_fontsize=11)
    
    plt.tight_layout(rect=[0, 0, 0.82, 1])
    
    # Save outputs
    for ext in ['png', 'svg']:
        output_file = f"{output_prefix}.{ext}"
        dpi = 300 if ext == 'png' else None
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"[INFO] Saved: {output_file}")
    
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Plot BUSCO locations colored by Merian element assignments'
    )
    parser.add_argument('-f', '--file', required=True,
                       help='Path to all_location.tsv file')
    parser.add_argument('-l', '--lengths', default=None,
                       help='Path to chrom_lengths.tsv file (optional)')
    parser.add_argument('-p', '--prefix', default='buscopainter',
                       help='Output prefix for plot files')
    parser.add_argument('-m', '--minimum', type=int, default=3,
                       help='Minimum BUSCOs per chromosome (default: 3)')
    parser.add_argument('--palette', choices=['categorical', 'spectrum', 'merianbow', 'merianbow4'], 
                       default='categorical',
                       help='Color palette: categorical (default), spectrum (turbo), merianbow (original), or merianbow4 (CVD-optimized)')
    parser.add_argument('--label-threshold', type=int, default=5,
                       help='Minimum BUSCOs for a Merian to appear in chromosome label (default: 5)')
    parser.add_argument('--panel-size', type=int, default=40,
                       help='Max chromosomes per panel before splitting into columns (default: 40)')
    
    args = parser.parse_args()
    
    print("[INFO] Loading data...")
    locations, chrom_lengths = load_data(args.file, args.lengths)
    
    print(f"[INFO] Plotting {len(locations)} BUSCOs across {len(chrom_lengths)} chromosomes...")
    plot_merian_chromosomes(locations, chrom_lengths, args.prefix,
                           args.minimum, args.palette, args.label_threshold,
                           args.panel_size)
    
    print("[INFO] Done.")

if __name__ == '__main__':
    main()
