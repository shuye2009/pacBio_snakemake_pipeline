#!/usr/bin/env python
"""
Plot heatmap of DSS-called differentially methylated regions.
Reads DSS DMR BED file directly (no filtering) and extracts methylation
values from pb-cpg-tools combined.bed.gz files for each region.
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pysam
import os


def load_dmr_regions_from_bed(bed_file):
    """Load DMR regions from a BED file.
    
    Returns list of region dicts with chrom, start, end, name.
    """
    regions = []
    with open(bed_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = f"{chrom}:{start}-{end}"
            regions.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': name
            })
    print(f"Loaded {len(regions)} DSS DMR regions from BED file")
    return regions


def extract_methylation_for_regions(bed_gz_file, regions):
    """Extract mean methylation for each region from a combined.bed.gz file."""
    region_means = {}
    try:
        tbx = pysam.TabixFile(bed_gz_file)
    except Exception as e:
        print(f"Warning: Could not open tabix file {bed_gz_file}: {e}")
        return {r['name']: np.nan for r in regions}
    try:
        for r in regions:
            meth_values = []
            try:
                for row in tbx.fetch(r['chrom'], r['start'], r['end']):
                    if row.startswith('#'):
                        continue
                    fields = row.split('\t')
                    if len(fields) >= 4:
                        try:
                            meth = float(fields[3]) / 100.0
                            meth_values.append(meth)
                        except ValueError:
                            continue
            except ValueError:
                pass
            if meth_values:
                region_means[r['name']] = np.mean(meth_values)
            else:
                region_means[r['name']] = np.nan
    finally:
        tbx.close()
    return region_means


def main():
    parser = argparse.ArgumentParser(description='Generate heatmap for DSS DMRs')
    parser.add_argument('--dmr-bed', required=True, help='Path to DSS DMR BED file')
    parser.add_argument('--bed-dir', required=True, help='Directory containing combined.bed.gz files')
    parser.add_argument('--case-samples', required=True, help='Comma-separated case sample names')
    parser.add_argument('--control-samples', required=True, help='Comma-separated control sample names')
    parser.add_argument('--output-png', required=True, help='Output PNG file path')
    parser.add_argument('--output-pdf', required=True, help='Output PDF file path')
    args = parser.parse_args()

    case_samples = args.case_samples.split(",")
    control_samples = args.control_samples.split(",")

    regions = load_dmr_regions_from_bed(args.dmr_bed)

    if len(regions) == 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, 'No DSS DMRs found',
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()
        return

    meth_data = {}
    for sample in case_samples + control_samples:
        bed_file = f"{args.bed_dir}/{sample}.combined.bed.gz"
        if os.path.exists(bed_file):
            print(f"Processing {sample}...")
            meth_data[sample] = extract_methylation_for_regions(bed_file, regions)
        else:
            print(f"Warning: BED file not found for {sample}: {bed_file}")

    if not meth_data:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, 'No methylation data available', ha='center', va='center', transform=ax.transAxes)
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()
        return

    meth_matrix = pd.DataFrame(meth_data)
    meth_matrix = meth_matrix.dropna(how='all')
    nan_count = meth_matrix.isna().sum().sum()
    if nan_count > 0:
        print(f"Warning: {nan_count} missing methylation values imputed to 0.5 "
              f"({nan_count / (len(meth_matrix) * len(meth_matrix.columns)) * 100:.1f}% of matrix)")
    meth_matrix = meth_matrix.fillna(0.5)

    print(f"Created matrix with {len(meth_matrix)} regions x {len(meth_matrix.columns)} samples")

    meth_matrix.index = [str(idx)[:50] if len(str(idx)) > 50 else str(idx) for idx in meth_matrix.index]

    sample_colors = pd.Series(
        ['#e74c3c' if s in case_samples else '#3498db' for s in meth_matrix.columns],
        index=meth_matrix.columns,
        name='Group'
    )

    n_regions = len(meth_matrix)
    fig_height = max(8, min(50, n_regions * 0.3))

    if len(meth_matrix) > 1 and len(meth_matrix.columns) > 1:
        g = sns.clustermap(
            meth_matrix,
            cmap='RdYlBu_r',
            center=0.5,
            vmin=0, vmax=1,
            col_colors=sample_colors,
            row_cluster=True,
            col_cluster=True,
            figsize=(max(14, len(meth_matrix.columns) * 1.5 + 2), fig_height + 2),
            dendrogram_ratio=(0.15, 0.2),
            cbar_kws={'label': 'Methylation Level'},
            xticklabels=True,
            yticklabels=True,
            tree_kws={'linewidths': 1.5},
        )
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
        g.fig.suptitle(f'DSS DMRs (n={n_regions})', y=1.02)
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='#e74c3c', label='Case'),
                         Patch(facecolor='#3498db', label='Control')]
        g.ax_heatmap.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1.15))
        g.savefig(args.output_png, dpi=300, bbox_inches='tight')
        g.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()
    else:
        fig, ax = plt.subplots(figsize=(max(12, len(meth_matrix.columns) * 1.5), fig_height))
        sns.heatmap(meth_matrix, cmap='RdYlBu_r', center=0.5,
                   vmin=0, vmax=1, ax=ax,
                   xticklabels=True, yticklabels=True,
                   cbar_kws={'label': 'Methylation Level'})
        plt.title(f'DSS DMRs (n={n_regions})')
        plt.xlabel('Samples (Red=Case, Blue=Control)')
        plt.ylabel('DMR Regions')
        plt.xticks(rotation=45, ha='right')
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()

    print(f"Saved heatmap with {n_regions} DSS DMR regions")


if __name__ == '__main__':
    main()
