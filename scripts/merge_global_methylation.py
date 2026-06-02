#!/usr/bin/env python3
"""
Merge per-sample global methylation TSVs into one matrix and generate
a heatmap and PCA plot.

Inputs: one TSV per sample with columns chrom, bin_start, bin_end, mean_mod_score
Outputs: merged TSV, heatmap (PNG/PDF), PCA plot (PNG/PDF)
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
import os


def merge_samples(tsv_paths, sample_names):
    """Merge per-sample TSVs into a single DataFrame."""
    merged = None
    for f, sample in zip(tsv_paths, sample_names):
        df = pd.read_csv(f, sep='\t')
        df = df.rename(columns={"mean_mod_score": sample})
        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on=["chrom", "bin_start", "bin_end"], how="outer")
    # Natural sort by chromosome: numeric chromosomes first, then non-numeric
    def chrom_sort_key(c):
        c_stripped = c.replace("chr", "")
        if c_stripped.isdigit():
            return (0, int(c_stripped), "")
        return (1, 0, c_stripped)

    merged["_chrom_key"] = merged["chrom"].map(chrom_sort_key)
    merged = merged.sort_values(["_chrom_key", "bin_start"]).reset_index(drop=True)
    merged = merged.drop(columns=["_chrom_key"])
    return merged


def plot_heatmap(merged, sample_names, heatmap_png, heatmap_pdf, bin_size=50000,
                 bin_info=None):
    """Generate per-chromosome heatmaps composited into one figure.

    Each chromosome gets its own subplot with width proportional to its
    genomic length (number of bins). Rows are hierarchically clustered.
    If bin_info is provided, annotation tracks are plotted in separate
    rows below the sample heatmap, each with its own colorbar and scale.
    """
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    from matplotlib.gridspec import GridSpec

    mat = merged[sample_names].copy()
    chroms = merged["chrom"].values

    # Drop bins with all NaN, fill remaining
    valid_mask = ~mat.isna().all(axis=1)
    mat = mat.loc[valid_mask].fillna(mat.loc[valid_mask].mean())
    chroms = chroms[valid_mask]

    # Align bin_info if provided
    annot_cols = []
    annot_data = None
    if bin_info is not None:
        annot_cols = [c for c in bin_info.columns
                      if c not in ("chrom", "bin_start", "bin_end")]
        if annot_cols:
            merged_idx = merged.loc[valid_mask, ["chrom", "bin_start", "bin_end"]].copy()
            merged_idx["_order"] = range(len(merged_idx))
            bi = bin_info.merge(merged_idx, on=["chrom", "bin_start", "bin_end"], how="inner")
            bi = bi.sort_values("_order")
            annot_data = bi[annot_cols].T.values  # n_annotations × n_bins

    # Hierarchical clustering on samples (rows) using genome-wide data
    data = mat.T.values  # shape: (n_samples, n_bins)
    if len(sample_names) > 1:
        dist = pdist(data, metric='euclidean')
        Z = linkage(dist, method='ward')
        row_order = leaves_list(Z)
    else:
        row_order = [0]

    sample_labels = [sample_names[i] for i in row_order]

    # Group bin indices by chromosome
    chrom_order = []
    chrom_bin_counts = []
    chrom_indices = {}
    for i, c in enumerate(chroms):
        if c not in chrom_indices:
            chrom_indices[c] = []
            chrom_order.append(c)
        chrom_indices[c].append(i)
    for c in chrom_order:
        chrom_bin_counts.append(len(chrom_indices[c]))

    # Build GridSpec: width ratios proportional to bin count
    total_bins = sum(chrom_bin_counts)
    width_ratios = [max(1, int(c / total_bins * 100)) for c in chrom_bin_counts]

    n_annot = len(annot_cols)
    n_rows = 1 + n_annot  # sample row + one per annotation
    n_cols = len(chrom_order) + 1  # chromosomes + colorbar column
    height_ratios = [len(sample_names)] + [1] * n_annot
    width_ratios_cbar = width_ratios + [max(1, sum(width_ratios) // 90)]

    fig = plt.figure(figsize=(max(16, len(chrom_order) * 1.8 + 0.8),
                               max(4, sum(height_ratios) * 0.45)))
    gs = GridSpec(n_rows, n_cols, figure=fig,
                  width_ratios=width_ratios_cbar, height_ratios=height_ratios,
                  wspace=0.02, hspace=0.04)

    # Row 0: sample methylation heatmap
    vmin, vmax = 0, 100
    for idx, chrom in enumerate(chrom_order):
        ax = fig.add_subplot(gs[0, idx])
        col_idx = chrom_indices[chrom]
        chrom_data = data[row_order][:, col_idx]

        im = ax.imshow(chrom_data, aspect='auto', cmap='RdYlBu_r',
                       vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_title(chrom, fontsize=7, pad=2)
        ax.set_xticks([])
        if idx == 0:
            ax.set_yticks(range(len(sample_labels)))
            ax.set_yticklabels(sample_labels, fontsize=7)
        else:
            ax.set_yticks([])

    # Methylation colorbar in rightmost column
    cbar_ax = fig.add_subplot(gs[0, -1])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar.set_label("Mod Score", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # Annotation rows with colorbars in rightmost column
    annot_cmaps = ['Greens', 'Oranges', 'Purples', 'Blues', 'Reds', 'Greys']
    for a_idx, col_name in enumerate(annot_cols):
        row = 1 + a_idx
        cmap = annot_cmaps[a_idx % len(annot_cmaps)]
        a_data = annot_data[a_idx]
        a_vmin = 0
        a_vmax = max(a_data.max(), 1)

        for idx, chrom in enumerate(chrom_order):
            ax = fig.add_subplot(gs[row, idx])
            col_idx = chrom_indices[chrom]
            chrom_annot = a_data[col_idx].reshape(1, -1)

            im_a = ax.imshow(chrom_annot, aspect='auto', cmap=cmap,
                             vmin=a_vmin, vmax=a_vmax, interpolation='nearest')
            ax.set_xticks([])
            ax.set_yticks([])
            if idx == 0:
                label = col_name.replace("_count", "")
                ax.set_ylabel(label, fontsize=7, rotation=0,
                              ha='right', va='center', labelpad=20)

        # Colorbar in rightmost column for this annotation
        cbar_ax = fig.add_subplot(gs[row, -1])
        cbar_a = fig.colorbar(im_a, cax=cbar_ax, orientation='vertical')
        cbar_a.ax.tick_params(labelsize=6)

    fig.suptitle(f"Global Methylation Heatmap by Chromosome (mean mod_score per {bin_size // 1000}kb bin)",
                 fontsize=11, y=0.98)
    plt.savefig(heatmap_png, dpi=150, bbox_inches='tight')
    plt.savefig(heatmap_pdf, bbox_inches='tight')
    plt.close()


def plot_pca(merged, sample_names, case_samples, control_samples, pca_png, pca_pdf, bin_size=50000):
    """Generate a PCA plot of samples based on binned methylation."""
    mat = merged[sample_names].copy()

    # Drop bins with all NaN, fill remaining NaN
    valid = mat.dropna(how="all")
    valid = valid.fillna(valid.mean())

    # Transpose: samples as rows, bins as features
    pca_input = valid.T
    pca = PCA(n_components=min(2, len(sample_names)))
    coords = pca.fit_transform(pca_input)

    case_set = set(case_samples)

    fig, ax = plt.subplots(figsize=(7, 6))
    for i, sample in enumerate(sample_names):
        color = "tab:red" if sample in case_set else "tab:blue"
        y = coords[i, 1] if coords.shape[1] > 1 else 0
        ax.scatter(coords[i, 0], y, color=color, s=80, zorder=3)
        ax.annotate(sample, (coords[i, 0], y), fontsize=7, ha='left', va='bottom')

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    if pca.n_components_ > 1:
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    else:
        ax.set_ylabel("PC2")
    ax.set_title(f"PCA of Global Methylation ({bin_size // 1000}kb bins)")

    legend_elements = [
        Patch(facecolor='tab:red', label='Case'),
        Patch(facecolor='tab:blue', label='Control'),
    ]
    ax.legend(handles=legend_elements, loc='best')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(pca_png, dpi=150, bbox_inches='tight')
    plt.savefig(pca_pdf, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-sample global methylation and generate heatmap + PCA"
    )
    parser.add_argument("--tsvs", required=True, nargs='+', help="Per-sample TSV files")
    parser.add_argument("--samples", required=True, nargs='+', help="Sample names (same order as tsvs)")
    parser.add_argument("--case-samples", required=True, nargs='+', help="Case sample names")
    parser.add_argument("--control-samples", required=True, nargs='+', help="Control sample names")
    parser.add_argument("--output-tsv", required=True, help="Merged output TSV")
    parser.add_argument("--output-heatmap-png", required=True)
    parser.add_argument("--output-heatmap-pdf", required=True)
    parser.add_argument("--output-pca-png", required=True)
    parser.add_argument("--output-pca-pdf", required=True)
    parser.add_argument("--bin-size", type=int, default=50000, help="Bin size in bp")
    parser.add_argument("--bin-info", required=False, help="Path to bin_annotation.tsv")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    merged = merge_samples(args.tsvs, args.samples)
    merged.to_csv(args.output_tsv, sep='\t', index=False)

    bin_info = None
    if args.bin_info and os.path.exists(args.bin_info):
        bin_info = pd.read_csv(args.bin_info, sep='\t')

    plot_heatmap(merged, args.samples, args.output_heatmap_png, args.output_heatmap_pdf,
                 args.bin_size, bin_info)
    plot_pca(merged, args.samples, args.case_samples, args.control_samples,
             args.output_pca_png, args.output_pca_pdf, args.bin_size)


if __name__ == "__main__":
    main()
