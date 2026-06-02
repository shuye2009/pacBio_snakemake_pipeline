#!/usr/bin/env python3
"""
Plot density distributions of mod_score, cov, and est_mod_count.

Produces:
  - Per-sample: one PNG/PDF with 3 subplots (one per feature)
  - Composite: one PNG/PDF with a 3×N grid (features × samples) overlaid
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde


def load_bed(path):
    """Load pb-CpG-tools combined BED, skipping comment lines."""
    df = pd.read_csv(
        path, sep='\t', comment='#', header=None,
        names=["chrom", "begin", "end", "mod_score", "type",
               "cov", "est_mod_count", "est_unmod_count", "discretized_mod_score"],
        compression='gzip'
    )
    return df


def plot_per_sample(bed_path, sample_name, output_png, output_pdf):
    """Generate density plot for a single sample with 3 subplots."""
    df = load_bed(bed_path)
    features = ["mod_score", "cov", "est_mod_count"]
    labels = ["Modification Score", "Coverage", "Estimated Modified Count"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle(f"Density Distributions — {sample_name}", fontsize=13)

    for ax, feat, label in zip(axes, features, labels):
        data = df[feat].dropna()
        if feat in ("cov", "est_mod_count"):
            data = np.log2(data + 1)
            xlabel = f"log2({label})"
        else:
            xlabel = label
        if len(data) > 10:
            kde = gaussian_kde(data, bw_method=0.3)
            x_grid = np.linspace(data.min(), data.max(), 300)
            ax.plot(x_grid, kde(x_grid), linewidth=1.5, color="steelblue")
            ax.fill_between(x_grid, kde(x_grid), alpha=0.2, color="steelblue")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Density")
        ax.set_title(label)
        # Add median line
        med = data.median()
        ax.axvline(med, color='red', linestyle='--', linewidth=1, label=f"median={med:.1f}")
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close()


def plot_composite(bed_paths, sample_names, case_samples, control_samples,
                   output_png, output_pdf):
    """Generate composite density plot overlaying all samples for each feature."""
    features = ["mod_score", "cov", "est_mod_count"]
    labels = ["Modification Score", "Coverage", "Estimated Modified Count"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("Density Distributions — All Samples", fontsize=13)

    case_set = set(case_samples)
    control_set = set(control_samples)

    for ax, feat, label in zip(axes, features, labels):
        all_vals = []
        for bed_path, sample in zip(bed_paths, sample_names):
            df = load_bed(bed_path)
            data = df[feat].dropna()
            if feat in ("cov", "est_mod_count"):
                data = np.log2(data + 1)
            all_vals.append(data)
        # Shared x range across samples
        x_min = min(d.min() for d in all_vals)
        x_max = max(d.max() for d in all_vals)
        x_grid = np.linspace(x_min, x_max, 300)
        for data, sample in zip(all_vals, sample_names):
            color = "tab:red" if sample in case_set else "tab:blue"
            if len(data) > 10:
                kde = gaussian_kde(data, bw_method=0.3)
                ax.plot(x_grid, kde(x_grid), linewidth=1.2, color=color,
                        alpha=0.7, label=sample)
        if feat in ("cov", "est_mod_count"):
            ax.set_xlabel(f"log2({label})")
        else:
            ax.set_xlabel(label)
        ax.set_ylabel("Density")
        ax.set_title(label)

    # Create legend with group labels
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='tab:red', linewidth=2, label='Case'),
        Line2D([0], [0], color='tab:blue', linewidth=2, label='Control'),
    ]
    axes[-1].legend(handles=legend_elements, loc='upper right', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close()


def plot_per_chromosome(bed_path, sample_name, output_png, output_pdf):
    """Generate a line plot (KDE) of mod_score for each chromosome."""
    df = load_bed(bed_path)

    # Get sorted chromosome list (natural sort: chr1, chr2, ..., chr22, chrX, chrY)
    chroms = df["chrom"].unique()
    def chrom_sort_key(c):
        c_stripped = c.replace("chr", "")
        if c_stripped.isdigit():
            return (0, int(c_stripped))
        return (1, c_stripped)
    chroms = sorted(chroms, key=chrom_sort_key)

    # Determine grid layout
    n_chroms = len(chroms)
    ncols = 6
    nrows = int(np.ceil(n_chroms / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 2.5), squeeze=False)
    fig.suptitle(f"mod_score Density per Chromosome \u2014 {sample_name}", fontsize=13)

    x_grid = np.linspace(0, 100, 200)

    for i, chrom in enumerate(chroms):
        row, col = divmod(i, ncols)
        ax = axes[row][col]
        data = df.loc[df["chrom"] == chrom, "mod_score"].dropna()
        if len(data) > 10:
            kde = gaussian_kde(data, bw_method=0.3)
            ax.plot(x_grid, kde(x_grid), linewidth=1.2, color="steelblue")
            ax.fill_between(x_grid, kde(x_grid), alpha=0.2, color="steelblue")
        ax.set_title(chrom, fontsize=9)
        ax.set_xlim(0, 100)
        ax.set_xlabel("mod_score", fontsize=7)
        ax.set_ylabel("Density", fontsize=7)
        ax.tick_params(labelsize=6)

    # Hide unused axes
    for i in range(n_chroms, nrows * ncols):
        row, col = divmod(i, ncols)
        axes[row][col].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Plot density distributions of CpG metrics")
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # Per-sample mode
    sp = subparsers.add_parser("per-sample")
    sp.add_argument("--bed", required=True, help="Path to combined.bed.gz")
    sp.add_argument("--sample", required=True, help="Sample name")
    sp.add_argument("--output-png", required=True)
    sp.add_argument("--output-pdf", required=True)

    # Per-chromosome mode
    pc = subparsers.add_parser("per-chromosome")
    pc.add_argument("--bed", required=True, help="Path to combined.bed.gz")
    pc.add_argument("--sample", required=True, help="Sample name")
    pc.add_argument("--output-png", required=True)
    pc.add_argument("--output-pdf", required=True)

    # Composite mode
    cp = subparsers.add_parser("composite")
    cp.add_argument("--beds", required=True, nargs='+', help="Paths to combined.bed.gz files")
    cp.add_argument("--samples", required=True, nargs='+', help="Sample names (same order as beds)")
    cp.add_argument("--case-samples", required=True, nargs='+', help="Case sample names")
    cp.add_argument("--control-samples", required=True, nargs='+', help="Control sample names")
    cp.add_argument("--output-png", required=True)
    cp.add_argument("--output-pdf", required=True)

    args = parser.parse_args()

    if args.mode == "per-sample":
        plot_per_sample(args.bed, args.sample, args.output_png, args.output_pdf)
    elif args.mode == "per-chromosome":
        plot_per_chromosome(args.bed, args.sample, args.output_png, args.output_pdf)
    elif args.mode == "composite":
        plot_composite(args.beds, args.samples, args.case_samples, args.control_samples,
                       args.output_png, args.output_pdf)


if __name__ == "__main__":
    main()
