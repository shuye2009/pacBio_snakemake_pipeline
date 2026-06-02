#!/usr/bin/env python
"""
Plot methylation distribution for each sample as violin plots.
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import numpy as np
import os


def main():
    parser = argparse.ArgumentParser(description='Generate methylation distribution plot')
    parser.add_argument('--bed-dir', required=True, help='Directory containing BED.gz files')
    parser.add_argument('--case-samples', required=True, help='Comma-separated list of case sample names')
    parser.add_argument('--control-samples', required=True, help='Comma-separated list of control sample names')
    parser.add_argument('--output-png', required=True, help='Output PNG file path')
    parser.add_argument('--output-pdf', required=True, help='Output PDF file path')
    args = parser.parse_args()

    case_samples = args.case_samples.split(",")
    control_samples = args.control_samples.split(",")

    # Sample methylation values from each sample
    meth_values = []
    sample_labels = []
    group_labels = []

    for sample in case_samples + control_samples:
        bed_file = f"{args.bed_dir}/{sample}.combined.bed.gz"
        if os.path.exists(bed_file):
            # Read a sample of methylation values
            values = []
            with gzip.open(bed_file, 'rt') as f:
                for i, line in enumerate(f):
                    if i >= 100000:  # Sample first 100k sites
                        break
                    fields = line.strip().split('\t')
                    if len(fields) >= 4:
                        try:
                            values.append(float(fields[3]))
                        except ValueError:
                            continue
            
            meth_values.extend(values)
            sample_labels.extend([sample] * len(values))
            group = 'Case' if sample in case_samples else 'Control'
            group_labels.extend([group] * len(values))

    if meth_values:
        df = pd.DataFrame({
            'Methylation': meth_values,
            'Sample': sample_labels,
            'Group': group_labels
        })
        
        # Create violin plot
        fig, ax = plt.subplots(figsize=(12, 6))
        palette = {'Case': '#e74c3c', 'Control': '#3498db'}
        sns.violinplot(data=df, x='Sample', y='Methylation', hue='Group',
                      palette=palette, ax=ax, inner='box')
        
        ax.set_xlabel('Sample')
        ax.set_ylabel('Methylation Level')
        ax.set_title('Methylation Distribution by Sample')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()
    else:
        # Create empty plot if no data
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, 'No methylation data available', ha='center', va='center', transform=ax.transAxes)
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    main()
