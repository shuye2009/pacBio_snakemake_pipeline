#!/usr/bin/env python
"""
Plot PCA of samples based on methylation profiles.
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import os


def main():
    parser = argparse.ArgumentParser(description='Generate sample PCA plot')
    parser.add_argument('--profile-dir', required=True, help='Directory containing profile TSV files')
    parser.add_argument('--case-samples', required=True, help='Comma-separated list of case sample names')
    parser.add_argument('--control-samples', required=True, help='Comma-separated list of control sample names')
    parser.add_argument('--output-png', required=True, help='Output PNG file path')
    parser.add_argument('--output-pdf', required=True, help='Output PDF file path')
    args = parser.parse_args()

    case_samples = args.case_samples.split(",")
    control_samples = args.control_samples.split(",")

    # Load methylation profiles
    meth_data = {}
    for sample in case_samples + control_samples:
        profile_file = f"{args.profile_dir}/{sample}.region.profile.tsv"
        if os.path.exists(profile_file):
            try:
                df = pd.read_csv(profile_file, sep='\t', comment='#')
            except Exception:
                df = pd.read_csv(profile_file, sep='\t', comment='#', on_bad_lines='skip')
            
            # Find the index column (region identifier)
            index_col = None
            for col in ['region', 'region_id', 'name', 'cpg_label']:
                if col in df.columns:
                    index_col = col
                    break
            if index_col is None and len(df.columns) > 0:
                index_col = df.columns[0]
            
            # Find the methylation value column
            if index_col and len(df.columns) > 1:
                if 'mean_combined_methyl' in df.columns:
                    meth_data[sample] = df.set_index(index_col)['mean_combined_methyl']
                elif 'methylation' in df.columns:
                    meth_data[sample] = df.set_index(index_col)['methylation']
                else:
                    numeric_cols = df.select_dtypes(include=['float64', 'int64']).columns
                    if len(numeric_cols) > 0:
                        meth_data[sample] = df.set_index(index_col)[numeric_cols[0]]

    if len(meth_data) >= 2:
        # Create matrix and handle missing values
        meth_matrix = pd.DataFrame(meth_data).T.fillna(0.5)
        
        # Perform PCA
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(meth_matrix)
        pca = PCA(n_components=min(2, len(meth_matrix)))
        pca_result = pca.fit_transform(scaled_data)
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        for i, sample in enumerate(meth_matrix.index):
            color = '#e74c3c' if sample in case_samples else '#3498db'
            marker = 'o' if sample in case_samples else 's'
            ax.scatter(pca_result[i, 0], pca_result[i, 1] if pca_result.shape[1] > 1 else 0,
                      c=color, marker=marker, s=100, alpha=0.7)
            ax.annotate(sample, (pca_result[i, 0], pca_result[i, 1] if pca_result.shape[1] > 1 else 0),
                       fontsize=8, alpha=0.8)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        if pca_result.shape[1] > 1:
            ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        ax.set_title('PCA of Methylation Profiles\n(Red circles=Case, Blue squares=Control)')
        plt.tight_layout()
        
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()
    else:
        # Create empty plot if insufficient data
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, 'Insufficient data for PCA (need at least 2 samples)', ha='center', va='center', transform=ax.transAxes)
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    main()
