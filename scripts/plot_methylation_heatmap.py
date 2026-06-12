#!/usr/bin/env python
"""
Plot methylation heatmap across samples at significant regions from cohort comparison.
Filters regions based on min_delta and adjusted p-value thresholds.
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests


def main():
    parser = argparse.ArgumentParser(description='Generate methylation heatmap for significant regions')
    parser.add_argument('--cohort-comparison', required=True, help='Path to cohort_comparison.tsv file')
    parser.add_argument('--profile-dir', required=True, help='Directory containing profile TSV files')
    parser.add_argument('--case-samples', required=True, help='Comma-separated list of case sample names')
    parser.add_argument('--control-samples', required=True, help='Comma-separated list of control sample names')
    parser.add_argument('--min-delta', type=float, required=True, help='Minimum delta threshold')
    parser.add_argument('--pvalue-cutoff', type=float, default=0.05, help='Adjusted p-value cutoff for significance')
    parser.add_argument('--output-png', required=True, help='Output PNG file path')
    parser.add_argument('--output-pdf', required=True, help='Output PDF file path')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file with p-values')
    parser.add_argument('--output-bed', required=True, help='Output BED file with significant regions')
    parser.add_argument('--output-igv-tsv', required=True, help='Output TSV file with significant regions for igv-reports')
    parser.add_argument('--output-zscore-dist-png', required=False, help='Output PNG for z-score distribution')
    parser.add_argument('--output-zscore-dist-pdf', required=False, help='Output PDF for z-score distribution')
    parser.add_argument('--top-n', type=int, default=500, help='Top N regions for heatmap')
    parser.add_argument('--igv-top-n', type=int, default=100, help='Top N regions for IGV TSV')
    args = parser.parse_args()

    case_samples = args.case_samples.split(",")
    control_samples = args.control_samples.split(",")

    # Load cohort comparison to get significant regions
    try:
        comparison = pd.read_csv(args.cohort_comparison, sep='\t', comment='#')
    except Exception:
        comparison = pd.read_csv(args.cohort_comparison, sep='\t', comment='#', on_bad_lines='skip')

    print(f"Loaded {len(comparison)} regions from cohort_comparison")
    print(f"Columns: {list(comparison.columns)}")

    # Find delta and zscore columns from MethBat compare output
    # Prefer combined methylation columns over phased/haplotype columns
    delta_col = None
    zscore_col = None
    
    # Check for MethBat compare output columns first
    if 'delta_avg_combined_methyls' in comparison.columns:
        delta_col = 'delta_avg_combined_methyls'
    if 'zscore_avg_combined_methyls' in comparison.columns:
        zscore_col = 'zscore_avg_combined_methyls'
    
    # Fallback to other column names
    if delta_col is None or zscore_col is None:
        for col in comparison.columns:
            col_lower = col.lower()
            if delta_col is None:
                if col_lower == 'delta' or col_lower == 'diff':
                    delta_col = col
                elif 'delta' in col_lower and 'abs' not in col_lower:
                    delta_col = col
            if zscore_col is None:
                if col_lower == 'z_score' or col_lower == 'zscore':
                    zscore_col = col

    # Calculate delta from means if still not present
    if delta_col is None and 'case_mean' in comparison.columns and 'control_mean' in comparison.columns:
        comparison['delta'] = comparison['case_mean'] - comparison['control_mean']
        delta_col = 'delta'
    elif delta_col is None and 'mean1' in comparison.columns and 'mean2' in comparison.columns:
        comparison['delta'] = comparison['mean1'] - comparison['mean2']
        delta_col = 'delta'

    print(f"Using delta_col: {delta_col}, zscore_col: {zscore_col}")

    # Compute p-values from z-scores and apply FDR correction
    if zscore_col:
        zscores = comparison[zscore_col].values
        
        # Handle NaN values in z-scores
        valid_mask = ~np.isnan(zscores)
        print(f"Valid z-scores (non-NaN): {valid_mask.sum()} out of {len(zscores)}")
        
        # Initialize p-value arrays with NaN
        pvalues = np.full(len(zscores), np.nan)
        adj_pvalues = np.full(len(zscores), np.nan)
        
        if valid_mask.sum() > 0:
            # Compute two-tailed p-values from z-scores (normal distribution)
            valid_zscores = zscores[valid_mask]
            filtered_zscores = valid_zscores[np.abs(valid_zscores) < 10]
            # Fit normal distribution to the data
            mu, std = scipy_stats.norm.fit(filtered_zscores)
                       
            # Plot z-score distribution with fitted normal curve
            if args.output_zscore_dist_png or args.output_zscore_dist_pdf:
                fig_dist, ax_dist = plt.subplots(figsize=(10, 6))
                
                # Plot histogram as density
                ax_dist.hist(filtered_zscores, bins=100, density=True, 
                             alpha=0.7, color='steelblue', edgecolor='white',
                             label='Observed z-scores')
                 
                # Generate x values for the fitted curve
                x = np.linspace(filtered_zscores.min(), filtered_zscores.max(), 1000)
                fitted_pdf = scipy_stats.norm.pdf(x, mu, std)
                
                # Plot fitted normal distribution
                ax_dist.plot(x, fitted_pdf, 'r-', linewidth=2, 
                            label=f'Fitted Normal (μ={mu:.3f}, σ={std:.3f})')
                
                # Plot standard normal for comparison
                standard_pdf = scipy_stats.norm.pdf(x, 0, 1)
                ax_dist.plot(x, standard_pdf, 'g--', linewidth=2, 
                            label='Standard Normal (μ=0, σ=1)')
                
                # Perform Shapiro-Wilk test (on a sample if too large)
                if len(filtered_zscores) > 5000:
                    sample_zscores = np.random.choice(filtered_zscores, 5000, replace=False)
                else:
                    sample_zscores = filtered_zscores
                shapiro_stat, shapiro_p = scipy_stats.shapiro(sample_zscores)
                
                # Perform Kolmogorov-Smirnov test against fitted normal
                ks_stat, ks_p = scipy_stats.kstest(filtered_zscores, 'norm', args=(mu, std))
                
                # Add statistics text box
                stats_text = (f'N = {len(filtered_zscores):,}\n'
                             f'Fitted μ = {mu:.4f}\n'
                             f'Fitted σ = {std:.4f}\n'
                             f'Shapiro-Wilk p = {shapiro_p:.2e}\n'
                             f'K-S test p = {ks_p:.2e}')
                ax_dist.text(0.02, 0.98, stats_text, transform=ax_dist.transAxes, 
                            fontsize=10, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                
                ax_dist.set_xlabel('Z-score', fontsize=12)
                ax_dist.set_ylabel('Density', fontsize=12)
                ax_dist.set_title('Z-score Distribution with Fitted Normal Curve', fontsize=14)
                ax_dist.legend(loc='upper right')
                ax_dist.grid(True, alpha=0.3)
                
                plt.tight_layout()
                
                if args.output_zscore_dist_png:
                    fig_dist.savefig(args.output_zscore_dist_png, dpi=300, bbox_inches='tight')
                    print(f"Saved z-score distribution plot to {args.output_zscore_dist_png}")
                if args.output_zscore_dist_pdf:
                    fig_dist.savefig(args.output_zscore_dist_pdf, bbox_inches='tight')
                    print(f"Saved z-score distribution plot to {args.output_zscore_dist_pdf}")
                plt.close(fig_dist)
                
                print(f"Z-score distribution: μ={mu:.4f}, σ={std:.4f}")
                print(f"Shapiro-Wilk test: statistic={shapiro_stat:.4f}, p-value={shapiro_p:.2e}")
                print(f"K-S test against fitted normal: statistic={ks_stat:.4f}, p-value={ks_p:.2e}")
            else:
                # Fit normal distribution even if not plotting
                mu, std = scipy_stats.norm.fit(filtered_zscores)
                print(f"Fitted z-score distribution: μ={mu:.4f}, σ={std:.4f}") 
            
            # Compute p-values using fitted normal distribution parameters
            valid_pvalues = 2 * scipy_stats.norm.sf(abs(valid_zscores), loc=0, scale=std)
            # Handle edge cases
            valid_pvalues = np.clip(valid_pvalues, 1e-300, 1.0)
            pvalues[valid_mask] = valid_pvalues
            
            # Adjust p-values for multiple testing (Benjamini-Hochberg FDR)
            _, valid_adj_pvalues, _, _ = multipletests(valid_pvalues, method='fdr_bh')
            adj_pvalues[valid_mask] = valid_adj_pvalues
        
        comparison['pvalue'] = pvalues
        comparison['adj_pvalue'] = adj_pvalues
        
        print(f"Computed p-values for {valid_mask.sum()} regions with valid z-scores")
    else:
        print("Warning: No z-score column found, cannot compute p-values")
        comparison['pvalue'] = np.nan
        comparison['adj_pvalue'] = np.nan

    # Filter for significant regions using p-value cutoff (fall back to delta-only if no valid p-values)
    if delta_col and 'adj_pvalue' in comparison.columns:
        has_valid_pvalues = comparison['adj_pvalue'].notna().any()
        if has_valid_pvalues:
            comparison['significant'] = (abs(comparison[delta_col]) >= args.min_delta) & (comparison['adj_pvalue'] <= args.pvalue_cutoff)
        else:
            print("Warning: No valid p-values, falling back to delta-only significance")
            comparison['significant'] = abs(comparison[delta_col]) >= args.min_delta
        sig_regions = comparison[comparison['significant']].copy()
        print(f"Found {len(sig_regions)} significant regions (|delta| >= {args.min_delta}, adj_pvalue <= {args.pvalue_cutoff})")
    else:
        print("Warning: Could not filter by significance, using all regions")
        comparison['significant'] = True
        sig_regions = comparison
    
    # Save cohort comparison with p-values to TSV
    comparison.to_csv(args.output_tsv, sep='\t', index=False)
    print(f"Saved cohort comparison with p-values to {args.output_tsv}")
    
    # Save significant regions as BED file (standard format, no header) and TSV (with header for igv-reports)
    if 'chrom' in sig_regions.columns and 'start' in sig_regions.columns and 'end' in sig_regions.columns:
        # BED file (no header)
        bed_df = pd.DataFrame()
        bed_df['chrom'] = sig_regions['chrom'].values
        bed_df['start'] = sig_regions['start'].values
        bed_df['end'] = sig_regions['end'].values
        # Add name column (region label if available, otherwise coordinates)
        if 'cpg_label' in sig_regions.columns and sig_regions['cpg_label'].notna().any():
            bed_df['name'] = sig_regions['cpg_label'].values
        else:
            bed_df['name'] = [f"{c}:{s}-{e}" for c, s, e in zip(bed_df['chrom'], bed_df['start'], bed_df['end'])]
        # Add score (use delta scaled to 0-1000)
        if delta_col:
            bed_df['score'] = (abs(sig_regions[delta_col]) * 1000).astype(int).clip(0, 1000)
        else:
            bed_df['score'] = 0
        bed_df['strand'] = '.'
        bed_df.to_csv(args.output_bed, sep='\t', index=False, header=False)
        print(f"Saved {len(bed_df)} significant regions to BED file: {args.output_bed}")
        
        # TSV file (with header for igv-reports) - limit to top 100 by significance
        tsv_df = pd.DataFrame()
        tsv_df['CHROM'] = sig_regions['chrom'].values
        tsv_df['START'] = sig_regions['start'].values
        tsv_df['END'] = sig_regions['end'].values
        tsv_df['NAME'] = bed_df['name'].values
        if delta_col:
            tsv_df['DELTA'] = sig_regions[delta_col].values
        else:
            tsv_df['DELTA'] = 0
        if 'adj_pvalue' in sig_regions.columns:
            tsv_df['ADJ_PVALUE'] = sig_regions['adj_pvalue'].values
        else:
            tsv_df['ADJ_PVALUE'] = np.nan
        
        # Sort by significance and limit to top 100 for igv-reports
        tsv_df = tsv_df.sort_values('ADJ_PVALUE', ascending=True).head(args.igv_top_n)
        tsv_df.to_csv(args.output_igv_tsv, sep='\t', index=False, header=True)
        print(f"Saved top {len(tsv_df)} significant regions to TSV file: {args.output_igv_tsv}")

    # Create region identifiers from chrom:start-end for matching
    # Also create a mapping from coordinates to region labels for display
    coord_to_label = {}
    if 'chrom' in sig_regions.columns and 'start' in sig_regions.columns and 'end' in sig_regions.columns:
        sig_regions['region_key'] = sig_regions['chrom'] + ':' + sig_regions['start'].astype(str) + '-' + sig_regions['end'].astype(str)
    else:
        print("Warning: Could not find chrom/start/end columns in cohort_comparison")
        sig_regions['region_key'] = sig_regions.index.astype(str)
    
    significant_region_keys = set(sig_regions['region_key'].tolist())
    print(f"Significant region keys: {len(significant_region_keys)}")

    # Load methylation profiles for each sample (filter to significant regions to save memory)
    meth_data = {}
    coord_to_label = {}  # Map coordinates to region labels
    
    for sample in case_samples + control_samples:
        profile_file = f"{args.profile_dir}/{sample}.region.profile.tsv"
        if os.path.exists(profile_file):
            try:
                df = pd.read_csv(profile_file, sep='\t', comment='#')
            except Exception:
                df = pd.read_csv(profile_file, sep='\t', comment='#', on_bad_lines='skip')
            
            # Create region key from chrom:start-end for matching
            if 'chrom' in df.columns and 'start' in df.columns and 'end' in df.columns:
                df['region_key'] = df['chrom'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
            elif 'cpg_label' in df.columns:
                df['region_key'] = df['cpg_label']
            else:
                df['region_key'] = df.iloc[:, 0].astype(str)
            
            # Filter to significant regions only to save memory
            if significant_region_keys:
                df = df[df['region_key'].isin(significant_region_keys)]
            
            # Build coordinate to label mapping (use cpg_label if available)
            if 'cpg_label' in df.columns and not coord_to_label:
                for _, row in df.iterrows():
                    coord_to_label[row['region_key']] = row['cpg_label']
            
            # Find the methylation value column
            meth_col = None
            if 'mean_combined_methyl' in df.columns:
                meth_col = 'mean_combined_methyl'
            elif 'methylation' in df.columns:
                meth_col = 'methylation'
            else:
                numeric_cols = df.select_dtypes(include=['float64', 'int64']).columns.tolist()
                # Exclude coordinate columns
                numeric_cols = [c for c in numeric_cols if c not in ['start', 'end']]
                if numeric_cols:
                    meth_col = numeric_cols[0]
                    print(f"Warning: Using column {meth_col} for sample {sample}")
            
            if meth_col:
                meth_data[sample] = df.set_index('region_key')[meth_col]

    if meth_data:
        # Create methylation matrix
        meth_matrix = pd.DataFrame(meth_data)
        
        # Filter to only significant regions
        if significant_region_keys:
            matching_indices = [idx for idx in meth_matrix.index if idx in significant_region_keys]
            if matching_indices:
                meth_matrix = meth_matrix.loc[matching_indices]
                print(f"Filtered to {len(meth_matrix)} significant regions")
            else:
                print("Warning: No matching regions found between profiles and cohort_comparison")
                meth_matrix = meth_matrix.iloc[:0]
        else:
            print("Warning: No significant regions, skipping heatmap")
            meth_matrix = meth_matrix.iloc[:0]
        
        # Sort by significance and limit to top_n for heatmap
        if 'region_key' in sig_regions.columns and 'adj_pvalue' in sig_regions.columns:
            sig_order = sig_regions.set_index('region_key')['adj_pvalue'].sort_values()
            common_keys = [k for k in sig_order.index if k in meth_matrix.index]
            if common_keys:
                meth_matrix = meth_matrix.loc[common_keys]
                meth_matrix = meth_matrix.iloc[:args.top_n]
                print(f"Limited heatmap to top {len(meth_matrix)} regions by significance")
        
        # Convert coordinate-based index to region labels for y-axis display
        if coord_to_label:
            new_index = [coord_to_label.get(idx, idx) for idx in meth_matrix.index]
            meth_matrix.index = new_index
            print(f"Converted {sum(1 for idx in meth_matrix.index if idx in coord_to_label.values())} region labels")
        
        # Truncate row labels (y-axis) to 50 characters
        meth_matrix.index = [str(idx)[:50] if len(str(idx)) > 50 else str(idx) for idx in meth_matrix.index]
        
        # Create sample annotations for column colors
        sample_colors = pd.Series(
            ['#e74c3c' if s in case_samples else '#3498db' for s in meth_matrix.columns],
            index=meth_matrix.columns,
            name='Group'
        )
        
        # Plot clustered heatmap
        n_regions = len(meth_matrix)
        
        if n_regions == 0:
            print("No significant regions to plot, creating empty heatmap")
            fig, ax = plt.subplots(figsize=(10, 8))
            ax.text(0.5, 0.5, 'No significant regions found', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Significant Regions (n=0, |delta|>={args.min_delta}, adj_pval<={args.pvalue_cutoff})')
            plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
            plt.savefig(args.output_pdf, bbox_inches='tight')
            plt.close()
            print(f"Saved empty heatmap (0 significant regions)")
            return
        
        fig_height = max(8, min(50, n_regions * 0.3))
        
        # Drop rows/columns with all NaN values for clustering
        meth_matrix_clean = meth_matrix.dropna(how='all', axis=0).dropna(how='all', axis=1)
        # Fill remaining NaN with row mean for clustering
        meth_matrix_filled = meth_matrix_clean.apply(lambda x: x.fillna(x.mean()), axis=1)
        
        if len(meth_matrix_filled) > 1 and len(meth_matrix_filled.columns) > 1:
            # Use clustermap for hierarchical clustering
            g = sns.clustermap(
                meth_matrix_filled,
                cmap='RdYlBu_r',
                center=0.5,
                vmin=0, vmax=1,
                col_colors=sample_colors[meth_matrix_filled.columns],
                row_cluster=True,
                col_cluster=True,
                figsize=(max(14, len(meth_matrix_filled.columns) * 1.5 + 2), fig_height + 2),
                dendrogram_ratio=(0.15, 0.2),
                cbar_kws={'label': 'Methylation Level'},
                xticklabels=True,
                yticklabels=True,
                tree_kws={'linewidths': 1.5},
            )
            
            # Rotate x-axis labels
            plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
            
            # Add title
            g.fig.suptitle(f'Significant Regions (n={n_regions}, |delta|>={args.min_delta}, adj_pval<={args.pvalue_cutoff})', y=1.02)
            
            # Add legend for sample groups
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor='#e74c3c', label='Case'),
                             Patch(facecolor='#3498db', label='Control')]
            g.ax_heatmap.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1.15))
            
            g.savefig(args.output_png, dpi=300, bbox_inches='tight')
            g.savefig(args.output_pdf, bbox_inches='tight')
            plt.close()
        else:
            # Fallback to regular heatmap if not enough data for clustering
            fig, ax = plt.subplots(figsize=(max(12, len(meth_matrix.columns) * 1.5), fig_height))
            sns.heatmap(meth_matrix, cmap='RdYlBu_r', center=0.5, 
                       vmin=0, vmax=1, ax=ax,
                       xticklabels=True, yticklabels=True,
                       cbar_kws={'label': 'Methylation Level'})
            plt.title(f'Significant Regions (n={n_regions}, |delta|>={args.min_delta}, adj_pval<={args.pvalue_cutoff})')
            plt.xlabel('Samples (Red=Case, Blue=Control)')
            plt.ylabel('Regions')
            plt.xticks(rotation=45, ha='right')
            plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
            plt.savefig(args.output_pdf, bbox_inches='tight')
            plt.close()
        
        print(f"Saved heatmap with {n_regions} significant regions")
    else:
        # Create empty plot if no data
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, 'No methylation data available', ha='center', va='center', transform=ax.transAxes)
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    main()
