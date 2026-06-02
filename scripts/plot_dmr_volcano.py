#!/usr/bin/env python
"""
Plot volcano plot of differentially methylated regions.
Uses -log10(adjusted p-value) on y-axis computed from z-scores.
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests


def main():
    parser = argparse.ArgumentParser(description='Generate DMR volcano plot')
    parser.add_argument('--stats', required=True, help='Path to signature stats TSV file')
    parser.add_argument('--min-delta', type=float, required=True, help='Minimum delta threshold')
    parser.add_argument('--pvalue-cutoff', type=float, default=0.05, help='Adjusted p-value cutoff for significance')
    parser.add_argument('--output-png', required=True, help='Output PNG file path')
    parser.add_argument('--output-pdf', required=True, help='Output PDF file path')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file with p-values')
    parser.add_argument('--output-zscore-dist-png', required=False, help='Output PNG for z-score distribution')
    parser.add_argument('--output-zscore-dist-pdf', required=False, help='Output PDF for z-score distribution')
    args = parser.parse_args()

    # Significance cutoff for adjusted p-value
    PVALUE_CUTOFF = args.pvalue_cutoff

    # Load stats
    try:
        stats_df = pd.read_csv(args.stats, sep='\t', comment='#')
    except Exception:
        stats_df = pd.read_csv(args.stats, sep='\t', comment='#', on_bad_lines='skip')

    # Print columns for debugging
    print(f"Available columns: {list(stats_df.columns)}")

    # Create volcano plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # MethBat signature_stats.tsv columns:
    # chrom, start, end, z_score, r_value, t_value, num_samples1, num_samples2, mean1, mean2, median1, median2
    # mean1 = case mean, mean2 = control mean (based on --baseline-category control --compare-category case)
    
    # Find z_score column
    zscore_col = None
    for col in stats_df.columns:
        col_lower = col.lower()
        if col_lower == 'z_score' or col_lower == 'zscore' or col_lower == 't_value':
            zscore_col = col
            break
    
    # Calculate delta from mean1 - mean2 if not present
    delta_col = None
    for col in stats_df.columns:
        col_lower = col.lower()
        if col_lower == 'delta' or col_lower == 'diff':
            delta_col = col
            break
    
    # If no delta column, calculate from mean1 - mean2
    if delta_col is None and 'mean1' in stats_df.columns and 'mean2' in stats_df.columns:
        stats_df['delta'] = stats_df['mean1'] - stats_df['mean2']
        delta_col = 'delta'
        print("Calculated delta from mean1 - mean2")

    print(f"Detected delta_col: {delta_col}, zscore_col: {zscore_col}")

    if delta_col and zscore_col:
        # Compute p-values from z-scores (two-tailed test, normal distribution)
        zscores = stats_df[zscore_col].values
        
        # Remove NaN values and outliers for distribution fitting
        zscores_clean = zscores[~np.isnan(zscores) & (np.abs(zscores) < 10)]
        
        # Fit normal distribution to the data
        mu, std = scipy_stats.norm.fit(zscores_clean)
        
        # Plot z-score distribution with fitted normal curve
        if args.output_zscore_dist_png or args.output_zscore_dist_pdf:
            fig_dist, ax_dist = plt.subplots(figsize=(10, 6))
            
            # Plot histogram as density
            ax_dist.hist(zscores_clean, bins=100, density=True, 
                         alpha=0.7, color='steelblue', edgecolor='white',
                         label='Observed z-scores')
            
            # Generate x values for the fitted curve
            x = np.linspace(zscores_clean.min(), zscores_clean.max(), 1000)
            fitted_pdf = scipy_stats.norm.pdf(x, mu, std)
            
            # Plot fitted normal distribution
            ax_dist.plot(x, fitted_pdf, 'r-', linewidth=2, 
                        label=f'Fitted Normal (μ={mu:.3f}, σ={std:.3f})')
            
            # Plot standard normal for comparison
            standard_pdf = scipy_stats.norm.pdf(x, 0, 1)
            ax_dist.plot(x, standard_pdf, 'g--', linewidth=2, 
                        label='Standard Normal (μ=0, σ=1)')
            
            # Perform Shapiro-Wilk test (on a sample if too large)
            if len(zscores_clean) > 5000:
                sample_zscores = np.random.choice(zscores_clean, 5000, replace=False)
            else:
                sample_zscores = zscores_clean
            shapiro_stat, shapiro_p = scipy_stats.shapiro(sample_zscores)
            
            # Perform Kolmogorov-Smirnov test against standard normal
            ks_stat, ks_p = scipy_stats.kstest(zscores_clean, 'norm', args=(mu, std))
            
            # Add statistics text box
            stats_text = (f'N = {len(zscores_clean):,}\n'
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
            print(f"Fitted z-score distribution: μ={mu:.4f}, σ={std:.4f}")
        
        # Compute p-values using fitted normal distribution parameters
        pvalues = 2 * scipy_stats.norm.sf(zscores, loc=mu, scale=std)  # two-tailed p-value
        
        # Handle edge cases (very small p-values)
        pvalues = np.clip(pvalues, 1e-300, 1.0)
        
        # Add raw p-value to dataframe
        stats_df['pvalue'] = pvalues
        
        # Adjust p-values for multiple testing (Benjamini-Hochberg FDR)
        _, adj_pvalues, _, _ = multipletests(pvalues, method='fdr_bh')
        
        # Compute -log10(adjusted p-value)
        neg_log10_adj_pval = -np.log10(adj_pvalues)
        stats_df['neg_log10_adj_pval'] = neg_log10_adj_pval
        stats_df['adj_pvalue'] = adj_pvalues
        
        # Mark significant regions
        stats_df['significant'] = (abs(stats_df[delta_col]) >= args.min_delta) & (stats_df['adj_pvalue'] <= PVALUE_CUTOFF)
        
        # Save stats with p-values to TSV
        stats_df.to_csv(args.output_tsv, sep='\t', index=False)
        print(f"Saved stats with p-values to {args.output_tsv}")
        
        # Significance threshold line
        neg_log10_cutoff = -np.log10(PVALUE_CUTOFF)
        
        # Determine significant points
        significant = (abs(stats_df[delta_col]) >= args.min_delta) & (stats_df['adj_pvalue'] <= PVALUE_CUTOFF)
        
        print(f"Total regions: {len(stats_df)}")
        print(f"Significant regions (|delta| >= {args.min_delta}, adj_pval <= {PVALUE_CUTOFF}): {significant.sum()}")
        
        # Plot non-significant points
        ax.scatter(stats_df.loc[~significant, delta_col], 
                  stats_df.loc[~significant, 'neg_log10_adj_pval'],
                  c='gray', alpha=0.5, s=20, label='Not significant')
        
        # Plot significant points
        hyper = significant & (stats_df[delta_col] > 0)
        hypo = significant & (stats_df[delta_col] < 0)
        
        ax.scatter(stats_df.loc[hyper, delta_col], 
                  stats_df.loc[hyper, 'neg_log10_adj_pval'],
                  c='#e74c3c', alpha=0.7, s=30, label=f'Hypermethylated (n={hyper.sum()})')
        ax.scatter(stats_df.loc[hypo, delta_col], 
                  stats_df.loc[hypo, 'neg_log10_adj_pval'],
                  c='#3498db', alpha=0.7, s=30, label=f'Hypomethylated (n={hypo.sum()})')
        
        # Add threshold lines
        ax.axhline(y=neg_log10_cutoff, color='black', linestyle='--', alpha=0.5, 
                   label=f'adj. p-value = {PVALUE_CUTOFF}')
        ax.axvline(x=args.min_delta, color='black', linestyle='--', alpha=0.3)
        ax.axvline(x=-args.min_delta, color='black', linestyle='--', alpha=0.3)
        
        ax.set_xlabel('Methylation Difference (Case - Control)')
        ax.set_ylabel('-log₁₀(adjusted p-value)')
    else:
        # Fallback: plot first two numeric columns
        numeric_cols = stats_df.select_dtypes(include=['float64', 'int64']).columns.tolist()
        if len(numeric_cols) >= 2:
            ax.scatter(stats_df[numeric_cols[0]], stats_df[numeric_cols[1]], c='gray', alpha=0.5, s=20)
            ax.set_xlabel(numeric_cols[0])
            ax.set_ylabel(numeric_cols[1])
        else:
            ax.text(0.5, 0.5, 'Insufficient data for volcano plot', ha='center', va='center', transform=ax.transAxes)

    ax.set_title('Differentially Methylated Regions')
    ax.legend(loc='upper right')
    plt.tight_layout()

    plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
    plt.savefig(args.output_pdf, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    main()
