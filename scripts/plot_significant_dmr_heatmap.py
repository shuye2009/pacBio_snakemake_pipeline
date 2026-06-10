#!/usr/bin/env python
"""
Plot heatmap of significantly differentially methylated regions only.
Extracts methylation values from pb-cpg-tools combined.bed.gz files for DMR coordinates.
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


def load_dmr_regions(dmr_stats_tsv):
    """Load significant DMR regions from dmr_stats_with_pvalues.tsv file.
    
    Filters for regions where significant=True.
    Returns tuple of (regions list, filtered dataframe for BED output).
    """
    df = pd.read_csv(dmr_stats_tsv, sep='\t')
    
    # Filter for significant regions
    if 'significant' in df.columns:
        df = df[df['significant'] == True].copy()
    
    print(f"Found {len(df)} significant regions out of total")
    
    regions = []
    for _, row in df.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])
        name = f"{chrom}:{start}-{end}"
        regions.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'name': name
        })
    return regions, df


def extract_methylation_for_regions(bed_gz_file, regions):
    """Extract mean methylation for each region from a combined.bed.gz file.
    
    Uses tabix index (.tbi) for fast random access to specific genomic regions.
    pb-cpg-tools combined.bed.gz format:
    #chrom  begin   end     mod_score       type    cov     est_mod_count   est_unmod_count discretized_mod_score
    mod_score is in column index 3, values are 0-100 (need to convert to 0-1)
    """
    region_means = {}
    
    # Open tabix-indexed BED file
    try:
        tbx = pysam.TabixFile(bed_gz_file)
    except Exception as e:
        print(f"Warning: Could not open tabix file {bed_gz_file}: {e}")
        return {r['name']: np.nan for r in regions}
    
    try:
        for r in regions:
            chrom = r['chrom']
            start = r['start']
            end = r['end']
            name = r['name']
            
            meth_values = []
            try:
                # Query the specific region using tabix index
                for row in tbx.fetch(chrom, start, end):
                    # Skip header/comment lines
                    if row.startswith('#'):
                        continue
                    fields = row.split('\t')
                    # Column 3 is mod_score (0-100 scale)
                    if len(fields) >= 4:
                        try:
                            # Convert from 0-100 to 0-1 scale
                            meth = float(fields[3]) / 100.0
                            meth_values.append(meth)
                        except ValueError:
                            continue
            except ValueError:
                # Region not found in file (e.g., chromosome not present)
                pass
            
            if meth_values:
                region_means[name] = np.mean(meth_values)
            else:
                region_means[name] = np.nan
    finally:
        tbx.close()
    
    return region_means


def main():
    parser = argparse.ArgumentParser(description='Generate heatmap for significant DMRs')
    parser.add_argument('--dmr-stats', required=True, help='Path to dmr_stats_with_pvalues.tsv file')
    parser.add_argument('--bed-dir', required=True, help='Directory containing combined.bed.gz files')
    parser.add_argument('--case-samples', required=True, help='Comma-separated list of case sample names')
    parser.add_argument('--control-samples', required=True, help='Comma-separated list of control sample names')
    parser.add_argument('--output-png', required=True, help='Output PNG file path')
    parser.add_argument('--output-pdf', required=True, help='Output PDF file path')
    parser.add_argument('--output-bed', required=True, help='Output BED file with significant DMRs')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file with significant DMRs for igv-reports')
    args = parser.parse_args()

    case_samples = args.case_samples.split(",")
    control_samples = args.control_samples.split(",")

    # Load significant DMR regions
    regions, dmr_df = load_dmr_regions(args.dmr_stats)
    print(f"Loaded {len(regions)} significant DMR regions")
    
    # Save significant DMRs as BED file (standard format, no header)
    if regions and len(dmr_df) > 0:
        bed_df = pd.DataFrame()
        bed_df['chrom'] = dmr_df['chrom'].values
        bed_df['start'] = dmr_df['start'].values
        bed_df['end'] = dmr_df['end'].values
        bed_df['name'] = [f"{c}:{s}-{e}" for c, s, e in zip(bed_df['chrom'], bed_df['start'], bed_df['end'])]
        # Add score (use delta scaled to 0-1000)
        if 'delta' in dmr_df.columns:
            bed_df['score'] = (abs(dmr_df['delta']) * 1000).astype(int).clip(0, 1000)
        else:
            bed_df['score'] = 0
        bed_df['strand'] = '.'
        bed_df.to_csv(args.output_bed, sep='\t', index=False, header=False)
        print(f"Saved {len(bed_df)} significant DMRs to BED file: {args.output_bed}")
        
        # Save significant DMRs as TSV file (with header for igv-reports) - limit to top 100
        tsv_df = pd.DataFrame()
        tsv_df['CHROM'] = dmr_df['chrom'].values
        tsv_df['START'] = dmr_df['start'].values
        tsv_df['END'] = dmr_df['end'].values
        tsv_df['NAME'] = bed_df['name'].values
        if 'delta' in dmr_df.columns:
            tsv_df['DELTA'] = dmr_df['delta'].values
        else:
            tsv_df['DELTA'] = 0
        if 'adj_pvalue' in dmr_df.columns:
            tsv_df['ADJ_PVALUE'] = dmr_df['adj_pvalue'].values
        else:
            tsv_df['ADJ_PVALUE'] = np.nan
        
        tsv_df.to_csv(args.output_tsv, sep='\t', index=False, header=True)
        print(f"Saved {len(tsv_df)} significant DMRs to TSV file: {args.output_tsv}")
    else:
        # Create empty BED file
        open(args.output_bed, 'w').close()
        print(f"Created empty BED file (no significant DMRs): {args.output_bed}")
        # Create empty TSV file with header
        with open(args.output_tsv, 'w') as f:
            f.write("CHROM\tSTART\tEND\tNAME\tDELTA\tADJ_PVALUE\n")
        print(f"Created empty TSV file (no significant DMRs): {args.output_tsv}")

    if len(regions) == 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, 'No significant DMRs found in dmr_stats_with_pvalues.tsv', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()
        return

    # Extract methylation values for each sample
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

    # Create methylation matrix
    meth_matrix = pd.DataFrame(meth_data)
    
    # Remove regions with all NaN values
    meth_matrix = meth_matrix.dropna(how='all')
    
    # Fill remaining NaN with 0.5 (neutral)
    meth_matrix = meth_matrix.fillna(0.5)

    print(f"Created matrix with {len(meth_matrix)} regions x {len(meth_matrix.columns)} samples")

    # Truncate row labels to 50 characters
    meth_matrix.index = [str(idx)[:50] if len(str(idx)) > 50 else str(idx) for idx in meth_matrix.index]

    # Create sample annotations for column colors
    sample_colors = pd.Series(
        ['#e74c3c' if s in case_samples else '#3498db' for s in meth_matrix.columns],
        index=meth_matrix.columns,
        name='Group'
    )

    # Plot clustered heatmap
    n_regions = len(meth_matrix)
    fig_height = max(8, min(50, n_regions * 0.3))
    
    if len(meth_matrix) > 1 and len(meth_matrix.columns) > 1:
        # Use clustermap for hierarchical clustering
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
        
        # Rotate x-axis labels
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
        
        # Add title
        g.fig.suptitle(f'Significant DMRs (n={n_regions})', y=1.02)
        
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
        plt.title(f'Significant DMRs (n={n_regions})')
        plt.xlabel('Samples (Red=Case, Blue=Control)')
        plt.ylabel('DMR Regions')
        plt.xticks(rotation=45, ha='right')
        plt.savefig(args.output_png, dpi=300, bbox_inches='tight')
        plt.savefig(args.output_pdf, bbox_inches='tight')
        plt.close()
    
    print(f"Saved heatmap with {n_regions} significant DMR regions")


if __name__ == '__main__':
    main()
