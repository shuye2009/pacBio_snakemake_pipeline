#!/usr/bin/env python3
"""
Extract significant ASM (allele-specific methylation) regions.

Filters ASM haplotype comparison results by avg_abs_meth_deltas threshold
and outputs BED and TSV files for IGV visualization.

Input: asm_haplotype_comparison.tsv (filtered cohort profile with ASM regions)
"""

import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser(description='Extract significant ASM regions')
    parser.add_argument('--asm-comparison', required=True, help='Path to asm_haplotype_comparison.tsv')
    parser.add_argument('--min-delta', type=float, default=0.2, help='Minimum avg_abs_meth_deltas threshold')
    parser.add_argument('--output-bed', required=True, help='Output BED file')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file for igv-reports')
    args = parser.parse_args()

    # Handle empty files
    if os.path.getsize(args.asm_comparison) == 0:
        print(f"Warning: Input file {args.asm_comparison} is empty")
        open(args.output_bed, 'w').close()
        with open(args.output_tsv, 'w') as f:
            f.write("CHROM\tSTART\tEND\tNAME\tAVG_DELTA\tAVG_METHYL\n")
        return
    
    # Load ASM comparison data, skip comment lines
    df = pd.read_csv(args.asm_comparison, sep='\t', comment='#')
    print(f"Loaded {len(df)} rows from {args.asm_comparison}")
    
    if len(df) == 0:
        print("No data in input file")
        open(args.output_bed, 'w').close()
        with open(args.output_tsv, 'w') as f:
            f.write("CHROM\tSTART\tEND\tNAME\tAVG_DELTA\tAVG_METHYL\n")
        return
    
    # Filter by avg_abs_meth_deltas threshold
    if 'avg_abs_meth_deltas' in df.columns:
        sig_regions = df[df['avg_abs_meth_deltas'] >= args.min_delta].copy()
        print(f"Found {len(sig_regions)} regions with avg_abs_meth_deltas >= {args.min_delta}")
    else:
        print("Warning: avg_abs_meth_deltas column not found, using all regions")
        sig_regions = df.copy()
    
    if len(sig_regions) == 0:
        print("No significant ASM regions found")
        open(args.output_bed, 'w').close()
        with open(args.output_tsv, 'w') as f:
            f.write("CHROM\tSTART\tEND\tNAME\tAVG_DELTA\tAVG_METHYL\n")
        return
    
    # Sort by avg_abs_meth_deltas descending
    if 'avg_abs_meth_deltas' in sig_regions.columns:
        sig_regions = sig_regions.sort_values('avg_abs_meth_deltas', ascending=False)
    
    # Create BED file (no header)
    bed_df = pd.DataFrame()
    bed_df['chrom'] = sig_regions['chrom'].values
    bed_df['start'] = sig_regions['start'].values
    bed_df['end'] = sig_regions['end'].values
    bed_df['name'] = [f"{c}:{s}-{e}" for c, s, e in zip(bed_df['chrom'], bed_df['start'], bed_df['end'])]
    if 'avg_abs_meth_deltas' in sig_regions.columns:
        bed_df['score'] = (sig_regions['avg_abs_meth_deltas'] * 1000).astype(int).clip(0, 1000).values
    else:
        bed_df['score'] = 0
    bed_df['strand'] = '.'
    bed_df.to_csv(args.output_bed, sep='\t', index=False, header=False)
    print(f"Saved {len(bed_df)} regions to BED file: {args.output_bed}")
    
    # Create TSV file (with header for igv-reports) - limit to top 100
    tsv_df = pd.DataFrame()
    tsv_df['CHROM'] = sig_regions['chrom'].values
    tsv_df['START'] = sig_regions['start'].values
    tsv_df['END'] = sig_regions['end'].values
    tsv_df['NAME'] = bed_df['name'].values
    
    if 'avg_abs_meth_deltas' in sig_regions.columns:
        tsv_df['AVG_DELTA'] = sig_regions['avg_abs_meth_deltas'].round(4).values
    else:
        tsv_df['AVG_DELTA'] = 0
    
    if 'avg_combined_methyls' in sig_regions.columns:
        tsv_df['AVG_METHYL'] = sig_regions['avg_combined_methyls'].round(4).values
    else:
        tsv_df['AVG_METHYL'] = 0
    
    # Limit to top 100 for igv-reports
    tsv_df = tsv_df.head(100)
    tsv_df.to_csv(args.output_tsv, sep='\t', index=False, header=True)
    print(f"Saved top {len(tsv_df)} regions to TSV file: {args.output_tsv}")


if __name__ == "__main__":
    main()
