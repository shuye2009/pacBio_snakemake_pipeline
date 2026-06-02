#!/usr/bin/env python3
"""
Extract differential ASM (allele-specific methylation) regions.

Identifies regions where ASM differs between case and control groups
based on avg_abs_meth_deltas difference. Outputs BED and TSV files
for IGV visualization.

Input: asm_cohort.profile.tsv from methbat build
Differential ASM: |case_avg_abs_meth_deltas - control_avg_abs_meth_deltas| >= min_delta
"""

import argparse
import pandas as pd
import numpy as np
import os


def main():
    parser = argparse.ArgumentParser(description='Extract differential ASM regions')
    parser.add_argument('--asm-comparison', required=True, help='Path to asm_cohort.profile.tsv')
    parser.add_argument('--min-delta', type=float, default=0.2, help='Minimum difference in avg_abs_meth_deltas between case and control')
    parser.add_argument('--pvalue-cutoff', type=float, default=0.05, help='P-value cutoff (not used currently)')
    parser.add_argument('--output-bed', required=True, help='Output BED file')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file for igv-reports')
    args = parser.parse_args()

    # Handle empty files
    if os.path.getsize(args.asm_comparison) == 0:
        print(f"Warning: Input file {args.asm_comparison} is empty")
        open(args.output_bed, 'w').close()
        with open(args.output_tsv, 'w') as f:
            f.write("CHROM\tSTART\tEND\tNAME\tCASE_DELTA\tCONTROL_DELTA\tDIFF_DELTA\n")
        return
    
    # Load cohort profile, skip comment lines
    df = pd.read_csv(args.asm_comparison, sep='\t', comment='#')
    print(f"Loaded {len(df)} rows from {args.asm_comparison}")
    
    # Separate case and control rows
    case_df = df[df['data_category'] == 'case'].copy()
    control_df = df[df['data_category'] == 'control'].copy()
    
    print(f"Found {len(case_df)} case rows and {len(control_df)} control rows")
    
    if len(case_df) == 0 or len(control_df) == 0:
        print("Warning: Missing case or control data")
        open(args.output_bed, 'w').close()
        with open(args.output_tsv, 'w') as f:
            f.write("CHROM\tSTART\tEND\tNAME\tCASE_DELTA\tCONTROL_DELTA\tDIFF_DELTA\n")
        return
    
    # Merge case and control on region coordinates
    case_df = case_df.rename(columns={'avg_abs_meth_deltas': 'case_avg_abs_meth_deltas'})
    control_df = control_df.rename(columns={'avg_abs_meth_deltas': 'control_avg_abs_meth_deltas'})
    
    merged = pd.merge(
        case_df[['chrom', 'start', 'end', 'cpg_label', 'case_avg_abs_meth_deltas']],
        control_df[['chrom', 'start', 'end', 'control_avg_abs_meth_deltas']],
        on=['chrom', 'start', 'end'],
        how='inner'
    )
    
    print(f"Merged {len(merged)} regions with both case and control data")
    
    # Calculate difference in avg_abs_meth_deltas between case and control
    merged['diff_delta'] = merged['case_avg_abs_meth_deltas'] - merged['control_avg_abs_meth_deltas']
    merged['abs_diff_delta'] = merged['diff_delta'].abs()
    
    # Filter for differential ASM regions
    sig_regions = merged[merged['abs_diff_delta'] >= args.min_delta].copy()
    
    print(f"Found {len(sig_regions)} differential ASM regions (|case - control| >= {args.min_delta})")
    
    if len(sig_regions) == 0:
        print("No differential ASM regions found")
        open(args.output_bed, 'w').close()
        with open(args.output_tsv, 'w') as f:
            f.write("CHROM\tSTART\tEND\tNAME\tCASE_DELTA\tCONTROL_DELTA\tDIFF_DELTA\n")
        return
    
    # Sort by absolute difference
    sig_regions = sig_regions.sort_values('abs_diff_delta', ascending=False)
    
    # Create BED file (no header)
    bed_df = pd.DataFrame()
    bed_df['chrom'] = sig_regions['chrom'].values
    bed_df['start'] = sig_regions['start'].values
    bed_df['end'] = sig_regions['end'].values
    bed_df['name'] = [f"{c}:{s}-{e}" for c, s, e in zip(bed_df['chrom'], bed_df['start'], bed_df['end'])]
    bed_df['score'] = (sig_regions['abs_diff_delta'] * 1000).astype(int).clip(0, 1000).values
    bed_df['strand'] = '.'
    bed_df.to_csv(args.output_bed, sep='\t', index=False, header=False)
    print(f"Saved {len(bed_df)} regions to BED file: {args.output_bed}")
    
    # Create TSV file (with header for igv-reports) - limit to top 100
    tsv_df = pd.DataFrame()
    tsv_df['CHROM'] = sig_regions['chrom'].values
    tsv_df['START'] = sig_regions['start'].values
    tsv_df['END'] = sig_regions['end'].values
    tsv_df['NAME'] = bed_df['name'].values
    tsv_df['CASE_DELTA'] = sig_regions['case_avg_abs_meth_deltas'].round(4).values
    tsv_df['CONTROL_DELTA'] = sig_regions['control_avg_abs_meth_deltas'].round(4).values
    tsv_df['DIFF_DELTA'] = sig_regions['diff_delta'].round(4).values
    
    # Limit to top 100 for igv-reports
    tsv_df = tsv_df.head(100)
    tsv_df.to_csv(args.output_tsv, sep='\t', index=False, header=True)
    print(f"Saved top {len(tsv_df)} regions to TSV file: {args.output_tsv}")


if __name__ == "__main__":
    main()
