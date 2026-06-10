#!/usr/bin/env python3
"""
Filter DSS DMRs for IGV reporting.
Sorts by abs(diff.Methy) descending, takes top N, writes filtered TSV and BED.
"""

import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser(
        description="Filter DSS DMRs by absolute methylation difference for IGV"
    )
    parser.add_argument("--input-tsv", required=True, help="DSS dmr_results.tsv")
    parser.add_argument("--input-bed", required=True, help="DSS dmr_results.bed")
    parser.add_argument("--top-n", type=int, default=500, help="Number of top DMRs to keep")
    parser.add_argument("--output-tsv", required=True, help="Filtered TSV output")
    parser.add_argument("--output-bed", required=True, help="Filtered BED output")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    tsv = pd.read_csv(args.input_tsv, sep='\t')
    bed = pd.read_csv(args.input_bed, sep='\t', header=None,
                      names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

    diff_col = None
    for col in ['diff.Methy', 'diff_Methy', 'diff']:
        if col in tsv.columns:
            diff_col = col
            break

    if diff_col is None:
        print(f"Error: methylation difference column not found. Available: {list(tsv.columns)}")
        raise SystemExit(1)

    tsv['_abs_diff'] = tsv[diff_col].abs()
    tsv = tsv.sort_values('_abs_diff', ascending=False).head(args.top_n)
    tsv = tsv.drop(columns=['_abs_diff'])

    tsv.to_csv(args.output_tsv, sep='\t', index=False)

    region_ids = set()
    for _, row in tsv.iterrows():
        chrom = row.get('chr', row.get('CHROM', ''))
        start = row.get('start', row.get('START', 0))
        end = row.get('end', row.get('END', 0))
        region_ids.add(f"{chrom}:{int(start)}-{int(end)}")

    bed_filtered = bed[bed['name'].isin(region_ids)]
    bed_filtered.to_csv(args.output_bed, sep='\t', header=False, index=False)

    print(f"Filtered to top {len(tsv)} DMRs by abs(methylation difference)")


if __name__ == "__main__":
    main()
