#!/usr/bin/env python3
"""
Filter methbat significant DMRs for IGV reporting.
Sorts by ADJ_PVALUE ascending, takes top N, writes filtered TSV and BED.
"""

import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser(
        description="Filter methbat DMRs by adjusted p-value for IGV"
    )
    parser.add_argument("--input-tsv", required=True, help="significant_dmrs.tsv")
    parser.add_argument("--input-bed", required=True, help="significant_dmrs.bed")
    parser.add_argument("--top-n", type=int, default=500, help="Number of top DMRs to keep")
    parser.add_argument("--output-tsv", required=True, help="Filtered TSV output")
    parser.add_argument("--output-bed", required=True, help="Filtered BED output")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    tsv = pd.read_csv(args.input_tsv, sep='\t')
    bed = pd.read_csv(args.input_bed, sep='\t', header=None,
                      names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

    if 'ADJ_PVALUE' not in tsv.columns:
        print("Error: ADJ_PVALUE column not found in TSV")
        raise SystemExit(1)

    tsv = tsv.sort_values('ADJ_PVALUE', ascending=True).head(args.top_n)

    tsv.to_csv(args.output_tsv, sep='\t', index=False)

    region_ids = set(
        tsv['CHROM'] + ':' + tsv['START'].astype(str) + '-' + tsv['END'].astype(str)
    )
    bed_filtered = bed[bed['name'].isin(region_ids)]
    bed_filtered.to_csv(args.output_bed, sep='\t', header=False, index=False)

    print(f"Filtered to top {len(tsv)} DMRs by ADJ_PVALUE")


if __name__ == "__main__":
    main()
