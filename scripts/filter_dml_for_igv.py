#!/usr/bin/env python3
"""
Filter DSS DML results to the top N most significant loci for IGV reporting.

Reads the full DML TSV, sorts by p-value (ascending), takes the top N rows,
and writes filtered TSV and BED files.
"""

import argparse
import pandas as pd
import numpy as np
import os


def main():
    parser = argparse.ArgumentParser(
        description="Filter DML results to top N by p-value for IGV"
    )
    parser.add_argument("--dml-tsv", required=True, help="Full DML TSV from DSS")
    parser.add_argument("--top-n", type=int, default=500,
                        help="Number of top DMLs to keep")
    parser.add_argument("--output-tsv", required=True, help="Filtered TSV output")
    parser.add_argument("--output-bed", required=True, help="Filtered BED output")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    dml = pd.read_csv(args.dml_tsv, sep='\t')

    if 'pval' not in dml.columns:
        print(f"Error: 'pval' column not found. Available columns: {list(dml.columns)}")
        raise SystemExit(1)

    dml = dml.sort_values('pval', ascending=True).head(args.top_n)

    dml.to_csv(args.output_tsv, sep='\t', index=False)

    bed = dml[['chr', 'pos']].copy()
    bed['start'] = bed['pos'] - 1
    bed['end'] = bed['pos']
    bed['name'] = bed['chr'] + ':' + bed['pos'].astype(str)
    bed['score'] = np.clip((-np.log10(dml['pval'].values) * 100).astype(int), 0, 1000)
    bed['strand'] = '.'
    bed = bed[['chr', 'start', 'end', 'name', 'score', 'strand']]
    bed.to_csv(args.output_bed, sep='\t', header=False, index=False)


if __name__ == "__main__":
    main()
