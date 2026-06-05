#!/usr/bin/env python3
"""
Convert DSS DMR TSV to the column format expected by annotate_dmrs.R.

DSS callDMR output columns: chr, start, end, length, nCG, meanMethy1,
    meanMethy2, diff.Methy, areaStat
Expected by annotate_dmrs.R: CHROM, START, END, DELTA
"""

import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser(
        description="Convert DSS DMR TSV for annotation pipeline"
    )
    parser.add_argument("--input-tsv", required=True, help="DSS DMR TSV")
    parser.add_argument("--output-tsv", required=True, help="Converted TSV")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    df = pd.read_csv(args.input_tsv, sep='\t')
    df = df.rename(columns={
        'chr': 'CHROM',
        'start': 'START',
        'end': 'END',
        'diff.Methy': 'DELTA',
    })
    df[['CHROM', 'START', 'END', 'DELTA']].to_csv(
        args.output_tsv, sep='\t', index=False
    )


if __name__ == "__main__":
    main()
