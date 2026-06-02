#!/usr/bin/env python3
"""
Average CpG mod_score per genomic bin for a single sample.

Reads a pb-CpG-tools combined BED file, tiles each chromosome into
fixed-size bins using the FASTA index, and outputs mean mod_score per bin.
"""

import argparse
import pandas as pd
import numpy as np
import os


def main():
    parser = argparse.ArgumentParser(
        description="Average mod_score per genomic bin for a single sample"
    )
    parser.add_argument("--bed", required=True, help="Path to combined.bed.gz")
    parser.add_argument("--fai", required=True, help="Path to genome FASTA index (.fai)")
    parser.add_argument("--bin-size", type=int, default=50000, help="Bin size in bp (default: 50000)")
    parser.add_argument("--output-tsv", required=True, help="Output TSV path")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    # Load chromosome sizes from FASTA index
    fai = pd.read_csv(args.fai, sep='\t', header=None,
                      names=["chrom", "length", "offset", "linebases", "linewidth"])

    # pb-CpG-tools combined BED columns (header line starts with #chrom):
    # chrom, begin, end, mod_score, type, cov, est_mod_count, est_unmod_count, discretized_mod_score
    df = pd.read_csv(args.bed, sep='\t', comment='#', header=None,
                     names=["chrom", "begin", "end", "mod_score", "type",
                            "cov", "est_mod_count", "est_unmod_count", "discretized_mod_score"],
                     compression='gzip')

    # Assign each CpG to a bin
    bin_size = args.bin_size
    df["bin_start"] = (df["begin"] // bin_size) * bin_size

    # Average mod_score per chromosome per bin
    avg = df.groupby(["chrom", "bin_start"])["mod_score"].mean().reset_index()
    avg.columns = ["chrom", "bin_start", "mean_mod_score"]
    avg["bin_end"] = avg["bin_start"] + bin_size

    # Clip bin_end to chromosome length
    chrom_len = fai.set_index("chrom")["length"].to_dict()
    avg["bin_end"] = avg.apply(
        lambda r: min(r["bin_end"], chrom_len.get(r["chrom"], r["bin_end"])), axis=1
    )

    # Reorder columns
    avg = avg[["chrom", "bin_start", "bin_end", "mean_mod_score"]]
    avg = avg.sort_values(["chrom", "bin_start"]).reset_index(drop=True)
    avg.to_csv(args.output_tsv, sep='\t', index=False)


if __name__ == "__main__":
    main()
