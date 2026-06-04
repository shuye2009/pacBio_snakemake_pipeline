#!/usr/bin/env python3
"""
Compute annotation overlap counts per genomic bin.

Tiles the genome into fixed-size bins and counts, for each bin, the
number of annotation intervals (genes, repeats, centromeres, enhancers,
CGIs, significant DMRs) that intersect the bin range.
"""

import argparse
import pandas as pd
import numpy as np
import os


def load_bed(path):
    """Load a BED file, returning DataFrame with chrom, start, end."""
    df = pd.read_csv(path, sep='\t', header=None, comment='#',
                     usecols=[0, 1, 2], names=["chrom", "start", "end"])
    return df


def load_gtf_genes(path):
    """Load gene regions from a GTF file."""
    df = pd.read_csv(path, sep='\t', header=None, comment='#',
                     usecols=[0, 2, 3, 4],
                     names=["chrom", "feature", "start", "end"])
    df = df[df["feature"] == "gene"].copy()
    df["start"] = df["start"] - 1  # GTF is 1-based, convert to 0-based
    return df[["chrom", "start", "end"]]


def count_overlaps(bin_start, bin_end, starts, ends):
    """Count how many intervals intersect [bin_start, bin_end).

    Uses vectorized range intersection: interval overlaps if
    interval.start < bin_end AND interval.end > bin_start.
    """
    return np.sum((starts < bin_end) & (ends > bin_start))


def main():
    parser = argparse.ArgumentParser(
        description="Compute annotation overlap counts per genomic bin"
    )
    parser.add_argument("--fai", required=True, help="Genome FASTA index (.fai)")
    parser.add_argument("--bin-size", type=int, default=500000, help="Bin size in bp")
    parser.add_argument("--genes", required=False, help="GTF file for gene regions")
    parser.add_argument("--repeats", required=False, help="BED file for repeats")
    parser.add_argument("--centromeres", required=False, help="BED file for centromeres")
    parser.add_argument("--enhancers", required=False, help="BED file for enhancers")
    parser.add_argument("--cgi", required=False, help="BED file for CpG islands")
    parser.add_argument("--dmrs", required=False, help="BED file for significant DMRs")
    parser.add_argument("--dml", required=False, help="BED file for DSS DML results")
    parser.add_argument("--dmr", required=False, help="BED file for DSS DMR results")
    parser.add_argument("--output-tsv", required=True, help="Output TSV path")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    # Load chromosome sizes
    fai = pd.read_csv(args.fai, sep='\t', header=None,
                      names=["chrom", "length", "offset", "linebases", "linewidth"])

    # Load annotations and pre-group by chromosome into numpy arrays
    annotations = {}
    if args.genes:
        df = load_gtf_genes(args.genes)
        annotations["gene_count"] = {
            c: (g["start"].values, g["end"].values)
            for c, g in df.groupby("chrom")
        }
    if args.repeats:
        df = load_bed(args.repeats)
        annotations["repeat_count"] = {
            c: (g["start"].values, g["end"].values)
            for c, g in df.groupby("chrom")
        }
    if args.centromeres:
        df = load_bed(args.centromeres)
        annotations["centromere_count"] = {
            c: (g["start"].values, g["end"].values)
            for c, g in df.groupby("chrom")
        }
    if args.enhancers:
        df = load_bed(args.enhancers)
        annotations["enhancer_count"] = {
            c: (g["start"].values, g["end"].values)
            for c, g in df.groupby("chrom")
        }
    if args.cgi:
        df = load_bed(args.cgi)
        annotations["cgi_count"] = {
            c: (g["start"].values, g["end"].values)
            for c, g in df.groupby("chrom")
        }
    if args.dmrs:
        df = load_bed(args.dmrs)
        annotations["dmr_count"] = {
            c: (g["start"].values, g["end"].values)
            for c, g in df.groupby("chrom")
        }
    if args.dml:
        df = load_bed(args.dml)
        annotations["dml_count"] = {
            c: (g["start"].values, g["end"].values)
            for c, g in df.groupby("chrom")
        }
    if args.dmr:
        df = load_bed(args.dmr)
        annotations["dss_dmr_count"] = {
            c: (g["start"].values, g["end"].values)
            for c, g in df.groupby("chrom")
        }

    if not annotations:
        print("Warning: No annotation tracks provided, output will be empty.")
        bins = []
        for _, row in fai.iterrows():
            chrom = row["chrom"]
            length = row["length"]
            for start in range(0, length, args.bin_size):
                bins.append({
                    "chrom": chrom,
                    "bin_start": start,
                    "bin_end": min(start + args.bin_size, length)
                })
        bins_df = pd.DataFrame(bins)
        bins_df.to_csv(args.output_tsv, sep='\t', index=False)
        return

    # Build bins with vectorized overlap counts
    rows = []
    for _, chrom_row in fai.iterrows():
        chrom = chrom_row["chrom"]
        length = chrom_row["length"]
        for start in range(0, length, args.bin_size):
            end = min(start + args.bin_size, length)
            row = {"chrom": chrom, "bin_start": start, "bin_end": end}
            for col_name, chrom_dict in annotations.items():
                arrs = chrom_dict.get(chrom)
                if arrs is None:
                    row[col_name] = 0
                else:
                    row[col_name] = count_overlaps(start, end, arrs[0], arrs[1])
            rows.append(row)

    result = pd.DataFrame(rows)
    cols = ["chrom", "bin_start", "bin_end"] + list(annotations.keys())
    result = result[cols]
    result.to_csv(args.output_tsv, sep='\t', index=False)


if __name__ == "__main__":
    main()
