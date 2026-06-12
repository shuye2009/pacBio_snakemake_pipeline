#!/usr/bin/env python3
"""
DSS-based region methylation comparison using bedtools for speed.

For each region BED file, computes:
  - Mean methylation difference (case - control) across CpGs in the region
  - Number of CpGs in the region
  - Number of DSS DMRs overlapping the region (total, hyper, hypo)
  - Number of DSS DMLs overlapping the region (total, hyper, hypo)

Uses bedtools map for per-sample per-region methylation aggregation
and bedtools intersect for overlap counting.
"""

import argparse
import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed


def bedtools_map(bed_a, bed_b, columns=5, operations="mean,count"):
    """Run bedtools map and return parsed DataFrame.

    Returns DataFrame with columns: chrom, start, end, name, [op1, op2, ...]
    """
    cmd = [
        "bedtools", "map",
        "-a", bed_a,
        "-b", bed_b,
        "-c", str(columns),
        "-o", operations,
        "-null", "NA",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"bedtools map failed: {result.stderr}", file=sys.stderr)
        return None

    lines = [l for l in result.stdout.strip().split('\n') if l.strip()]
    if not lines:
        return None

    op_names = operations.split(',')
    col_names = ['chrom', 'start', 'end', 'name'] + op_names
    rows = []
    for line in lines:
        fields = line.split('\t')
        rows.append(fields[:len(col_names)])

    df = pd.DataFrame(rows, columns=col_names)
    for op in op_names:
        df[op] = pd.to_numeric(df[op], errors='coerce')
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    return df


def bedtools_intersect_count(bed_a, bed_b):
    """Count overlaps of bed_b features in each region of bed_a.

    Returns list of counts (same order as bed_a).
    """
    cmd = [
        "bedtools", "intersect",
        "-a", bed_a,
        "-b", bed_b,
        "-c",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"bedtools intersect failed: {result.stderr}", file=sys.stderr)
        return []

    counts = []
    for line in result.stdout.strip().split('\n'):
        if line.strip():
            fields = line.split('\t')
            counts.append(int(fields[-1]) if fields[-1].isdigit() else 0)
    return counts


def compute_region_methylation_bedtools(regions_df, case_files, control_files, threads=1):
    """Compute per-region methylation stats using bedtools map.

    Runs bedtools map per sample in parallel, then aggregates across case/control groups.
    Returns DataFrame with: n_cpgs, case_mean_meth, control_mean_meth, meth_diff
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        regions_df[['chrom', 'start', 'end', 'name']].to_csv(
            f, sep='\t', index=False, header=False)
        regions_tmp = f.name

    try:
        all_files = list(case_files) + list(control_files)
        n_case = len(case_files)

        # Run bedtools map in parallel for all samples
        results_by_file = {}
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(bedtools_map, regions_tmp, f): f for f in all_files}
            for future in as_completed(futures):
                bed_file = futures[future]
                try:
                    df = future.result()
                    if df is not None:
                        results_by_file[bed_file] = df
                except Exception as e:
                    print(f"bedtools map failed for {bed_file}: {e}", file=sys.stderr)

        # Collect per-sample means
        case_means = []
        control_means = []
        all_cpg_counts = []

        for bed_file in case_files:
            if bed_file in results_by_file:
                df = results_by_file[bed_file]
                case_means.append(df.set_index('name')['mean'])
                all_cpg_counts.append(df.set_index('name')['count'])

        for bed_file in control_files:
            if bed_file in results_by_file:
                df = results_by_file[bed_file]
                control_means.append(df.set_index('name')['mean'])
                all_cpg_counts.append(df.set_index('name')['count'])

        # Aggregate
        case_mean = pd.concat(case_means, axis=1).mean(axis=1) if case_means else pd.Series(dtype=float)
        control_mean = pd.concat(control_means, axis=1).mean(axis=1) if control_means else pd.Series(dtype=float)
        max_cpgs = pd.concat(all_cpg_counts, axis=1).max(axis=1) if all_cpg_counts else pd.Series(dtype=int)

        # Align to regions
        result = pd.DataFrame(index=regions_df['name'])
        result['n_cpgs'] = max_cpgs.reindex(result.index).fillna(0).astype(int)
        result['case_mean_meth'] = case_mean.reindex(result.index)
        result['control_mean_meth'] = control_mean.reindex(result.index)
        result['meth_diff'] = result['case_mean_meth'] - result['control_mean_meth']

        return result
    finally:
        os.unlink(regions_tmp)


def count_overlaps_bedtools(regions_df, features_df, feature_type='dmr'):
    """Count overlapping features per region using bedtools intersect.

    Returns DataFrame with overlap counts.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f_reg:
        regions_df[['chrom', 'start', 'end', 'name']].to_csv(
            f_reg, sep='\t', index=False, header=False)
        regions_tmp = f_reg.name

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f_feat:
        features_df[['chrom', 'start', 'end']].to_csv(
            f_feat, sep='\t', index=False, header=False)
        features_tmp = f_feat.name

    try:
        total_counts = bedtools_intersect_count(regions_tmp, features_tmp)

        # Hyper/hypo counts
        hyper = features_df[features_df['delta'] > 0] if 'delta' in features_df.columns else features_df.iloc[:0]
        hypo = features_df[features_df['delta'] < 0] if 'delta' in features_df.columns else features_df.iloc[:0]

        if not hyper.empty:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f_h:
                hyper[['chrom', 'start', 'end']].to_csv(
                    f_h, sep='\t', index=False, header=False)
                hyper_tmp = f_h.name
            hyper_counts = bedtools_intersect_count(regions_tmp, hyper_tmp)
            os.unlink(hyper_tmp)
        else:
            hyper_counts = [0] * len(regions_df)

        if not hypo.empty:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f_h:
                hypo[['chrom', 'start', 'end']].to_csv(
                    f_h, sep='\t', index=False, header=False)
                hypo_tmp = f_h.name
            hypo_counts = bedtools_intersect_count(regions_tmp, hypo_tmp)
            os.unlink(hypo_tmp)
        else:
            hypo_counts = [0] * len(regions_df)

        result = pd.DataFrame({
            f'n_{feature_type}s': total_counts,
            f'n_hyper_{feature_type}s': hyper_counts,
            f'n_hypo_{feature_type}s': hypo_counts,
        })
        return result
    finally:
        os.unlink(regions_tmp)
        os.unlink(features_tmp)


def load_bed_regions(bed_file):
    """Load regions from a BED file."""
    if not os.path.exists(bed_file) or os.path.getsize(bed_file) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'name'])

    df = pd.read_csv(bed_file, sep='\t', header=None, comment='#')
    cols = min(4, df.shape[1])
    names = ['chrom', 'start', 'end', 'name'][:cols]
    df = df.iloc[:, :cols]
    df.columns = names
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    if 'name' not in df.columns:
        df['name'] = df['chrom'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    return df


def load_dmr_bed(bed_file, tsv_file):
    """Load DSS DMRs with delta values."""
    bed = load_bed_regions(bed_file)
    if bed.empty:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'delta'])

    tsv = pd.read_csv(tsv_file, sep='\t')
    diff_col = None
    for col in ['diff.Methy', 'diff_Methy', 'diff', 'DELTA']:
        if col in tsv.columns:
            diff_col = col
            break

    if diff_col:
        bed['delta'] = 0.0
        for i, row in bed.iterrows():
            chrom, start, end = row['chrom'], row['start'], row['end']
            match = tsv[
                (tsv.get('chr', tsv.get('CHROM', '')) == chrom) &
                (tsv.get('start', tsv.get('START', 0)) == start + 1) &
                (tsv.get('end', tsv.get('END', 0)) == end)
            ]
            if not match.empty:
                bed.at[i, 'delta'] = float(match[diff_col].values[0])
    else:
        bed['delta'] = 0.0

    return bed[['chrom', 'start', 'end', 'delta']]


def load_dml_bed(bed_file, tsv_file):
    """Load DSS DMLs with delta values."""
    bed = load_bed_regions(bed_file)
    if bed.empty:
        return pd.DataFrame(columns=['chrom', 'pos', 'delta'])

    tsv = pd.read_csv(tsv_file, sep='\t')
    diff_col = None
    for col in ['diff', 'diff.Methy', 'diff_Methy', 'DELTA']:
        if col in tsv.columns:
            diff_col = col
            break

    if diff_col:
        bed['delta'] = 0.0
        for i, row in bed.iterrows():
            chrom, pos = row['chrom'], row['start'] + 1
            match = tsv[
                (tsv.get('chr', tsv.get('CHROM', '')) == chrom) &
                (tsv.get('pos', tsv.get('POS', 0)) == pos)
            ]
            if not match.empty:
                bed.at[i, 'delta'] = float(match[diff_col].values[0])
    else:
        bed['delta'] = 0.0

    return bed[['chrom', 'start', 'delta']].rename(columns={'start': 'pos'})


def main():
    parser = argparse.ArgumentParser(
        description="DSS-based region methylation comparison using bedtools"
    )
    parser.add_argument('--case-beds', required=True,
                        help='Comma-separated case combined.bed.gz files')
    parser.add_argument('--control-beds', required=True,
                        help='Comma-separated control combined.bed.gz files')
    parser.add_argument('--dmr-bed', required=True, help='DSS DMR BED file')
    parser.add_argument('--dmr-tsv', required=True, help='DSS DMR TSV file')
    parser.add_argument('--dml-bed', required=True, help='DSS DML BED file')
    parser.add_argument('--dml-tsv', required=True, help='DSS DML TSV file')
    parser.add_argument('--region-bed', required=True, action='append', default=[],
                        help='Region BED file (can specify multiple)')
    parser.add_argument('--region-name', required=True, action='append', default=[],
                        help='Region name for each --region-bed')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel bedtools calls')
    args = parser.parse_args()

    if len(args.region_bed) != len(args.region_name):
        print("Error: number of --region-bed and --region-name must match", file=sys.stderr)
        raise SystemExit(1)

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    case_files = args.case_beds.split(',')
    control_files = args.control_beds.split(',')

    print(f"Loading DSS DMRs and DMLs...")
    dmrs = load_dmr_bed(args.dmr_bed, args.dmr_tsv)
    dmls = load_dml_bed(args.dml_bed, args.dml_tsv)
    print(f"  DMRs: {len(dmrs)}, DMLs: {len(dmls)}")

    all_results = []

    for region_bed, region_name in zip(args.region_bed, args.region_name):
        print(f"Processing {region_name}: {region_bed}")
        regions = load_bed_regions(region_bed)
        if regions.empty:
            print(f"  Skipping {region_name}: no regions found")
            continue

        print(f"  Loaded {len(regions)} regions")

        print(f"  Computing methylation stats via bedtools map ({len(case_files)} case, {len(control_files)} control samples, {args.threads} threads)...")
        meth_stats = compute_region_methylation_bedtools(regions, case_files, control_files, args.threads)

        print(f"  Counting DMR overlaps via bedtools intersect...")
        dmr_counts = count_overlaps_bedtools(regions, dmrs, 'dmr') if not dmrs.empty else pd.DataFrame(
            {'n_dmrs': 0, 'n_hyper_dmrs': 0, 'n_hypo_dmrs': 0}, index=regions.index)

        print(f"  Counting DML overlaps via bedtools intersect...")
        dml_counts = count_overlaps_bedtools(regions, dmls, 'dml') if not dmls.empty else pd.DataFrame(
            {'n_dmls': 0, 'n_hyper_dmls': 0, 'n_hypo_dmls': 0}, index=regions.index)

        combined = pd.concat([
            regions[['chrom', 'start', 'end', 'name']].reset_index(drop=True),
            meth_stats.reset_index(drop=True),
            dmr_counts.reset_index(drop=True),
            dml_counts.reset_index(drop=True),
        ], axis=1)

        combined['region_type'] = region_name
        all_results.append(combined)

    if not all_results:
        print("No regions processed. Creating empty output.")
        pd.DataFrame(columns=[
            'chrom', 'start', 'end', 'name', 'region_type',
            'n_cpgs', 'case_mean_meth', 'control_mean_meth', 'meth_diff',
            'n_dmrs', 'n_hyper_dmrs', 'n_hypo_dmrs',
            'n_dmls', 'n_hyper_dmls', 'n_hypo_dmls',
        ]).to_csv(args.output_tsv, sep='\t', index=False)
        return

    final = pd.concat(all_results, ignore_index=True)
    final.to_csv(args.output_tsv, sep='\t', index=False, float_format='%.6g')
    print(f"Wrote {len(final)} region comparisons to {args.output_tsv}")


if __name__ == "__main__":
    main()
