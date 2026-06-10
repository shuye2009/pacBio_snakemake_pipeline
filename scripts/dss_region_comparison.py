#!/usr/bin/env python3
"""
DSS-based region methylation comparison using DSS results.

For each region BED file, computes:
  - Mean methylation difference (case - control) across CpGs in the region
  - Number of CpGs in the region
  - Number of DSS DMRs overlapping the region (total, hyper, hypo)
  - Number of DSS DMLs overlapping the region (total, hyper, hypo)

Supports multiple case and control samples (CpGs pooled per group).
Outputs a single TSV with one row per region.
"""

import argparse
import pandas as pd
import numpy as np
import pysam
import os
import sys


def load_cpg_methylation_multi(bed_gz_files):
    """Load and pool CpG methylation from multiple combined.bed.gz files.

    Returns DataFrame with columns: chrom, pos, meth (0-1 scale)
    """
    frames = []
    for bed_gz_file in bed_gz_files:
        records = []
        try:
            tbx = pysam.TabixFile(bed_gz_file)
        except Exception as e:
            print(f"Error opening {bed_gz_file}: {e}", file=sys.stderr)
            continue
        try:
            for contig in tbx.contigs:
                try:
                    for row in tbx.fetch(contig):
                        if row.startswith('#'):
                            continue
                        fields = row.split('\t')
                        if len(fields) >= 4:
                            try:
                                meth = float(fields[3]) / 100.0
                                records.append({
                                    'chrom': fields[0],
                                    'pos': int(fields[1]) + 1,
                                    'meth': meth
                                })
                            except (ValueError, IndexError):
                                continue
                except ValueError:
                    continue
        finally:
            tbx.close()
        if records:
            frames.append(pd.DataFrame(records))
    if not frames:
        return pd.DataFrame(columns=['chrom', 'pos', 'meth'])
    return pd.concat(frames, ignore_index=True)


def load_bed_regions(bed_file):
    """Load regions from a BED file.

    Returns DataFrame with columns: chrom, start, end, name
    """
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
    """Load DSS DMRs with delta values.

    Returns DataFrame with: chrom, start, end, delta
    """
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
    """Load DSS DMLs with delta values.

    Returns DataFrame with: chrom, pos, delta
    """
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


def count_overlaps(regions_df, features_df, feature_type='dmr'):
    """Count overlapping features per region.

    features_df must have chrom, start, end, delta columns.
    Returns DataFrame with overlap counts.
    """
    results = []
    for _, region in regions_df.iterrows():
        r_chrom = region['chrom']
        r_start = region['start']
        r_end = region['end']

        overlaps = features_df[
            (features_df['chrom'] == r_chrom) &
            (features_df['start'] < r_end) &
            (features_df['end'] > r_start)
        ]

        total = len(overlaps)
        hyper = int((overlaps['delta'] > 0).sum()) if 'delta' in overlaps.columns else 0
        hypo = int((overlaps['delta'] < 0).sum()) if 'delta' in overlaps.columns else 0

        results.append({
            f'n_{feature_type}s': total,
            f'n_hyper_{feature_type}s': hyper,
            f'n_hypo_{feature_type}s': hypo,
        })

    return pd.DataFrame(results)


def compute_region_methylation(regions_df, case_cpg, control_cpg):
    """Compute mean methylation difference per region.

    Returns DataFrame with: n_cpgs, case_mean_meth, control_mean_meth, meth_diff
    """
    results = []
    for _, region in regions_df.iterrows():
        r_chrom = region['chrom']
        r_start = region['start']
        r_end = region['end']

        case_in = case_cpg[
            (case_cpg['chrom'] == r_chrom) &
            (case_cpg['pos'] > r_start) &
            (case_cpg['pos'] <= r_end)
        ]
        ctrl_in = control_cpg[
            (control_cpg['chrom'] == r_chrom) &
            (control_cpg['pos'] > r_start) &
            (control_cpg['pos'] <= r_end)
        ]

        n_cpgs = max(len(case_in), len(ctrl_in))
        case_mean = case_in['meth'].mean() if not case_in.empty else np.nan
        ctrl_mean = ctrl_in['meth'].mean() if not ctrl_in.empty else np.nan
        meth_diff = case_mean - ctrl_mean if not (np.isnan(case_mean) or np.isnan(ctrl_mean)) else np.nan

        results.append({
            'n_cpgs': n_cpgs,
            'case_mean_meth': case_mean,
            'control_mean_meth': ctrl_mean,
            'meth_diff': meth_diff,
        })

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="DSS-based region methylation comparison (1 case vs 1 control)"
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
    args = parser.parse_args()

    if len(args.region_bed) != len(args.region_name):
        print("Error: number of --region-bed and --region-name must match", file=sys.stderr)
        raise SystemExit(1)

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    case_files = args.case_beds.split(',')
    control_files = args.control_beds.split(',')

    print(f"Loading CpG methylation data ({len(case_files)} case, {len(control_files)} control samples)...")
    case_cpg = load_cpg_methylation_multi(case_files)
    control_cpg = load_cpg_methylation_multi(control_files)
    print(f"  Case: {len(case_cpg)} CpGs, Control: {len(control_cpg)} CpGs")

    print("Loading DSS DMRs and DMLs...")
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

        meth_stats = compute_region_methylation(regions, case_cpg, control_cpg)

        dmr_counts = count_overlaps(regions, dmrs, 'dmr') if not dmrs.empty else pd.DataFrame(
            {'n_dmrs': 0, 'n_hyper_dmrs': 0, 'n_hypo_dmrs': 0}, index=regions.index)

        dml_counts = count_overlaps(regions, dmls, 'dml') if not dmls.empty else pd.DataFrame(
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


if __name__ == '__main__':
    main()
