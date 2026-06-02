#!/usr/bin/env python3
"""
Compare ASM regions between case and control groups.

Computes overlaps between joint-segmented regions from case and control groups
and creates a confusion matrix showing overlaps by summary_label.
Uses bedtools intersect for fast overlap computation with Jaccard filtering.
"""

import argparse
import pandas as pd
import subprocess
import os
from collections import defaultdict


def run_bedtools_intersect(case_bed, control_bed):
    """
    Run bedtools intersect to find overlapping regions.
    Returns DataFrame with case regions, matched control regions, and overlap info.
    """
    # bedtools intersect -a case -b control -wao
    # -wao: Write all original A entries plus overlap info, including unmatched
    result = subprocess.run(
        ['bedtools', 'intersect', '-a', case_bed, '-b', control_bed, '-wao'],
        capture_output=True, text=True, check=True
    )
    
    # Parse output: case_chrom, case_start, case_end, case_label, ctrl_chrom, ctrl_start, ctrl_end, ctrl_label, overlap
    lines = result.stdout.strip().split('\n')
    rows = []
    for line in lines:
        if not line:
            continue
        fields = line.split('\t')
        rows.append({
            'case_chrom': fields[0],
            'case_start': int(fields[1]),
            'case_end': int(fields[2]),
            'case_label': fields[3],
            'control_chrom': fields[4] if fields[4] != '.' else None,
            'control_start': int(fields[5]) if fields[5] != '-1' else None,
            'control_end': int(fields[6]) if fields[6] != '-1' else None,
            'control_label': fields[7] if fields[7] != '.' else None,
            'overlap_bp': int(fields[8])
        })
    
    return pd.DataFrame(rows)


def compute_jaccard(case_start, case_end, ctrl_start, ctrl_end, overlap):
    """Compute Jaccard index from overlap info."""
    if overlap <= 0 or ctrl_start is None:
        return 0.0
    
    union = max(case_end, ctrl_end) - min(case_start, ctrl_start)
    return overlap / union if union > 0 else 0.0


def main():
    parser = argparse.ArgumentParser(description='Compare ASM regions between case and control groups')
    parser.add_argument('--case-bed', required=True, help='BED file with case ASM regions')
    parser.add_argument('--control-bed', required=True, help='BED file with control ASM regions')
    parser.add_argument('--output-comparison', required=True, help='Output TSV with detailed comparison')
    parser.add_argument('--output-confusion-matrix', required=True, help='Output TSV with confusion matrix')
    parser.add_argument('--output-case-specific', required=True, help='Output BED with case-specific regions')
    parser.add_argument('--output-control-specific', required=True, help='Output BED with control-specific regions')
    parser.add_argument('--output-shared', required=True, help='Output BED with shared regions')
    parser.add_argument('--min-jaccard', type=float, default=0.5, help='Minimum Jaccard index to count as shared (0-1)')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_comparison), exist_ok=True)

    print(f"Running bedtools intersect...")
    
    # Run bedtools intersect
    intersect_df = run_bedtools_intersect(args.case_bed, args.control_bed)
    
    print(f"Processing {len(intersect_df)} intersection results...")

    # Compute Jaccard index for each overlap
    intersect_df['jaccard_index'] = intersect_df.apply(
        lambda row: compute_jaccard(
            row['case_start'], row['case_end'],
            row['control_start'], row['control_end'],
            row['overlap_bp']
        ), axis=1
    )

    # For each case region, keep only the best overlap (highest Jaccard)
    # Group by case region and get the row with max Jaccard
    case_key_cols = ['case_chrom', 'case_start', 'case_end', 'case_label']
    
    # Sort by Jaccard descending and drop duplicates keeping first (best)
    intersect_df_sorted = intersect_df.sort_values('jaccard_index', ascending=False)
    best_matches = intersect_df_sorted.drop_duplicates(subset=case_key_cols, keep='first').copy()

    # Classify regions
    best_matches['status'] = best_matches.apply(
        lambda row: 'shared' if row['jaccard_index'] >= args.min_jaccard else 'case_specific',
        axis=1
    )
    best_matches['label_changed'] = best_matches.apply(
        lambda row: row['case_label'] != row['control_label'] if row['status'] == 'shared' else True,
        axis=1
    )

    # Round Jaccard for output
    best_matches['jaccard_index'] = best_matches['jaccard_index'].round(4)

    # Find control-specific regions (not matched by any case region meeting threshold)
    matched_control_keys = set()
    shared_df = best_matches[best_matches['status'] == 'shared']
    for _, row in shared_df.iterrows():
        if row['control_chrom'] is not None:
            matched_control_keys.add((row['control_chrom'], row['control_start'], row['control_end']))

    # Load control BED to find unmatched
    control_df = pd.read_csv(args.control_bed, sep='\t', header=None, 
                             names=['chrom', 'start', 'end', 'label'], comment='#')
    
    control_specific_rows = []
    for _, ctrl_row in control_df.iterrows():
        key = (ctrl_row['chrom'], ctrl_row['start'], ctrl_row['end'])
        if key not in matched_control_keys:
            control_specific_rows.append({
                'case_chrom': 'N/A',
                'case_start': 'N/A',
                'case_end': 'N/A',
                'case_label': 'N/A',
                'control_chrom': ctrl_row['chrom'],
                'control_start': ctrl_row['start'],
                'control_end': ctrl_row['end'],
                'control_label': ctrl_row['label'],
                'overlap_bp': 0,
                'jaccard_index': 0.0,
                'status': 'control_specific',
                'label_changed': True
            })

    # Combine results
    comparison_df = pd.concat([best_matches, pd.DataFrame(control_specific_rows)], ignore_index=True)

    # Count statistics
    n_shared = len(comparison_df[comparison_df['status'] == 'shared'])
    n_case_specific = len(comparison_df[comparison_df['status'] == 'case_specific'])
    n_control_specific = len(comparison_df[comparison_df['status'] == 'control_specific'])
    
    print(f"Shared: {n_shared}, Case-specific: {n_case_specific}, Control-specific: {n_control_specific}")

    # Build confusion matrix for shared regions
    confusion = defaultdict(lambda: defaultdict(int))
    label_types = ['Methylated', 'Unmethylated', 'AlleleSpecificMethylation', 'Uncategorized', 'NoData']

    for _, row in shared_df.iterrows():
        case_label = row['case_label']
        control_label = row['control_label']
        if control_label is not None:
            confusion[control_label][case_label] += 1

    # Write confusion matrix
    with open(args.output_confusion_matrix, 'w') as f:
        f.write("# Confusion matrix: rows=control labels, columns=case labels\n")
        f.write(f"# Based on Jaccard index >= {args.min_jaccard}\n")
        f.write("control_label\t" + "\t".join(label_types) + "\ttotal\n")
        for ctrl_label in label_types:
            row_counts = [str(confusion[ctrl_label][case_label]) for case_label in label_types]
            row_total = sum(confusion[ctrl_label].values())
            f.write(f"{ctrl_label}\t" + "\t".join(row_counts) + f"\t{row_total}\n")
        # Column totals
        col_totals = [str(sum(confusion[ctrl][case_label] for ctrl in label_types)) for case_label in label_types]
        f.write("total\t" + "\t".join(col_totals) + f"\t{n_shared}\n")

    # Write detailed comparison TSV
    comparison_df.to_csv(args.output_comparison, sep='\t', index=False)

    # Write BED files for each category
    # Case-specific
    case_specific_df = comparison_df[comparison_df['status'] == 'case_specific']
    with open(args.output_case_specific, 'w') as f:
        for _, row in case_specific_df.iterrows():
            f.write(f"{row['case_chrom']}\t{row['case_start']}\t{row['case_end']}\t{row['case_label']}\n")
    
    # Control-specific
    control_specific_df = comparison_df[comparison_df['status'] == 'control_specific']
    with open(args.output_control_specific, 'w') as f:
        for _, row in control_specific_df.iterrows():
            f.write(f"{row['control_chrom']}\t{row['control_start']}\t{row['control_end']}\t{row['control_label']}\n")
    
    # Shared (use case coordinates)
    with open(args.output_shared, 'w') as f:
        for _, row in shared_df.iterrows():
            f.write(f"{row['case_chrom']}\t{row['case_start']}\t{row['case_end']}\t{row['case_label']}\n")

    print(f"Wrote comparison to {args.output_comparison}")
    print(f"Wrote confusion matrix to {args.output_confusion_matrix}")


if __name__ == "__main__":
    main()
