#!/usr/bin/env python3
"""
Generate IGV track configuration JSON for igv-reports.

Creates a JSON file with track configurations for BAM files,
enabling haplotype grouping and methylation coloring.
"""

import argparse
import json
import os


def main():
    parser = argparse.ArgumentParser(description='Generate IGV track config JSON')
    parser.add_argument('--bams', nargs='+', required=True, help='BAM file paths')
    parser.add_argument('--output', required=True, help='Output JSON file path')
    parser.add_argument('--phasing-enabled', action='store_true', help='Enable haplotype grouping by HP tag')
    parser.add_argument('--gtf', default=None, help='GTF annotation file path (bgzipped, with .tbi index)')
    parser.add_argument('--gtf-index', default=None, help='GTF tabix index path')
    args = parser.parse_args()

    tracks = []
    for bam in args.bams:
        sample_name = os.path.basename(bam).replace('.haplotagged.bam', '').replace('.aligned.bam', '')
        track = {
            "name": sample_name,
            "url": bam,
            "indexURL": bam + ".bai",
            "format": "bam",
            "type": "alignment",
            "colorBy": "basemod",
            "basemodColorScheme": "5mC",
            "showSoftClips": False,
            "viewAsPairs": False,
            "height": 300,
            "displayMode": "EXPANDED",
            "showMismatches": True,
            "showInsertions": True
        }
        if args.phasing_enabled:
            track["groupBy"] = "tag:HP"
        tracks.append(track)

    if args.gtf:
        gtf_track = {
            "name": "Genes",
            "url": args.gtf,
            "indexURL": args.gtf_index or (args.gtf + ".tbi"),
            "format": "gtf",
            "type": "annotation",
            "height": 50,
            "displayMode": "EXPANDED"
        }
        tracks.append(gtf_track)

    with open(args.output, 'w') as f:
        json.dump(tracks, f, indent=2)
    
    print(f"Generated track config with {len(tracks)} BAM tracks: {args.output}")


if __name__ == "__main__":
    main()
