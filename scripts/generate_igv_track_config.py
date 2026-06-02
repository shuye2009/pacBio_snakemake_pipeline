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
    args = parser.parse_args()

    tracks = []
    for bam in args.bams:
        sample_name = os.path.basename(bam).replace('.haplotagged.bam', '')
        track = {
            "name": sample_name,
            "url": bam,
            "indexURL": bam + ".bai",
            "format": "bam",
            "type": "alignment",
            "colorBy": "basemod",
            "basemodColorScheme": "5mC",
            "groupBy": "tag:HP",
            "showSoftClips": False,
            "viewAsPairs": False,
            "height": 300,
            "displayMode": "EXPANDED",
            "showMismatches": True,
            "showInsertions": True
        }
        tracks.append(track)

    with open(args.output, 'w') as f:
        json.dump(tracks, f, indent=2)
    
    print(f"Generated track config with {len(tracks)} BAM tracks: {args.output}")


if __name__ == "__main__":
    main()
