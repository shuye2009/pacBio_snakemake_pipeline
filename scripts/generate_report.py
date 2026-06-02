#!/usr/bin/env python
"""
Generate HTML summary report combining all visualizations.
"""

import argparse
import pandas as pd
import base64
from datetime import datetime
import os


def img_to_base64(img_path):
    """Encode image to base64 string."""
    with open(img_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def main():
    parser = argparse.ArgumentParser(description='Generate methylation analysis HTML report')
    parser.add_argument('--heatmap', required=True, help='Path to heatmap PNG')
    parser.add_argument('--dmr-heatmap', required=True, help='Path to significant DMR heatmap PNG')
    parser.add_argument('--volcano', required=True, help='Path to volcano plot PNG')
    parser.add_argument('--pca', required=True, help='Path to PCA plot PNG')
    parser.add_argument('--distribution', required=True, help='Path to distribution plot PNG')
    parser.add_argument('--signature-stats', required=True, help='Path to signature stats TSV')
    parser.add_argument('--dmr-stats', required=True, help='Path to dmr_stats_with_pvalues.tsv')
    parser.add_argument('--region-stats', required=True, help='Path to cohort_comparison_with_pvalues.tsv')
    parser.add_argument('--global-heatmap', required=False, help='Path to global methylation heatmap PNG')
    parser.add_argument('--global-pca', required=False, help='Path to global methylation PCA PNG')
    parser.add_argument('--global-density-composite', required=False, help='Path to composite density plot PNG')
    parser.add_argument('--global-density-per-sample', required=False, nargs='+', help='Paths to per-sample density PNGs')
    parser.add_argument('--global-density-per-sample-samples', required=False, nargs='+', help='Sample names for per-sample density plots')
    parser.add_argument('--global-density-per-chrom', required=False, nargs='+', help='Paths to per-chromosome density PNGs')
    parser.add_argument('--global-density-per-chrom-samples', required=False, nargs='+', help='Sample names for per-chromosome plots')
    parser.add_argument('--case-samples', required=True, help='Comma-separated list of case sample names')
    parser.add_argument('--control-samples', required=True, help='Comma-separated list of control sample names')
    parser.add_argument('--output', required=True, help='Output HTML file path')
    args = parser.parse_args()

    case_samples = args.case_samples.split(",")
    control_samples = args.control_samples.split(",")
    n_case = len(case_samples)
    n_control = len(control_samples)

    # Load DMR statistics and count significant DMRs
    try:
        dmr_stats = pd.read_csv(args.dmr_stats, sep='\t', comment='#')
    except Exception:
        dmr_stats = pd.read_csv(args.dmr_stats, sep='\t', on_bad_lines='skip')
    
    # Count significant DMRs (those shown in DMR heatmap)
    if 'significant' in dmr_stats.columns:
        n_sig_dmrs = dmr_stats['significant'].sum()
    else:
        n_sig_dmrs = len(dmr_stats)
    n_total_dmrs = len(dmr_stats)
    
    # Load predefined region statistics and count significant regions
    try:
        region_stats = pd.read_csv(args.region_stats, sep='\t', comment='#')
    except Exception:
        region_stats = pd.read_csv(args.region_stats, sep='\t', on_bad_lines='skip')
    
    # Count significant predefined regions (those shown in region heatmap)
    if 'significant' in region_stats.columns:
        n_sig_regions = region_stats['significant'].sum()
    else:
        n_sig_regions = len(region_stats)
    n_total_regions = len(region_stats)

    # Generate sample table rows
    case_rows = ''.join(f'<tr><td>{s}</td><td class="case">Case</td></tr>' for s in case_samples)
    control_rows = ''.join(f'<tr><td>{s}</td><td class="control">Control</td></tr>' for s in control_samples)

    # Generate HTML
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PacBio Methylation Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background-color: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .summary {{ background-color: #ecf0f1; padding: 20px; border-radius: 5px; margin: 20px 0; }}
        .summary-item {{ display: inline-block; margin-right: 40px; }}
        .summary-value {{ font-size: 24px; font-weight: bold; color: #2980b9; }}
        .summary-label {{ font-size: 14px; color: #7f8c8d; }}
        img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; margin: 10px 0; }}
        .figure {{ margin: 20px 0; }}
        .figure-caption {{ font-style: italic; color: #666; margin-top: 5px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        .case {{ color: #e74c3c; font-weight: bold; }}
        .control {{ color: #3498db; font-weight: bold; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>PacBio Methylation Analysis Report</h1>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <div class="summary">
            <div class="summary-item">
                <div class="summary-value">{n_case}</div>
                <div class="summary-label">Case Samples</div>
            </div>
            <div class="summary-item">
                <div class="summary-value">{n_control}</div>
                <div class="summary-label">Control Samples</div>
            </div>
            <div class="summary-item">
                <div class="summary-value">{n_sig_dmrs}</div>
                <div class="summary-label">Significant DMRs (of {n_total_dmrs})</div>
            </div>
            <div class="summary-item">
                <div class="summary-value">{n_sig_regions}</div>
                <div class="summary-label">Significant Regions (of {n_total_regions})</div>
            </div>
        </div>
        
        <h2>Sample Information</h2>
        <table>
            <tr><th>Sample</th><th>Group</th></tr>
            {case_rows}
            {control_rows}
        </table>
        
        <h2>Sample Clustering (PCA)</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.pca)}" alt="PCA Plot">
            <div class="figure-caption">Principal Component Analysis of methylation profiles showing sample clustering.</div>
        </div>
        
        <h2>Differentially Methylated Regions</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.volcano)}" alt="Volcano Plot">
            <div class="figure-caption">Volcano plot showing methylation differences vs. statistical significance.</div>
        </div>
        
        <h2>Methylation Heatmap (Predefined Regions)</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.heatmap)}" alt="Methylation Heatmap">
            <div class="figure-caption">Heatmap of methylation levels at significant predefined regions across all samples.</div>
        </div>
        
        <h2>Significant DMR Heatmap</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.dmr_heatmap)}" alt="Significant DMR Heatmap">
            <div class="figure-caption">Heatmap of methylation levels at significantly differentially methylated regions (from signature analysis).</div>
        </div>
        
        <h2>Methylation Distribution</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.distribution)}" alt="Methylation Distribution">
            <div class="figure-caption">Distribution of methylation levels for each sample.</div>
        </div>
        
        <h2>Global Methylation Analysis (50kb bins)</h2>
"""

    # Global heatmap
    if args.global_heatmap and os.path.exists(args.global_heatmap):
        html += f"""
        <h3>Global Methylation Heatmap</h3>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.global_heatmap)}" alt="Global Methylation Heatmap">
            <div class="figure-caption">Heatmap of mean modification scores across 50kb genomic bins for all samples.</div>
        </div>
"""

    # Global PCA
    if args.global_pca and os.path.exists(args.global_pca):
        html += f"""
        <h3>Global Methylation PCA</h3>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.global_pca)}" alt="Global Methylation PCA">
            <div class="figure-caption">PCA of samples based on genome-wide 50kb bin methylation profiles.</div>
        </div>
"""

    # Composite density
    if args.global_density_composite and os.path.exists(args.global_density_composite):
        html += f"""
        <h3>Composite Density Distribution</h3>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.global_density_composite)}" alt="Composite Density">
            <div class="figure-caption">Density distributions of mod_score, coverage, and estimated modified count for all samples overlaid.</div>
        </div>
"""

    # Per-sample density plots
    if args.global_density_per_sample and args.global_density_per_sample_samples:
        html += """
        <h3>Per-Sample Density Distributions</h3>
"""
        for png_path, sample in zip(args.global_density_per_sample, args.global_density_per_sample_samples):
            if os.path.exists(png_path):
                html += f"""
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(png_path)}" alt="Density - {sample}">
            <div class="figure-caption">Density distributions of mod_score, coverage, and estimated modified count for {sample}.</div>
        </div>
"""

    # Per-chromosome density plots
    if args.global_density_per_chrom and args.global_density_per_chrom_samples:
        html += """
        <h3>Per-Chromosome Density (mod_score)</h3>
"""
        for png_path, sample in zip(args.global_density_per_chrom, args.global_density_per_chrom_samples):
            if os.path.exists(png_path):
                html += f"""
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(png_path)}" alt="Per-Chromosome Density - {sample}">
            <div class="figure-caption">Per-chromosome mod_score density for {sample}.</div>
        </div>
"""

    html += """
    </div>
</body>
</html>
"""

    with open(args.output, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    main()
