#!/usr/bin/env python
"""
Generate HTML summary report combining all visualizations.
"""

import argparse
import pandas as pd
import base64
import html
from datetime import datetime
import os


def img_to_base64(img_path):
    """Encode image to base64 string."""
    with open(img_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def main():
    parser = argparse.ArgumentParser(description='Generate methylation analysis HTML report')
    parser.add_argument('--heatmap', required=False, help='Path to heatmap PNG')
    parser.add_argument('--dmr-heatmap', required=False, help='Path to significant DMR heatmap PNG')
    parser.add_argument('--volcano', required=False, help='Path to volcano plot PNG')
    parser.add_argument('--pca', required=True, help='Path to PCA plot PNG')
    parser.add_argument('--distribution', required=True, help='Path to distribution plot PNG')
    parser.add_argument('--signature-stats', required=False, help='Path to signature stats TSV')
    parser.add_argument('--dmr-stats', required=False, help='Path to dmr_stats_with_pvalues.tsv')
    parser.add_argument('--region-stats', required=False, help='Path to cohort_comparison_with_pvalues.tsv')
    parser.add_argument('--dss-dmr', required=False, help='Path to DSS DMR results TSV')
    parser.add_argument('--dss-dml', required=False, help='Path to DSS DML results TSV')
    parser.add_argument('--dss-dmr-heatmap', required=False, help='Path to DSS DMR heatmap PNG')
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

    # Load DMR statistics and count significant DMRs (only if available)
    n_sig_dmrs = 0
    n_total_dmrs = 0
    n_sig_regions = 0
    n_total_regions = 0
    has_dmr_stats = args.dmr_stats and os.path.exists(args.dmr_stats)
    has_region_stats = args.region_stats and os.path.exists(args.region_stats)

    if has_dmr_stats:
        try:
            dmr_stats = pd.read_csv(args.dmr_stats, sep='\t', comment='#')
        except pd.errors.ParserError:
            dmr_stats = pd.read_csv(args.dmr_stats, sep='\t', on_bad_lines='skip')
        if 'significant' in dmr_stats.columns:
            n_sig_dmrs = dmr_stats['significant'].sum()
        else:
            n_sig_dmrs = len(dmr_stats)
        n_total_dmrs = len(dmr_stats)

    if has_region_stats:
        try:
            region_stats = pd.read_csv(args.region_stats, sep='\t', comment='#')
        except pd.errors.ParserError:
            region_stats = pd.read_csv(args.region_stats, sep='\t', on_bad_lines='skip')
        if 'significant' in region_stats.columns:
            n_sig_regions = region_stats['significant'].sum()
        else:
            n_sig_regions = len(region_stats)
        n_total_regions = len(region_stats)

    # Load DSS results
    n_dss_dmrs = 0
    n_dss_dmls = 0
    has_dss_dmr = args.dss_dmr and os.path.exists(args.dss_dmr)
    has_dss_dml = args.dss_dml and os.path.exists(args.dss_dml)

    if has_dss_dmr:
        try:
            dss_dmr = pd.read_csv(args.dss_dmr, sep='\t', comment='#')
            n_dss_dmrs = len(dss_dmr)
        except Exception:
            n_dss_dmrs = 0

    if has_dss_dml:
        try:
            dss_dml = pd.read_csv(args.dss_dml, sep='\t', comment='#')
            n_dss_dmls = len(dss_dml)
        except Exception:
            n_dss_dmls = 0

    # Generate sample table rows
    case_rows = ''.join(f'<tr><td>{html.escape(s)}</td><td class="case">Case</td></tr>' for s in case_samples)
    control_rows = ''.join(f'<tr><td>{html.escape(s)}</td><td class="control">Control</td></tr>' for s in control_samples)

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
            {f'''<div class="summary-item">
                <div class="summary-value">{n_sig_dmrs}</div>
                <div class="summary-label">Significant DMRs (of {n_total_dmrs})</div>
            </div>''' if has_dmr_stats else ''}
            {f'''<div class="summary-item">
                <div class="summary-value">{n_sig_regions}</div>
                <div class="summary-label">Significant Regions (of {n_total_regions})</div>
            </div>''' if has_region_stats else ''}
            {f'''<div class="summary-item">
                <div class="summary-value">{n_dss_dmrs}</div>
                <div class="summary-label">DSS DMRs</div>
            </div>''' if has_dss_dmr else ''}
            {f'''<div class="summary-item">
                <div class="summary-value">{n_dss_dmls}</div>
                <div class="summary-label">DSS DMLs</div>
            </div>''' if has_dss_dml else ''}
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
        
        {f'''<h2>Differentially Methylated Regions</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.volcano)}" alt="Volcano Plot">
            <div class="figure-caption">Volcano plot showing methylation differences vs. statistical significance.</div>
        </div>''' if args.volcano and os.path.exists(args.volcano) else ''}
        
        {f'''<h2>Methylation Heatmap (Predefined Regions)</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.heatmap)}" alt="Methylation Heatmap">
            <div class="figure-caption">Heatmap of methylation levels at significant predefined regions across all samples.</div>
        </div>''' if args.heatmap and os.path.exists(args.heatmap) else ''}
        
        {f'''<h2>Significant DMR Heatmap</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.dmr_heatmap)}" alt="Significant DMR Heatmap">
            <div class="figure-caption">Heatmap of methylation levels at significantly differentially methylated regions (from signature analysis).</div>
        </div>''' if args.dmr_heatmap and os.path.exists(args.dmr_heatmap) else ''}
        
        {f'''<h2>DSS Analysis Results</h2>
        <p>DSS (Dispersion Shrinkage for Sequencing) identified {n_dss_dmrs} differentially methylated regions (DMRs) and {n_dss_dmls} differentially methylated loci (DMLs).</p>
        <p>See <a href="../igv_dss_dmrs.html">IGV DSS DMR Report</a> and <a href="../igv_dss_dmls.html">IGV DSS DML Report</a> for interactive visualization.</p>
        ''' if has_dss_dmr or has_dss_dml else ''}
        
        {f'''<h2>DSS DMR Heatmap</h2>
        <div class="figure">
            <img src="data:image/png;base64,{img_to_base64(args.dss_dmr_heatmap)}" alt="DSS DMR Heatmap">
            <div class="figure-caption">Heatmap of methylation levels at all DSS-called DMRs across samples.</div>
        </div>''' if args.dss_dmr_heatmap and os.path.exists(args.dss_dmr_heatmap) else ''}
        
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
