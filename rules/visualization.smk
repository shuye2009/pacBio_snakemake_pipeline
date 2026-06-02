# =============================================================================
# Visualization rules for PacBio methylation analysis results
# =============================================================================

import os

ALL_SAMPLES = config["samples"]["case"] + config["samples"]["control"]
SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")


# =============================================================================
# Index_gtf: Sort and tabix-index GTF file for igv-reports
#
# Large GTF files need to be sorted and indexed for efficient access.
# =============================================================================
rule Index_gtf:
    input:
        gtf=config["genome"]["gtf"],
    output:
        sorted_gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
        index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
    log:
        config["directory"]["output"] + "/logs/visualization/index_gtf.log",
    shell:
        """
        mkdir -p $(dirname {output.sorted_gtf})
        mkdir -p $(dirname {log})
        
        module load samtools/1.20
        
        # Sort GTF by chromosome and position, compress with bgzip, and index with tabix
        (grep "^#" {input.gtf} || true; grep -v "^#" {input.gtf} | sort -k1,1 -k4,4n) | \
            bgzip -c > {output.sorted_gtf}
        
        tabix -p gff {output.sorted_gtf}
        
        echo "Sorted and indexed GTF file" 2>&1 | tee {log}
        """


# =============================================================================
# Plot_methylation_heatmap: Generate heatmap of methylation across samples
#
# Creates a heatmap showing methylation levels at significant regions from
# cohort_comparison.tsv, filtered by min_delta and adjusted p-value thresholds.
# =============================================================================
rule Plot_methylation_heatmap_region:
    input:
        cohort_comparison=METHBAT_DIR + "/region_cohort_comparison.tsv",
        profiles=expand(METHBAT_DIR + "/profiles_region/{sample}.region.profile.tsv", sample=ALL_SAMPLES),
    output:
        heatmap=VIS_DIR + "/methylation_heatmap.png",
        heatmap_pdf=VIS_DIR + "/methylation_heatmap.pdf",
        stats_tsv=VIS_DIR + "/region_cohort_comparison_with_pvalues.tsv",
        sig_bed=VIS_DIR + "/significant_regions.bed",
        sig_tsv=VIS_DIR + "/significant_regions.tsv",
        zscore_dist_png=VIS_DIR + "/zscore_distribution_regions.png",
        zscore_dist_pdf=VIS_DIR + "/zscore_distribution_regions.pdf",
    params:
        case_samples=",".join(config["samples"]["case"]),
        control_samples=",".join(config["samples"]["control"]),
        profile_dir=METHBAT_DIR + "/profiles_region",
        min_delta=config["methbat"]["min_delta"],
        pvalue_cutoff=config["methbat"]["pvalue_cutoff"],
        script=os.path.join(SCRIPTS_DIR, "plot_methylation_heatmap.py"),
    log:
        config["directory"]["output"] + "/logs/visualization/heatmap.log",
    shell:
        """
        mkdir -p $(dirname {output.heatmap})
        mkdir -p $(dirname {log})
        
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --cohort-comparison {input.cohort_comparison} \
            --profile-dir {params.profile_dir} \
            --case-samples {params.case_samples} \
            --control-samples {params.control_samples} \
            --min-delta {params.min_delta} \
            --pvalue-cutoff {params.pvalue_cutoff} \
            --output-png {output.heatmap} \
            --output-pdf {output.heatmap_pdf} \
            --output-tsv {output.stats_tsv} \
            --output-bed {output.sig_bed} \
            --output-igv-tsv {output.sig_tsv} \
            --output-zscore-dist-png {output.zscore_dist_png} \
            --output-zscore-dist-pdf {output.zscore_dist_pdf} \
            2>&1 | tee {log}
        """


# =============================================================================
# Plot_dmr_volcano: Volcano plot of differentially methylated regions
#
# Shows the relationship between methylation difference (delta) and
# statistical significance (z-score) for each region.
# =============================================================================
rule Plot_dmr_volcano:
    input:
        stats=METHBAT_BASE + "/signature.signature_stats.tsv",
    output:
        volcano=VIS_BASE + "/dmr_volcano.png",
        volcano_pdf=VIS_BASE + "/dmr_volcano.pdf",
        stats_tsv=VIS_BASE + "/dmr_stats_with_pvalues.tsv",
        zscore_dist_png=VIS_BASE + "/zscore_distribution_dmrs.png",
        zscore_dist_pdf=VIS_BASE + "/zscore_distribution_dmrs.pdf",
    params:
        min_delta=config["methbat"]["min_delta"],
        pvalue_cutoff=config["methbat"]["pvalue_cutoff"],
        script=os.path.join(SCRIPTS_DIR, "plot_dmr_volcano.py"),
    log:
        config["directory"]["output"] + "/logs/visualization/volcano.log",
    shell:
        """
        mkdir -p $(dirname {output.volcano})
        mkdir -p $(dirname {log})
        
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --stats {input.stats} \
            --min-delta {params.min_delta} \
            --pvalue-cutoff {params.pvalue_cutoff} \
            --output-png {output.volcano} \
            --output-pdf {output.volcano_pdf} \
            --output-tsv {output.stats_tsv} \
            --output-zscore-dist-png {output.zscore_dist_png} \
            --output-zscore-dist-pdf {output.zscore_dist_pdf} \
            2>&1 | tee {log}
        """


# =============================================================================
# Plot_significant_dmr_heatmap: Heatmap of significant DMRs only
#
# Shows methylation levels at significantly differentially methylated regions.
# DMR coordinates from dmr_stats_with_pvalues.tsv (significant=True), methylation from combined.bed.gz.
# =============================================================================
rule Plot_significant_dmr_heatmap:
    input:
        dmr_stats=VIS_BASE + "/dmr_stats_with_pvalues.tsv",
        beds=expand(config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz", sample=ALL_SAMPLES),
    output:
        heatmap=VIS_BASE + "/significant_dmr_heatmap.png",
        heatmap_pdf=VIS_BASE + "/significant_dmr_heatmap.pdf",
        sig_bed=VIS_BASE + "/significant_dmrs.bed",
        sig_tsv=VIS_BASE + "/significant_dmrs.tsv",
    params:
        case_samples=",".join(config["samples"]["case"]),
        control_samples=",".join(config["samples"]["control"]),
        bed_dir=config["directory"]["output"] + "/pb_cpg_tools",
        script=os.path.join(SCRIPTS_DIR, "plot_significant_dmr_heatmap.py"),
    log:
        config["directory"]["output"] + "/logs/visualization/significant_dmr_heatmap.log",
    shell:
        """
        mkdir -p $(dirname {output.heatmap})
        mkdir -p $(dirname {log})
        
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --dmr-stats {input.dmr_stats} \
            --bed-dir {params.bed_dir} \
            --case-samples {params.case_samples} \
            --control-samples {params.control_samples} \
            --output-png {output.heatmap} \
            --output-pdf {output.heatmap_pdf} \
            --output-bed {output.sig_bed} \
            --output-tsv {output.sig_tsv} \
            2>&1 | tee {log}
        """


# =============================================================================
# Plot_sample_pca: PCA plot of samples based on methylation profiles
#
# Visualizes sample clustering based on methylation patterns at signature
# regions. Useful for identifying outliers and batch effects.
# =============================================================================
rule Plot_sample_pca_region:
    input:
        profiles=expand(METHBAT_DIR + "/profiles_region/{sample}.region.profile.tsv", sample=ALL_SAMPLES),
    output:
        pca=VIS_DIR + "/region_sample_pca.png",
        pca_pdf=VIS_DIR + "/region_sample_pca.pdf",
    params:
        case_samples=",".join(config["samples"]["case"]),
        control_samples=",".join(config["samples"]["control"]),
        profile_dir=METHBAT_DIR + "/profiles_region",
        script=os.path.join(SCRIPTS_DIR, "plot_sample_pca.py"),
    log:
        config["directory"]["output"] + "/logs/visualization/pca.log",
    shell:
        """
        mkdir -p $(dirname {output.pca})
        mkdir -p $(dirname {log})
        
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --profile-dir {params.profile_dir} \
            --case-samples {params.case_samples} \
            --control-samples {params.control_samples} \
            --output-png {output.pca} \
            --output-pdf {output.pca_pdf} \
            2>&1 | tee {log}
        """



# =============================================================================
# IGV_reports: Generate interactive IGV HTML reports for significant regions
#
# Creates browsable HTML reports with IGV.js for visualizing significant
# regions and DMRs with haplotagged BAM tracks.
# =============================================================================
rule IGV_reports_regions:
    input:
        tsv=VIS_DIR + "/significant_regions.tsv",
        bed=VIS_DIR + "/significant_regions.bed",
        fasta=config["genome"]["fasta"],
        gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
        gtf_index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
        bams=expand(config["directory"]["output"] + "/phased/{sample}.haplotagged.bam", sample=ALL_SAMPLES),
        bais=expand(config["directory"]["output"] + "/phased/{sample}.haplotagged.bam.bai", sample=ALL_SAMPLES),
    output:
        report=VIS_DIR + "/igv_significant_regions.html",
        track_config=VIS_DIR + "/igv_regions_track_config.json",
    params:
        bam_args=lambda wildcards, input: " ".join(input.bams),
        script=os.path.join(SCRIPTS_DIR, "generate_igv_track_config.py"),
    log:
        config["directory"]["output"] + "/logs/visualization/igv_regions.log",
    shell:
        """
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {log})
        
        # Generate track config for BAM files with haplotype grouping and methylation coloring
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --bams {params.bam_args} \
            --output {output.track_config}
        
        module load igv-reports
        
        create_report {input.tsv} \
            --fasta {input.fasta} \
            --sequence 1 --begin 2 --end 3 \
            --info-columns NAME DELTA ADJ_PVALUE \
            --tracks {input.bed} {input.gtf} \
            --track-config {output.track_config} \
            --output {output.report} \
            2>&1 | tee {log}
        """


rule IGV_reports_dmrs:
    input:
        tsv=VIS_BASE + "/significant_dmrs.tsv",
        bed=VIS_BASE + "/significant_dmrs.bed",
        fasta=config["genome"]["fasta"],
        gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
        gtf_index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
        bams=expand(config["directory"]["output"] + "/phased/{sample}.haplotagged.bam", sample=ALL_SAMPLES),
        bais=expand(config["directory"]["output"] + "/phased/{sample}.haplotagged.bam.bai", sample=ALL_SAMPLES),
    output:
        report=VIS_BASE + "/igv_significant_dmrs.html",
        track_config=VIS_BASE + "/igv_dmrs_track_config.json",
    params:
        bam_args=lambda wildcards, input: " ".join(input.bams),
        script=os.path.join(SCRIPTS_DIR, "generate_igv_track_config.py"),
    log:
        config["directory"]["output"] + "/logs/visualization/igv_dmrs.log",
    shell:
        """
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {log})
        
        # Generate track config for BAM files with haplotype grouping and methylation coloring
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --bams {params.bam_args} \
            --output {output.track_config}
        
        module load igv-reports
        
        create_report {input.tsv} \
            --fasta {input.fasta} \
            --sequence 1 --begin 2 --end 3 \
            --info-columns NAME DELTA ADJ_PVALUE \
            --tracks {input.bed} {input.gtf} \
            --track-config {output.track_config} \
            --output {output.report} \
            2>&1 | tee {log}
        """


# =============================================================================
# ASM Visualization Rules (when phasing is enabled)
# =============================================================================
if PHASING_ENABLED:
    # =========================================================================
    # Extract_significant_ASM: Extract significant ASM regions for visualization
    # =========================================================================
    rule Extract_significant_ASM:
        input:
            asm_comparison=METHBAT_BASE + "/asm_haplotype_comparison.tsv",
        output:
            bed=VIS_BASE + "/significant_asm.bed",
            tsv=VIS_BASE + "/significant_asm.tsv",
        params:
            min_delta=config["methbat"]["min_delta"],
            script=os.path.join(SCRIPTS_DIR, "extract_significant_asm.py"),
        log:
            config["directory"]["output"] + "/logs/visualization/extract_significant_asm.log",
        shell:
            """
            mkdir -p $(dirname {output.bed})
            mkdir -p $(dirname {log})
            
            /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
                --asm-comparison {input.asm_comparison} \
                --min-delta {params.min_delta} \
                --output-bed {output.bed} \
                --output-tsv {output.tsv} \
                2>&1 | tee {log}
            """


    # =========================================================================
    # IGV_reports_ASM: Generate IGV report for differential ASM regions
    #
    # Creates interactive HTML report showing regions with significant
    # allele-specific methylation differences between haplotypes.
    # =========================================================================
    rule IGV_reports_ASM:
        input:
            tsv=VIS_BASE + "/significant_asm.tsv",
            bed=VIS_BASE + "/significant_asm.bed",
            fasta=config["genome"]["fasta"],
            gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
            gtf_index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
            bams=expand(config["directory"]["output"] + "/phased/{sample}.haplotagged.bam", sample=ALL_SAMPLES),
            bais=expand(config["directory"]["output"] + "/phased/{sample}.haplotagged.bam.bai", sample=ALL_SAMPLES),
        output:
            report=VIS_BASE + "/igv_significant_asm.html",
            track_config=VIS_BASE + "/igv_significant_asm_config.json",
        params:
            bam_args=lambda wildcards, input: " ".join(input.bams),
            script=os.path.join(SCRIPTS_DIR, "generate_igv_track_config.py"),
        log:
            config["directory"]["output"] + "/logs/visualization/igv_significant_asm.log",
        shell:
            """
            mkdir -p $(dirname {output.report})
            mkdir -p $(dirname {log})
            
            # Generate track config for BAM files with haplotype grouping and methylation coloring
            /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
                --bams {params.bam_args} \
                --output {output.track_config}
            
            module load igv-reports
            
            create_report {input.tsv} \
                --fasta {input.fasta} \
                --sequence 1 --begin 2 --end 3 \
                --info-columns NAME AVG_DELTA AVG_METHYL \
                --tracks {input.bed} {input.gtf} \
                --track-config {output.track_config} \
                --output {output.report} \
                2>&1 | tee {log}
            """


    # =========================================================================
    # Extract_differential_ASM: Extract regions with differential ASM between groups
    #
    # Filters ASM group comparison (case vs control) to find regions where
    # radiation changed allele-specific methylation patterns.
    # =========================================================================
    rule Extract_differential_ASM:
        input:
            asm_comparison=METHBAT_BASE + "/asm_cohort.profile.tsv",
        output:
            bed=VIS_BASE + "/differential_asm.bed",
            tsv=VIS_BASE + "/differential_asm.tsv",
        params:
            min_delta=config["methbat"]["min_delta"],
            pvalue_cutoff=config["methbat"].get("pvalue_cutoff", 0.05),
            script=os.path.join(SCRIPTS_DIR, "extract_differential_asm.py"),
        log:
            config["directory"]["output"] + "/logs/visualization/extract_differential_asm.log",
        shell:
            """
            mkdir -p $(dirname {output.bed})
            mkdir -p $(dirname {log})
            
            /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
                --asm-comparison {input.asm_comparison} \
                --min-delta {params.min_delta} \
                --pvalue-cutoff {params.pvalue_cutoff} \
                --output-bed {output.bed} \
                --output-tsv {output.tsv} \
                2>&1 | tee {log}
            """


    # =========================================================================
    # IGV_reports_differential_ASM: IGV report for radiation-induced ASM changes
    #
    # Creates interactive HTML report showing regions where allele-specific
    # methylation differs between irradiated (case) and control groups.
    # =========================================================================
    rule IGV_reports_differential_ASM:
        input:
            tsv=VIS_BASE + "/differential_asm.tsv",
            bed=VIS_BASE + "/differential_asm.bed",
            fasta=config["genome"]["fasta"],
            gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
            gtf_index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
            bams=expand(config["directory"]["output"] + "/phased/{sample}.haplotagged.bam", sample=ALL_SAMPLES),
            bais=expand(config["directory"]["output"] + "/phased/{sample}.haplotagged.bam.bai", sample=ALL_SAMPLES),
        output:
            report=VIS_BASE + "/igv_differential_asm.html",
            track_config=VIS_BASE + "/igv_differential_asm_track_config.json",
        params:
            bam_args=lambda wildcards, input: " ".join(input.bams),
            script=os.path.join(SCRIPTS_DIR, "generate_igv_track_config.py"),
        log:
            config["directory"]["output"] + "/logs/visualization/igv_differential_asm.log",
        shell:
            """
            mkdir -p $(dirname {output.report})
            mkdir -p $(dirname {log})
            
            # Generate track config for BAM files with haplotype grouping and methylation coloring
            /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
                --bams {params.bam_args} \
                --output {output.track_config}
            
            module load igv-reports
            
            create_report {input.tsv} \
                --fasta {input.fasta} \
                --sequence 1 --begin 2 --end 3 \
                --info-columns NAME CASE_DELTA CONTROL_DELTA DIFF_DELTA \
                --tracks {input.bed} {input.gtf} \
                --track-config {output.track_config} \
                --output {output.report} \
                2>&1 | tee {log}
            """


# =============================================================================
# Generate_report: Create HTML summary report
#
# Combines all visualizations into a single HTML report with statistics.
# =============================================================================
rule Generate_report:
    input:
        heatmap=VIS_DIR + "/methylation_heatmap.png",
        dmr_heatmap=VIS_BASE + "/significant_dmr_heatmap.png",
        volcano=VIS_BASE + "/dmr_volcano.png",
        pca=VIS_DIR + "/region_sample_pca.png",
        dist=VIS_DIR + "/methylation_distribution.png",
        signature_stats=METHBAT_BASE + "/signature.signature_stats.tsv",
        dmr_stats=VIS_BASE + "/dmr_stats_with_pvalues.tsv",
        region_stats=VIS_DIR + "/region_cohort_comparison_with_pvalues.tsv",
        global_heatmap=config["directory"]["output"] + "/global/global_methylation_heatmap.png",
        global_pca=config["directory"]["output"] + "/global/global_methylation_pca.png",
        global_density_composite=config["directory"]["output"] + "/global/density/composite_density.png",
        global_density_per_sample=expand(
            config["directory"]["output"] + "/global/density/{sample}.density.png",
            sample=ALL_SAMPLES
        ),
        global_density_per_chrom=expand(
            config["directory"]["output"] + "/global/density/{sample}.per_chrom_density.png",
            sample=ALL_SAMPLES
        ),
    output:
        report=VIS_DIR + "/methylation_report.html",
    params:
        case_samples=",".join(config["samples"]["case"]),
        control_samples=",".join(config["samples"]["control"]),
        script=os.path.join(SCRIPTS_DIR, "generate_report.py"),
        samples=" ".join(ALL_SAMPLES),
    log:
        config["directory"]["output"] + "/logs/visualization/report.log",
    shell:
        """
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {log})
        
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --heatmap {input.heatmap} \
            --dmr-heatmap {input.dmr_heatmap} \
            --volcano {input.volcano} \
            --pca {input.pca} \
            --distribution {input.dist} \
            --signature-stats {input.signature_stats} \
            --dmr-stats {input.dmr_stats} \
            --region-stats {input.region_stats} \
            --global-heatmap {input.global_heatmap} \
            --global-pca {input.global_pca} \
            --global-density-composite {input.global_density_composite} \
            --global-density-per-sample {input.global_density_per_sample} \
            --global-density-per-sample-samples {params.samples} \
            --global-density-per-chrom {input.global_density_per_chrom} \
            --global-density-per-chrom-samples {params.samples} \
            --case-samples {params.case_samples} \
            --control-samples {params.control_samples} \
            --output {output.report} \
            2>&1 | tee {log}
        """
