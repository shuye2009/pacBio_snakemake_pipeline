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
        top_n=TOP_N,
        igv_top_n=IGV_TOP_N,
        script=os.path.join(SCRIPTS_DIR, "plot_methylation_heatmap.py"),
    resources:
        mem_mb=64000,
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
            --top-n {params.top_n} \
            --igv-top-n {params.igv_top_n} \
            --output-png {output.heatmap} \
            --output-pdf {output.heatmap_pdf} \
            --output-tsv {output.stats_tsv} \
            --output-bed {output.sig_bed} \
            --output-igv-tsv {output.sig_tsv} \
            --output-zscore-dist-png {output.zscore_dist_png} \
            --output-zscore-dist-pdf {output.zscore_dist_pdf} \
            2>&1 | tee {log}
        """


if ENOUGH_SAMPLES:

    # =========================================================================
    # Plot_dmr_volcano: Volcano plot of differentially methylated regions
    #
    # Shows the relationship between methylation difference (delta) and
    # statistical significance (z-score) for each region.
    # =========================================================================
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


    # =========================================================================
    # Plot_significant_dmr_heatmap: Heatmap of significant DMRs only
    #
    # Shows methylation levels at significantly differentially methylated regions.
    # DMR coordinates from dmr_stats_with_pvalues.tsv (significant=True), methylation from combined.bed.gz.
    # =========================================================================
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
# Plot_dss_dmr_heatmap: Heatmap of DSS-called DMRs
#
# Shows methylation levels at all DSS DMRs without filtering.
# Reads DSS DMR BED file directly.
# =============================================================================
rule Plot_dss_dmr_heatmap:
    input:
        dmr_bed=DSS_BASE + f"/dmr_results.top{TOP_N}.bed",
        beds=expand(config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz", sample=ALL_SAMPLES),
    output:
        heatmap=VIS_BASE + "/dss_dmr_heatmap.png",
        heatmap_pdf=VIS_BASE + "/dss_dmr_heatmap.pdf",
    params:
        case_samples=",".join(config["samples"]["case"]),
        control_samples=",".join(config["samples"]["control"]),
        bed_dir=config["directory"]["output"] + "/pb_cpg_tools",
        script=os.path.join(SCRIPTS_DIR, "plot_dss_dmr_heatmap.py"),
    log:
        config["directory"]["output"] + "/logs/visualization/dss_dmr_heatmap.log",
    shell:
        """
        mkdir -p $(dirname {output.heatmap})
        mkdir -p $(dirname {log})
        
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --dmr-bed {input.dmr_bed} \
            --bed-dir {params.bed_dir} \
            --case-samples {params.case_samples} \
            --control-samples {params.control_samples} \
            --output-png {output.heatmap} \
            --output-pdf {output.heatmap_pdf} \
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
# Uses haplotagged BAMs when phasing is enabled, aligned BAMs otherwise.
# =============================================================================
rule IGV_reports_regions:
    input:
        tsv=VIS_DIR + "/significant_regions.tsv",
        bed=VIS_DIR + "/significant_regions.bed",
        fasta=config["genome"]["fasta"],
        gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
        gtf_index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
        bams=expand(
            config["directory"]["output"] + ("/phased/{sample}.haplotagged.bam" if PHASING_ENABLED else "/aligned/{sample}.aligned.bam"),
            sample=ALL_SAMPLES
        ),
        bais=expand(
            config["directory"]["output"] + ("/phased/{sample}.haplotagged.bam.bai" if PHASING_ENABLED else "/aligned/{sample}.aligned.bam.bai"),
            sample=ALL_SAMPLES
        ),
    output:
        report=VIS_DIR + "/igv_significant_regions.html",
        track_config=VIS_DIR + "/igv_regions_track_config.json",
    params:
        bam_args=lambda wildcards, input: " ".join(input.bams),
        script=os.path.join(SCRIPTS_DIR, "generate_igv_track_config.py"),
        phasing_flag="--phasing-enabled" if PHASING_ENABLED else "",
    log:
        config["directory"]["output"] + "/logs/visualization/igv_regions.log",
    shell:
        """
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {log})
        
        # Generate track config for BAM files with haplotype grouping and methylation coloring
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --bams {params.bam_args} \
            --output {output.track_config} \
            --gtf {input.gtf} --gtf-index {input.gtf_index} \
            {params.phasing_flag}
        
        module load igv-reports
        
        create_report {input.tsv} \
            --fasta {input.fasta} \
            --sequence 1 --begin 2 --end 3 \
            --info-columns NAME DELTA ADJ_PVALUE \
            --tracks {input.bed} \
            --track-config {output.track_config} \
            --output {output.report} \
            2>&1 | tee {log}
        """


if ENOUGH_SAMPLES:

    # =========================================================================
    # Filter_dmrs_for_igv: Filter significant DMRs to top N by adjusted p-value
    #
    # Reduces the number of DMRs for IGV reporting to keep reports manageable.
    # =========================================================================
    rule Filter_dmrs_for_igv:
        input:
            tsv=VIS_BASE + "/significant_dmrs.tsv",
            bed=VIS_BASE + "/significant_dmrs.bed",
        output:
            tsv=VIS_BASE + f"/significant_dmrs.top{TOP_N}.tsv",
            bed=VIS_BASE + f"/significant_dmrs.top{TOP_N}.bed",
        params:
            script=os.path.join(SCRIPTS_DIR, "filter_dmrs_for_igv.py"),
            top_n=TOP_N,
        log:
            config["directory"]["output"] + "/logs/visualization/filter_dmrs_for_igv.log",
        shell:
            """
            mkdir -p $(dirname {output.tsv})
            mkdir -p $(dirname {log})

            /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
                --input-tsv {input.tsv} \
                --input-bed {input.bed} \
                --top-n {params.top_n} \
                --output-tsv {output.tsv} \
                --output-bed {output.bed} \
                2>&1 | tee {log}
            """


    rule IGV_reports_dmrs:
        input:
            tsv=VIS_BASE + f"/significant_dmrs.top{TOP_N}.tsv",
            bed=VIS_BASE + f"/significant_dmrs.top{TOP_N}.bed",
            fasta=config["genome"]["fasta"],
            gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
            gtf_index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
            bams=expand(
                config["directory"]["output"] + ("/phased/{sample}.haplotagged.bam" if PHASING_ENABLED else "/aligned/{sample}.aligned.bam"),
                sample=ALL_SAMPLES
            ),
            bais=expand(
                config["directory"]["output"] + ("/phased/{sample}.haplotagged.bam.bai" if PHASING_ENABLED else "/aligned/{sample}.aligned.bam.bai"),
                sample=ALL_SAMPLES
            ),
        output:
            report=VIS_BASE + "/igv_significant_dmrs.html",
            track_config=VIS_BASE + "/igv_dmrs_track_config.json",
        params:
            bam_args=lambda wildcards, input: " ".join(input.bams),
            script=os.path.join(SCRIPTS_DIR, "generate_igv_track_config.py"),
            phasing_flag="--phasing-enabled" if PHASING_ENABLED else "",
        log:
            config["directory"]["output"] + "/logs/visualization/igv_dmrs.log",
        shell:
            """
            mkdir -p $(dirname {output.report})
            mkdir -p $(dirname {log})
            
            # Generate track config for BAM files with haplotype grouping and methylation coloring
            /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
                --bams {params.bam_args} \
                --output {output.track_config} \
                --gtf {input.gtf} --gtf-index {input.gtf_index} \
                {params.phasing_flag}
            
            module load igv-reports
            
            create_report {input.tsv} \
                --fasta {input.fasta} \
                --sequence 1 --begin 2 --end 3 \
                --info-columns NAME DELTA ADJ_PVALUE \
                --tracks {input.bed} \
                --track-config {output.track_config} \
                --output {output.report} \
                2>&1 | tee {log}
            """


# =============================================================================
# Filter_dss_dmrs_for_igv: Filter DSS DMRs to top N by abs(methylation diff)
#
# Produces two filtered sets: TOP_N for motif enrichment/annotation,
# IGV_TOP_N for IGV browser reports.
# =============================================================================
rule Filter_dss_dmrs_for_igv:
    input:
        tsv=DSS_BASE + "/dmr_results.tsv",
        bed=DSS_BASE + "/dmr_results.bed",
    output:
        tsv=DSS_BASE + f"/dmr_results.top{TOP_N}.tsv",
        bed=DSS_BASE + f"/dmr_results.top{TOP_N}.bed",
        tsv_igv=DSS_BASE + f"/dmr_results.top{IGV_TOP_N}.tsv",
        bed_igv=DSS_BASE + f"/dmr_results.top{IGV_TOP_N}.bed",
    params:
        script=os.path.join(SCRIPTS_DIR, "filter_dss_dmrs_for_igv.py"),
        top_n=TOP_N,
        igv_top_n=IGV_TOP_N,
    log:
        config["directory"]["output"] + "/logs/visualization/filter_dss_dmrs_for_igv.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --input-tsv {input.tsv} \
            --input-bed {input.bed} \
            --top-n {params.top_n} \
            --output-tsv {output.tsv} \
            --output-bed {output.bed}

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --input-tsv {input.tsv} \
            --input-bed {input.bed} \
            --top-n {params.igv_top_n} \
            --output-tsv {output.tsv_igv} \
            --output-bed {output.bed_igv} \
            2>&1 | tee {log}
        """


# =============================================================================
# IGV_reports_dss_dmrs: Generate IGV report for DSS DMR results
#
# Creates interactive HTML report showing differentially methylated regions
# called by DSS (Dispersion Shrinkage for Sequencing).
# Uses haplotagged BAMs when phasing is enabled, aligned BAMs otherwise.
# =============================================================================
rule IGV_reports_dss_dmrs:
    input:
        tsv=DSS_BASE + f"/dmr_results.top{IGV_TOP_N}.tsv",
        bed=DSS_BASE + f"/dmr_results.top{IGV_TOP_N}.bed",
        fasta=config["genome"]["fasta"],
        gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
        gtf_index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
        bams=expand(
            config["directory"]["output"] + ("/phased/{sample}.haplotagged.bam" if PHASING_ENABLED else "/aligned/{sample}.aligned.bam"),
            sample=ALL_SAMPLES
        ),
        bais=expand(
            config["directory"]["output"] + ("/phased/{sample}.haplotagged.bam.bai" if PHASING_ENABLED else "/aligned/{sample}.aligned.bam.bai"),
            sample=ALL_SAMPLES
        ),
    output:
        report=VIS_BASE + "/igv_dss_dmrs.html",
        track_config=VIS_BASE + "/igv_dss_dmrs_track_config.json",
    params:
        bam_args=lambda wildcards, input: " ".join(input.bams),
        script=os.path.join(SCRIPTS_DIR, "generate_igv_track_config.py"),
        phasing_flag="--phasing-enabled" if PHASING_ENABLED else "",
    log:
        config["directory"]["output"] + "/logs/visualization/igv_dss_dmrs.log",
    shell:
        """
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {log})

        # Fix scientific notation in BED coordinates for igv-reports compatibility
        bed_fixed=$(dirname {output.report})/dmr_results_fixed.bed
        awk 'BEGIN{{OFS="\t"}} {{if(NF>=3){{$2=sprintf("%.0f",$2);$3=sprintf("%.0f",$3)}};print}}' {input.bed} > "$bed_fixed"

        # Generate track config for BAM files with haplotype grouping and methylation coloring
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --bams {params.bam_args} \
            --output {output.track_config} \
            --gtf {input.gtf} --gtf-index {input.gtf_index} \
            {params.phasing_flag}

        module load igv-reports

        create_report {input.tsv} \
            --fasta {input.fasta} \
            --sequence 1 --begin 2 --end 3 \
            --info-columns areaStat diff.Methy nCG length \
            --tracks "$bed_fixed" \
            --track-config {output.track_config} \
            --output {output.report} \
            2>&1 | tee {log}

        rm -f "$bed_fixed"
        """


# =============================================================================
# Filter_dml_for_igv: Extract top N most significant DMLs for IGV reporting
#
# Full DML results can contain hundreds of thousands of loci, which
# overwhelms igv-reports. Produces two filtered sets: TOP_N and IGV_TOP_N.
# =============================================================================
rule Filter_dml_for_igv:
    input:
        tsv=DSS_BASE + "/dml_results.tsv",
    output:
        tsv=DSS_BASE + f"/dml_results.top{TOP_N}.tsv",
        bed=DSS_BASE + f"/dml_results.top{TOP_N}.bed",
        tsv_igv=DSS_BASE + f"/dml_results.top{IGV_TOP_N}.tsv",
        bed_igv=DSS_BASE + f"/dml_results.top{IGV_TOP_N}.bed",
    params:
        script=os.path.join(SCRIPTS_DIR, "filter_dml_for_igv.py"),
        top_n=TOP_N,
        igv_top_n=IGV_TOP_N,
    log:
        config["directory"]["output"] + "/logs/visualization/filter_dml_for_igv.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --dml-tsv {input.tsv} \
            --top-n {params.top_n} \
            --output-tsv {output.tsv} \
            --output-bed {output.bed}

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --dml-tsv {input.tsv} \
            --top-n {params.igv_top_n} \
            --output-tsv {output.tsv_igv} \
            --output-bed {output.bed_igv} \
            2>&1 | tee {log}
        """


# =============================================================================
# IGV_reports_dss_dmls: Generate IGV report for DSS DML results
#
# Creates interactive HTML report showing the top differentially methylated
# loci (single CpG sites) called by DSS.
# Uses haplotagged BAMs when phasing is enabled, aligned BAMs otherwise.
# =============================================================================
rule IGV_reports_dss_dmls:
    input:
        tsv=DSS_BASE + f"/dml_results.top{IGV_TOP_N}.tsv",
        bed=DSS_BASE + f"/dml_results.top{IGV_TOP_N}.bed",
        fasta=config["genome"]["fasta"],
        gtf=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz",
        gtf_index=config["directory"]["output"] + "/visualization/genes.sorted.gtf.gz.tbi",
        bams=expand(
            config["directory"]["output"] + ("/phased/{sample}.haplotagged.bam" if PHASING_ENABLED else "/aligned/{sample}.aligned.bam"),
            sample=ALL_SAMPLES
        ),
        bais=expand(
            config["directory"]["output"] + ("/phased/{sample}.haplotagged.bam.bai" if PHASING_ENABLED else "/aligned/{sample}.aligned.bam.bai"),
            sample=ALL_SAMPLES
        ),
    output:
        report=VIS_BASE + "/igv_dss_dmls.html",
        track_config=VIS_BASE + "/igv_dss_dmls_track_config.json",
    params:
        bam_args=lambda wildcards, input: " ".join(input.bams),
        script=os.path.join(SCRIPTS_DIR, "generate_igv_track_config.py"),
        phasing_flag="--phasing-enabled" if PHASING_ENABLED else "",
    log:
        config["directory"]["output"] + "/logs/visualization/igv_dss_dmls.log",
    shell:
        """
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {log})

        # Fix scientific notation in BED coordinates for igv-reports compatibility
        bed_fixed=$(dirname {output.report})/dml_results_fixed.bed
        awk 'BEGIN{{OFS="\t"}} {{if(NF>=3){{$2=sprintf("%.0f",$2);$3=sprintf("%.0f",$3)}};print}}' {input.bed} > "$bed_fixed"

        # Generate track config for BAM files with haplotype grouping and methylation coloring
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --bams {params.bam_args} \
            --output {output.track_config} \
            --gtf {input.gtf} --gtf-index {input.gtf_index} \
            {params.phasing_flag}

        module load igv-reports

        create_report {input.tsv} \
            --fasta {input.fasta} \
            --sequence 1 --begin 2 --end 2 \
            --info-columns mu1 mu2 diff pval fdr \
            --tracks "$bed_fixed" \
            --track-config {output.track_config} \
            --output {output.report} \
            2>&1 | tee {log}

        rm -f "$bed_fixed"
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
            igv_top_n=IGV_TOP_N,
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
                --igv-top-n {params.igv_top_n} \
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
                --output {output.track_config} \
                --gtf {input.gtf} --gtf-index {input.gtf_index}
            
            module load igv-reports
            
            create_report {input.tsv} \
                --fasta {input.fasta} \
                --sequence 1 --begin 2 --end 3 \
                --info-columns NAME AVG_DELTA AVG_METHYL \
                --tracks {input.bed} \
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
            igv_top_n=IGV_TOP_N,
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
                --igv-top-n {params.igv_top_n} \
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
                --output {output.track_config} \
                --gtf {input.gtf} --gtf-index {input.gtf_index}
            
            module load igv-reports
            
            create_report {input.tsv} \
                --fasta {input.fasta} \
                --sequence 1 --begin 2 --end 3 \
                --info-columns NAME CASE_DELTA CONTROL_DELTA DIFF_DELTA \
                --tracks {input.bed} \
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
        pca=VIS_DIR + "/region_sample_pca.png",
        dist=config["directory"]["output"] + "/global/methylation_distribution.png",
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
        dss_dmr=DSS_BASE + "/dmr_results.tsv",
        dss_dml=DSS_BASE + "/dml_results.tsv",
        dss_dmr_heatmap=VIS_BASE + "/dss_dmr_heatmap.png",
    output:
        report=VIS_DIR + "/methylation_report.html",
    params:
        case_samples=",".join(config["samples"]["case"]),
        control_samples=",".join(config["samples"]["control"]),
        script=os.path.join(SCRIPTS_DIR, "generate_report.py"),
        samples=" ".join(ALL_SAMPLES),
    log:
        config["directory"]["output"] + "/logs/visualization/report.log",
    run:
        import os, subprocess

        os.makedirs(os.path.dirname(output.report), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        cmd = [
            "/cluster/home/t128737uhn/miniconda3/bin/python", params.script,
            "--pca", input.pca,
            "--distribution", input.dist,
            "--global-heatmap", input.global_heatmap,
            "--global-pca", input.global_pca,
            "--global-density-composite", input.global_density_composite,
            "--global-density-per-sample", " ".join(input.global_density_per_sample),
            "--global-density-per-sample-samples", params.samples,
            "--global-density-per-chrom", " ".join(input.global_density_per_chrom),
            "--global-density-per-chrom-samples", params.samples,
            "--case-samples", params.case_samples,
            "--control-samples", params.control_samples,
            "--dss-dmr", input.dss_dmr,
            "--dss-dml", input.dss_dml,
            "--dss-dmr-heatmap", input.dss_dmr_heatmap,
            "--output", output.report,
        ]

        cmd.extend([
            "--heatmap", VIS_DIR + "/methylation_heatmap.png",
        ])
        if ENOUGH_SAMPLES:
            cmd.extend([
                "--dmr-heatmap", VIS_BASE + "/significant_dmr_heatmap.png",
                "--volcano", VIS_BASE + "/dmr_volcano.png",
                "--signature-stats", METHBAT_BASE + "/signature.signature_stats.tsv",
                "--dmr-stats", VIS_BASE + "/dmr_stats_with_pvalues.tsv",
                "--region-stats", VIS_DIR + "/region_cohort_comparison_with_pvalues.tsv",
            ])

        with open(log[0], "w") as log_f:
            subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, check=True)
