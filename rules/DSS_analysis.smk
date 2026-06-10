# =============================================================================
# DSS (Dispersion Shrinkage for Sequencing) differential methylation analysis
#
# Uses the DSS R package to identify differentially methylated loci (DML)
# and differentially methylated regions (DMR) between case and control groups.
#
# DSS performs:
#   1. DMLtest  – dispersion shrinkage + Wald test at each CpG site
#   2. callDML  – filter significant CpG sites by delta and p-value
#   3. callDMR  – merge adjacent significant CpGs into DMRs
# =============================================================================

import os

SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")


# =============================================================================
# DSS_analysis: Run DSS DML test and DMR calling
#
# Loads all sample combined.bed.gz files, builds a BSseq object, runs
# DMLtest with smoothing, calls DMLs and DMRs. Outputs TSV and BED files.
# =============================================================================
rule DSS_analysis:
    input:
        beds=expand(
            config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
            sample=ALL_SAMPLES
        ),
    output:
        dml_tsv=DSS_BASE + "/dml_results.tsv",
        dml_bed=DSS_BASE + "/dml_results.bed",
        dmr_tsv=DSS_BASE + "/dmr_results.tsv",
        dmr_bed=DSS_BASE + "/dmr_results.bed",
    params:
        script=os.path.join(SCRIPTS_DIR, "run_dss_analysis.R"),
        beds=",".join(
            config["directory"]["output"] + "/pb_cpg_tools/" + s + ".combined.bed.gz"
            for s in ALL_SAMPLES
        ),
        samples=",".join(ALL_SAMPLES),
        groups=",".join(
            "case" if s in config["samples"]["case"] else "control"
            for s in ALL_SAMPLES
        ),
        delta=config.get("dss", {}).get("delta", 0.1),
        pvalue=config.get("dss", {}).get("pvalue", 1e-5),
        fdr=config.get("dss", {}).get("fdr", 0.05),
        min_cpg=config.get("dss", {}).get("min_cpg", 3),
        min_length=config.get("dss", {}).get("min_length", 50),
        merge_dist=config.get("dss", {}).get("merge_dist", 100),
        smoothing=config.get("dss", {}).get("smoothing", True),
        smoothing_span=config.get("dss", {}).get("smoothing_span", 500),
        threads=config.get("dss", {}).get("threads", 4),
    log:
        config["directory"]["output"] + "/logs/dss/dss_analysis.log",
    shell:
        """
        mkdir -p $(dirname {output.dml_tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/envs/R4.5.3/bin/Rscript {params.script} \
            --beds {params.beds} \
            --samples {params.samples} \
            --groups {params.groups} \
            --output-dml-tsv {output.dml_tsv} \
            --output-dml-bed {output.dml_bed} \
            --output-dmr-tsv {output.dmr_tsv} \
            --output-dmr-bed {output.dmr_bed} \
            --delta {params.delta} \
            --pvalue {params.pvalue} \
            --fdr {params.fdr} \
            --min-cpg {params.min_cpg} \
            --min-length {params.min_length} \
            --merge-dist {params.merge_dist} \
            --smoothing {params.smoothing} \
            --smoothing-span {params.smoothing_span} \
            --threads {params.threads} \
            2>&1 | tee {log}
        """


# =============================================================================
# DSS_region_comparison: Region-level methylation comparison using DSS results
#
# Computes per-region:
#   - Mean methylation difference (case - control) across CpGs (pooled across samples)
#   - CpG count per region
#   - DSS DMR/DML overlap counts (total, hyper, hypo)
#
# Supports any number of case/control samples. Processes the region specified
# in config target (e.g., enhancer, cgi).
# =============================================================================
rule DSS_region_comparison:
    input:
        case_beds=expand(
            config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
            sample=config["samples"]["case"]
        ),
        control_beds=expand(
            config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
            sample=config["samples"]["control"]
        ),
        dmr_bed=DSS_BASE + "/dmr_results.bed",
        dmr_tsv=DSS_BASE + "/dmr_results.tsv",
        dml_bed=DSS_BASE + "/dml_results.bed",
        dml_tsv=DSS_BASE + "/dml_results.tsv",
        region_bed=lambda wildcards: config["genome"][wildcards.region],
    output:
        tsv=DSS_BASE + "/{region}_comparison.tsv",
    params:
        script=os.path.join(SCRIPTS_DIR, "dss_region_comparison.py"),
        case_beds=lambda wildcards, input: ",".join(input.case_beds),
        control_beds=lambda wildcards, input: ",".join(input.control_beds),
    log:
        config["directory"]["output"] + "/logs/dss/{region}_comparison.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --case-beds {params.case_beds} \
            --control-beds {params.control_beds} \
            --dmr-bed {input.dmr_bed} \
            --dmr-tsv {input.dmr_tsv} \
            --dml-bed {input.dml_bed} \
            --dml-tsv {input.dml_tsv} \
            --region-bed {input.region_bed} --region-name {wildcards.region} \
            --output-tsv {output.tsv} \
            2>&1 | tee {log}
        """
