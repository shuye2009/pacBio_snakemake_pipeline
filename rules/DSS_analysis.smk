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
DSS_DIR     = config["directory"]["output"] + "/dss"


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
        dml_tsv=DSS_DIR + "/dml_results.tsv",
        dml_bed=DSS_DIR + "/dml_results.bed",
        dmr_tsv=DSS_DIR + "/dmr_results.tsv",
        dmr_bed=DSS_DIR + "/dmr_results.bed",
    params:
        script=os.path.join(SCRIPTS_DIR, "run_dss_analysis.R"),
        beds=",".join(
            config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz"
            for sample in ALL_SAMPLES
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
