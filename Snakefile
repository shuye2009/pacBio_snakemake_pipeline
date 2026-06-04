import os
os.environ["TMPDIR"] = "/cluster/projects/hardinggroup/Shuye/tmp"
os.environ["TEMPDIR"] = "/cluster/projects/hardinggroup/Shuye/tmp"

configfile: "config/config.yaml"

ALL_SAMPLES = config["samples"]["case"] + config["samples"]["control"]
PHASING_ENABLED = config.get("phasing", {}).get("enabled", False)

# Dynamically set methbat regions from target (e.g., enhancer, cgi, centromere, repeat)
TARGET = config.get("target", "enhancer")
if TARGET in config.get("genome", {}):
    config["methbat"]["regions"] = config["genome"][TARGET]

# Output subdirectories: base (target-independent) vs target-scoped
METHBAT_BASE = config["directory"]["output"] + "/methbat"
VIS_BASE = config["directory"]["output"] + "/visualization"
FUNC_BASE = config["directory"]["output"] + "/functional_analysis"
METHBAT_DIR = METHBAT_BASE + "/" + TARGET
VIS_DIR = VIS_BASE + "/" + TARGET
FUNC_DIR = FUNC_BASE + "/" + TARGET
DSS_BASE = config["directory"]["output"] + "/dss"

include: "rules/pbmm2.smk"
if PHASING_ENABLED:
    include: "rules/phasing.smk"
include: "rules/pb_cpg_tools.smk"
include: "rules/methbat.smk"
include: "rules/visualization.smk"
include: "rules/bio_insights.smk"
include: "rules/global_stats.smk"
include: "rules/DSS_analysis.smk"

def get_all_outputs():
    """Return all final outputs, including haplotype-specific when phasing is enabled."""
    outputs = [
        config["directory"]["output"] + "/global/global_methylation.tsv",
        config["directory"]["output"] + "/global/bin_annotation.tsv",
        config["directory"]["output"] + "/global/density/composite_density.png",
    ] + expand(config["directory"]["output"] + "/global/density/{sample}.per_chrom_density.png", sample=ALL_SAMPLES) + [
        METHBAT_BASE + "/signature.signature_regions.bed",
        METHBAT_BASE + "/signature.signature_stats.tsv",
        METHBAT_DIR + "/region_cohort_comparison.tsv",
        VIS_DIR + "/methylation_report.html",
        VIS_DIR + "/igv_significant_regions.html",
        VIS_BASE + "/igv_significant_dmrs.html",
        VIS_BASE + "/igv_dss_dmrs.html",
        VIS_BASE + "/igv_dss_dmls.html",
        FUNC_BASE + "/dmr_annotations.tsv",
        FUNC_BASE + "/dmr_annotation_distribution.png",
        FUNC_BASE + "/dmr_annotation_distribution.pdf",
        FUNC_BASE + "/dmr_gene_associations.tsv",
        FUNC_BASE + "/hyper_dmr_go_enrichment.tsv",
        FUNC_BASE + "/hypo_dmr_go_enrichment.tsv",
        FUNC_BASE + "/hyper_dmr_go_enrichment.png",
        FUNC_BASE + "/hypo_dmr_go_enrichment.png",
        FUNC_BASE + "/hyper_dmr_kegg_enrichment.tsv",
        FUNC_BASE + "/hypo_dmr_kegg_enrichment.tsv",
        FUNC_BASE + "/hyper_dmr_kegg_enrichment.png",
        FUNC_BASE + "/hypo_dmr_kegg_enrichment.png",
        FUNC_BASE + "/dmr_motif_enrichment.tsv",
        FUNC_BASE + "/dmr_motif_enrichment.png",
        FUNC_BASE + "/dmr_motif_enrichment.pdf",
        FUNC_BASE + "/dmr_motif_enrichment.html",
        DSS_BASE + "/dml_results.tsv",
        DSS_BASE + "/dml_results.bed",
        DSS_BASE + "/dmr_results.tsv",
        DSS_BASE + "/dmr_results.bed",
    ]
    if PHASING_ENABLED:
        outputs.extend([
            METHBAT_BASE + "/signature_hap1.signature_regions.bed",
            METHBAT_BASE + "/signature_hap1.signature_stats.tsv",
            METHBAT_BASE + "/signature_hap2.signature_regions.bed",
            METHBAT_BASE + "/signature_hap2.signature_stats.tsv",
            METHBAT_DIR + "/region_cohort_comparison_hap1.tsv",
            METHBAT_DIR + "/region_cohort_comparison_hap2.tsv",
            # ASM analysis outputs (separate for case and control groups)
            METHBAT_BASE + "/asm_case.meth_regions.bed",
            METHBAT_BASE + "/asm_control.meth_regions.bed",
            METHBAT_BASE + "/asm_haplotype_comparison.tsv",
            METHBAT_BASE + "/asm_group_comparison.tsv",
            # ASM IGV reports
            VIS_BASE + "/igv_significant_asm.html",
            VIS_BASE + "/igv_differential_asm.html",
        ])
    return outputs

rule all:
    input:
        get_all_outputs()
