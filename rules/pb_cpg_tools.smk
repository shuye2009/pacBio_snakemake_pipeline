# =============================================================================
# Helper functions for conditional input/output based on phasing configuration
# =============================================================================

def get_cpg_input_bam(wildcards):
    """Return haplotagged BAM if phasing enabled, otherwise aligned BAM."""
    if config.get("phasing", {}).get("enabled", False):
        return config["directory"]["output"] + f"/phased/{wildcards.sample}.haplotagged.bam"
    else:
        return config["directory"]["output"] + f"/aligned/{wildcards.sample}.aligned.bam"

def get_cpg_input_bai(wildcards):
    """Return BAM index corresponding to the input BAM."""
    if config.get("phasing", {}).get("enabled", False):
        return config["directory"]["output"] + f"/phased/{wildcards.sample}.haplotagged.bam.bai"
    else:
        return config["directory"]["output"] + f"/aligned/{wildcards.sample}.aligned.bam.bai"

def get_cpg_outputs():
    """
    Return output file dictionary. When phasing is enabled, includes
    haplotype-specific outputs (hap1/hap2) in addition to combined outputs.
    """
    base_outputs = {
        "combined_bed": config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
        "combined_bed_tbi": config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz.tbi",
        "combined_bw": config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bw",
    }
    if config.get("phasing", {}).get("enabled", False):
        base_outputs.update({
            "hap1_bed": config["directory"]["output"] + "/pb_cpg_tools/{sample}.hap1.bed.gz",
            "hap1_bed_tbi": config["directory"]["output"] + "/pb_cpg_tools/{sample}.hap1.bed.gz.tbi",
            "hap1_bw": config["directory"]["output"] + "/pb_cpg_tools/{sample}.hap1.bw",
            "hap2_bed": config["directory"]["output"] + "/pb_cpg_tools/{sample}.hap2.bed.gz",
            "hap2_bed_tbi": config["directory"]["output"] + "/pb_cpg_tools/{sample}.hap2.bed.gz.tbi",
            "hap2_bw": config["directory"]["output"] + "/pb_cpg_tools/{sample}.hap2.bw",
        })
    return base_outputs


# =============================================================================
# Aligned_bam_to_cpg_scores: Extract CpG methylation scores from aligned BAM
# 
# Uses pb-CpG-tools to read methylation probability tags (MM/ML) from aligned
# reads and outputs:
#   - BED files with per-CpG methylation scores
#   - BigWig files for genome browser visualization
# 
# When phasing is enabled, also outputs haplotype-specific methylation files.
# =============================================================================
rule Aligned_bam_to_cpg_scores:
    input:
        bam=get_cpg_input_bam,
        bai=get_cpg_input_bai,
    output:
        **get_cpg_outputs() # The ** is Python's dictionary unpacking operator.
    params:
        output_prefix=config["directory"]["output"] + "/pb_cpg_tools/{sample}",
        threads=config["pb_cpg_tools"]["threads"],
    log:
        config["directory"]["output"] + "/logs/pb_cpg_tools/{sample}.log",
    shell:
        """
        set +u
        source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
        conda activate pb-cpg-tools
        set -u
        mkdir -p $(dirname {params.output_prefix})
        mkdir -p $(dirname {log})
        aligned_bam_to_cpg_scores \
            --bam {input.bam} \
            --output-prefix {params.output_prefix} \
            --threads {params.threads} \
            --pileup-mode model \
            --modsites-mode denovo \
            2>&1 | tee {log}
        conda deactivate
        """
