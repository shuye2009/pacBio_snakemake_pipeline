# =============================================================================
# Pbmm2_index: Create a minimap2 index (.mmi) from the reference genome
# 
# This index is built once and reused for all sample alignments, improving
# performance compared to indexing on-the-fly during each alignment.
# =============================================================================
rule Pbmm2_index:
    input:
        ref=config["genome"]["fasta"],
    output:
        mmi=config["directory"]["output"] + "/reference/genome.mmi",
    params:
        preset=config["pbmm2"]["preset"],
    log:
        config["directory"]["output"] + "/logs/pbmm2/index.log",
    shell:
        """
        set +u
        source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
        conda activate pb-cpg-tools
        set -u
        mkdir -p $(dirname {output.mmi})
        mkdir -p $(dirname {log})
        pbmm2 index \
            {input.ref} \
            {output.mmi} \
            --preset {params.preset} \
            2>&1 | tee {log}
        conda deactivate
        """


# =============================================================================
# Pbmm2_align: Align PacBio HiFi reads to the reference genome
# 
# Takes unaligned BAM (uBAM) files with methylation tags (MM/ML) and produces
# sorted, indexed aligned BAM files. The methylation tags are preserved during
# alignment for downstream CpG analysis.
# =============================================================================
rule Pbmm2_align:
    input:
        ubam=config["directory"]["input"] + "/{sample}.bam",
        mmi=config["directory"]["output"] + "/reference/genome.mmi",
    output:
        bam=config["directory"]["output"] + "/aligned/{sample}.aligned.bam",
        bai=config["directory"]["output"] + "/aligned/{sample}.aligned.bam.bai",
    params:
        threads=config["pbmm2"]["threads"],
    log:
        config["directory"]["output"] + "/logs/pbmm2/{sample}.log",
    shell:
        """
        set +u
        source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
        conda activate pb-cpg-tools
        set -u
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        pbmm2 align \
            {input.mmi} \
            {input.ubam} \
            {output.bam} \
            --sort \
            --num-threads {params.threads} \
            --log-level INFO \
            2>&1 | tee {log}
        conda deactivate
        """
