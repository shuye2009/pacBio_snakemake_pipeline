# =============================================================================
# DeepVariant: Call variants (SNPs/indels) from aligned PacBio HiFi reads
# 
# Uses Google's DeepVariant via Singularity container. Variants are needed
# for haplotype phasing - heterozygous SNPs allow reads to be assigned to
# maternal or paternal chromosomes.
# =============================================================================
rule DeepVariant:
    input:
        bam=config["directory"]["output"] + "/aligned/{sample}.aligned.bam",
        bai=config["directory"]["output"] + "/aligned/{sample}.aligned.bam.bai",
        ref=config["genome"]["fasta"],
    output:
        vcf=config["directory"]["output"] + "/variants/{sample}.deepvariant.vcf.gz",
        gvcf=config["directory"]["output"] + "/variants/{sample}.deepvariant.g.vcf.gz",
    params:
        threads=config["deepvariant"]["threads"],
        model_type=config["deepvariant"]["model_type"],
        tmp_dir=config["directory"]["output"] + "/variants/{sample}_tmp",
        sif=config["deepvariant"]["singularity_image"],
    log:
        config["directory"]["output"] + "/logs/deepvariant/{sample}.log",
    shell:
        """
        mkdir -p {params.tmp_dir}
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        module load singularity
        # Check if GPU is available on this node
        if nvidia-smi &>/dev/null; then
            echo "GPU detected, running with --nv flag"
            GPU_FLAG="--nv"
        else
            echo "No GPU detected, running on CPU only"
            GPU_FLAG=""
        fi
        set +u
        singularity exec \
            $GPU_FLAG \
            --bind /cluster \
            --env TMPDIR={params.tmp_dir} \
            --env XLA_FLAGS=--xla_gpu_cuda_data_dir=/usr/local/cuda \
            --env LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH \
            {params.sif} \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model_type} \
            --ref={input.ref} \
            --reads={input.bam} \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --num_shards={params.threads} \
            --intermediate_results_dir={params.tmp_dir} \
            2>&1 | tee {log}
        rm -rf {params.tmp_dir}
        set -u
        """


# =============================================================================
# HiPhase: Phase variants and haplotag reads
# 
# Uses the variants from DeepVariant to determine haplotype blocks and assigns
# each read to haplotype 1 or 2 (HP tag). The haplotagged BAM enables
# haplotype-resolved methylation analysis in downstream steps.
# =============================================================================
rule HiPhase:
    input:
        bam=config["directory"]["output"] + "/aligned/{sample}.aligned.bam",
        bai=config["directory"]["output"] + "/aligned/{sample}.aligned.bam.bai",
        vcf=config["directory"]["output"] + "/variants/{sample}.deepvariant.vcf.gz",
        ref=config["genome"]["fasta"],
    output:
        phased_vcf=config["directory"]["output"] + "/phased/{sample}.phased.vcf.gz",
        haplotagged_bam=config["directory"]["output"] + "/phased/{sample}.haplotagged.bam",
        summary=config["directory"]["output"] + "/phased/{sample}.hiphase_stats.tsv",
        blocks=config["directory"]["output"] + "/phased/{sample}.hiphase_blocks.tsv",
    params:
        threads=config["hiphase"]["threads"],
    log:
        config["directory"]["output"] + "/logs/hiphase/{sample}.log",
    shell:
        """
        set +u
        source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
        conda activate hiphase
        set -u
        mkdir -p $(dirname {output.phased_vcf})
        mkdir -p $(dirname {log})
        hiphase \
            --threads {params.threads} \
            --bam {input.bam} \
            --vcf {input.vcf} \
            --reference {input.ref} \
            --output-vcf {output.phased_vcf} \
            --output-bam {output.haplotagged_bam} \
            --stats-file {output.summary} \
            --blocks-file {output.blocks} \
            2>&1 | tee {log}
        conda deactivate
        """


# =============================================================================
# Index_haplotagged_bam: Create BAM index for the haplotagged BAM file
# 
# Required for downstream tools to efficiently access the haplotagged BAM.
# =============================================================================
rule Index_haplotagged_bam:
    input:
        bam=config["directory"]["output"] + "/phased/{sample}.haplotagged.bam",
    output:
        bai=config["directory"]["output"] + "/phased/{sample}.haplotagged.bam.bai",
    shell:
        """
        module load samtools/1.20
        samtools index {input.bam}
        module unload samtools/1.20
        """
