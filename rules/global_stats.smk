# =============================================================================
# Global methylation statistics
# =============================================================================

import os

ALL_SAMPLES = config["samples"]["case"] + config["samples"]["control"]
GLOBAL_DIR = config["directory"]["output"] + "/global"
DSS_BASE = config["directory"]["output"] + "/dss"


# =============================================================================
# Global_methylation_per_sample: Average methylation per chromosome per sample
#
# Reads the combined BED file from pb-CpG-tools and computes the mean
# mod_score for each chromosome.
# =============================================================================
rule Global_methylation_per_sample:
    input:
        bed=config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
        fai=config["genome"]["fasta"] + ".fai",
    output:
        tsv=GLOBAL_DIR + "/per_sample/{sample}.global_methylation.tsv",
    params:
        script=os.path.join(workflow.basedir, "scripts", "global_methylation_per_sample.py"),
        bin_size=config["global"]["bin_size"],
    log:
        config["directory"]["output"] + "/logs/global/methylation_per_sample_{sample}.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --bed {input.bed} \
            --fai {input.fai} \
            --bin-size {params.bin_size} \
            --output-tsv {output.tsv} \
            2>&1 | tee {log}
        """


# =============================================================================
# Global_density_per_sample: Density plots of mod_score, cov, est_mod_count
#
# Produces one figure per sample with 3 density subplots.
# =============================================================================
rule Global_density_per_sample:
    input:
        bed=config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
    output:
        png=GLOBAL_DIR + "/density/{sample}.density.png",
        pdf=GLOBAL_DIR + "/density/{sample}.density.pdf",
    params:
        script=os.path.join(workflow.basedir, "scripts", "plot_global_density.py"),
    log:
        config["directory"]["output"] + "/logs/global/density_{sample}.log",
    shell:
        """
        mkdir -p $(dirname {output.png})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} per-sample \
            --bed {input.bed} \
            --sample {wildcards.sample} \
            --output-png {output.png} \
            --output-pdf {output.pdf} \
            2>&1 | tee {log}
        """


# =============================================================================
# Global_density_per_chromosome: KDE line plot of mod_score per chromosome
#
# Produces one figure per sample with a subplot for each chromosome to check for quality controls.
# =============================================================================
rule Global_density_per_chromosome:
    input:
        bed=config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
    output:
        png=GLOBAL_DIR + "/density/{sample}.per_chrom_density.png",
        pdf=GLOBAL_DIR + "/density/{sample}.per_chrom_density.pdf",
    params:
        script=os.path.join(workflow.basedir, "scripts", "plot_global_density.py"),
    log:
        config["directory"]["output"] + "/logs/global/density_per_chrom_{sample}.log",
    shell:
        """
        mkdir -p $(dirname {output.png})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} per-chromosome \
            --bed {input.bed} \
            --sample {wildcards.sample} \
            --output-png {output.png} \
            --output-pdf {output.pdf} \
            2>&1 | tee {log}
        """


# =============================================================================
# Global_density_composite: Composite density plot overlaying all samples
#
# Produces one figure with all samples overlaid, colored by case/control.
# =============================================================================
rule Global_density_composite:
    input:
        beds=expand(config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz", sample=ALL_SAMPLES),
    output:
        png=GLOBAL_DIR + "/density/composite_density.png",
        pdf=GLOBAL_DIR + "/density/composite_density.pdf",
    params:
        script=os.path.join(workflow.basedir, "scripts", "plot_global_density.py"),
        samples=" ".join(ALL_SAMPLES),
        case_samples=" ".join(config["samples"]["case"]),
        control_samples=" ".join(config["samples"]["control"]),
    log:
        config["directory"]["output"] + "/logs/global/density_composite.log",
    shell:
        """
        mkdir -p $(dirname {output.png})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} composite \
            --beds {input.beds} \
            --samples {params.samples} \
            --case-samples {params.case_samples} \
            --control-samples {params.control_samples} \
            --output-png {output.png} \
            --output-pdf {output.pdf} \
            2>&1 | tee {log}
        """


# =============================================================================
# Merge_global_methylation: Combine per-sample stats into one matrix
#
# Produces a single TSV with genomic bins as rows and samples as columns.
# Columns: chrom, bin_start, bin_end, sample1, sample2, ...
# =============================================================================
rule Merge_global_methylation:
    input:
        tsvs=expand(GLOBAL_DIR + "/per_sample/{sample}.global_methylation.tsv", sample=ALL_SAMPLES),
        bin_info=GLOBAL_DIR + "/bin_annotation.tsv",
    output:
        tsv=GLOBAL_DIR + "/global_methylation.tsv",
        heatmap_png=GLOBAL_DIR + "/global_methylation_heatmap.png",
        heatmap_pdf=GLOBAL_DIR + "/global_methylation_heatmap.pdf",
        pca_png=GLOBAL_DIR + "/global_methylation_pca.png",
        pca_pdf=GLOBAL_DIR + "/global_methylation_pca.pdf",
    params:
        script=os.path.join(workflow.basedir, "scripts", "merge_global_methylation.py"),
        samples=" ".join(ALL_SAMPLES),
        case_samples=" ".join(config["samples"]["case"]),
        control_samples=" ".join(config["samples"]["control"]),
        bin_size=config["global"]["bin_size"],
    log:
        config["directory"]["output"] + "/logs/global/merge_global_methylation.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --tsvs {input.tsvs} \
            --samples {params.samples} \
            --case-samples {params.case_samples} \
            --control-samples {params.control_samples} \
            --output-tsv {output.tsv} \
            --output-heatmap-png {output.heatmap_png} \
            --output-heatmap-pdf {output.heatmap_pdf} \
            --output-pca-png {output.pca_png} \
            --output-pca-pdf {output.pca_pdf} \
            --bin-size {params.bin_size} \
            --bin-info {input.bin_info} \
            2>&1 | tee {log}
        """


# =============================================================================
# Compute_bin_annotation: Fraction of each bin overlapping annotation tracks
#
# Computes, once for the whole genome, the proportion of each bin covered
# by genes, repeats, centromeres, enhancers, CGIs, and significant DMRs.
# =============================================================================
rule Compute_bin_annotation:
    input:
        fai=config["genome"]["fasta"] + ".fai",
        genes=config["genome"]["gtf"],
        repeats=config["genome"]["repeat"],
        centromeres=config["genome"]["centromere"],
        enhancers=config["genome"]["enhancer"],
        cgi=config["genome"]["cgi"],
        dmrs=VIS_BASE + "/significant_dmrs.bed",
        dml=DSS_BASE + "/dml_results.bed",
        dmr=DSS_BASE + "/dmr_results.bed",
    output:
        tsv=GLOBAL_DIR + "/bin_annotation.tsv",
    params:
        script=os.path.join(workflow.basedir, "scripts", "compute_bin_annotation.py"),
        bin_size=config["global"]["bin_size"],
    log:
        config["directory"]["output"] + "/logs/global/bin_annotation.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --fai {input.fai} \
            --bin-size {params.bin_size} \
            --genes {input.genes} \
            --repeats {input.repeats} \
            --centromeres {input.centromeres} \
            --enhancers {input.enhancers} \
            --cgi {input.cgi} \
            --dmrs {input.dmrs} \
            --dml {input.dml} \
            --dmr {input.dmr} \
            --output-tsv {output.tsv} \
            2>&1 | tee {log}
        """

# =============================================================================
# Plot_methylation_distribution: Distribution of methylation levels per sample
#
# Shows the overall methylation distribution for each sample as violin plots.
# =============================================================================
rule Plot_methylation_distribution:
    input:
        beds=expand(config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz", sample=ALL_SAMPLES),
    output:
        dist=GLOBAL_DIR + "/methylation_distribution.png",
        dist_pdf=GLOBAL_DIR + "/methylation_distribution.pdf",
    params:
        case_samples=",".join(config["samples"]["case"]),
        control_samples=",".join(config["samples"]["control"]),
        bed_dir=config["directory"]["output"] + "/pb_cpg_tools",
        script=os.path.join(SCRIPTS_DIR, "plot_methylation_distribution.py"),
    log:
        config["directory"]["output"] + "/logs/global/distribution.log",
    shell:
        """
        mkdir -p $(dirname {output.dist})
        mkdir -p $(dirname {log})
        
        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --bed-dir {params.bed_dir} \
            --case-samples {params.case_samples} \
            --control-samples {params.control_samples} \
            --output-png {output.dist} \
            --output-pdf {output.dist_pdf} \
            2>&1 | tee {log}
        """
