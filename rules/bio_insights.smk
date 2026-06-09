# =============================================================================
# Biological insights rules for DMR annotation and enrichment analysis
# =============================================================================

import os

SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")


if ENOUGH_SAMPLES:

    # =========================================================================
    # Annotate_dmrs: Annotate significant DMRs with custom T2T-CHM13 annotations
    #
    # Uses the annotatr R package to annotate DMRs against custom annotations
    # built from the GTF file (genes, exons, introns, CDS, UTRs, promoters) and
    # CGI regions from the config.
    # =========================================================================
    rule Annotate_dmrs:
        input:
            dmr_bed=VIS_BASE + "/significant_dmrs.bed",
            dmr_tsv=VIS_BASE + "/significant_dmrs.tsv",
            genome_fasta=config["genome"]["fasta"],
            gff3=config["genome"]["gff3"],
            cgi_bed=config["genome"]["cgi"],
            enhancer_bed=config["genome"]["enhancer"],
            eglink_bed=config["genome"]["eglink"],
        output:
            tsv=FUNC_BASE + "/dmr_annotations.tsv",
            png=FUNC_BASE + "/dmr_annotation_distribution.png",
            pdf=FUNC_BASE + "/dmr_annotation_distribution.pdf",
            gene_assoc=FUNC_BASE + "/dmr_gene_associations.tsv",
            hyper_go=FUNC_BASE + "/hyper_dmr_go_enrichment.tsv",
            hypo_go=FUNC_BASE + "/hypo_dmr_go_enrichment.tsv",
            hyper_go_png=FUNC_BASE + "/hyper_dmr_go_enrichment.png",
            hypo_go_png=FUNC_BASE + "/hypo_dmr_go_enrichment.png",
            hyper_kegg=FUNC_BASE + "/hyper_dmr_kegg_enrichment.tsv",
            hypo_kegg=FUNC_BASE + "/hypo_dmr_kegg_enrichment.tsv",
            hyper_kegg_png=FUNC_BASE + "/hyper_dmr_kegg_enrichment.png",
            hypo_kegg_png=FUNC_BASE + "/hypo_dmr_kegg_enrichment.png",
        params:
            script=os.path.join(SCRIPTS_DIR, "annotate_dmrs.R"),
        log:
            config["directory"]["output"] + "/logs/functional_analysis/annotate_dmrs.log",
        shell:
            """
            mkdir -p $(dirname {output.tsv})
            mkdir -p $(dirname {log})

            /cluster/home/t128737uhn/miniconda3/envs/R4.5.3/bin/Rscript {params.script} \
                --dmr-bed {input.dmr_bed} \
                --dmr-tsv {input.dmr_tsv} \
                --genome-fasta {input.genome_fasta} \
                --gff3 {input.gff3} \
                --cgi-bed {input.cgi_bed} \
                --enhancer-bed {input.enhancer_bed} \
                --eglink-bed {input.eglink_bed} \
                --output-tsv {output.tsv} \
                --output-png {output.png} \
                --output-pdf {output.pdf} \
                --gene-assoc-tsv {output.gene_assoc} \
                --hyper-go-tsv {output.hyper_go} \
                --hypo-go-tsv {output.hypo_go} \
                --hyper-go-png {output.hyper_go_png} \
                --hypo-go-png {output.hypo_go_png} \
                --hyper-kegg-tsv {output.hyper_kegg} \
                --hypo-kegg-tsv {output.hypo_kegg} \
                --hyper-kegg-png {output.hyper_kegg_png} \
                --hypo-kegg-png {output.hypo_kegg_png} \
                2>&1 | tee {log}
            """


if ENOUGH_SAMPLES:

    # =========================================================================
    # Motif_enrichment: TFBS motif enrichment analysis in DMRs
    #
    # Uses JASPAR2024 and TFBSTools to scan for transcription factor binding
    # site motifs in DMR sequences and tests for enrichment against
    # dinucleotide-shuffled background sequences.
    # =========================================================================
    rule Motif_enrichment:
        input:
            dmr_bed=VIS_BASE + "/significant_dmrs.bed",
            fasta=config["genome"]["fasta"],
        output:
            tsv=FUNC_BASE + "/dmr_motif_enrichment.tsv",
            png=FUNC_BASE + "/dmr_motif_enrichment.png",
            pdf=FUNC_BASE + "/dmr_motif_enrichment.pdf",
            html=FUNC_BASE + "/dmr_motif_enrichment.html",
        params:
            script=os.path.join(SCRIPTS_DIR, "motif_enrichment.R"),
            n_background=10,
            pvalue_cutoff=0.05,
            threads=16,
            jaspar_db=config.get("genome", {}).get("jaspar_db", ""),
        log:
            config["directory"]["output"] + "/logs/functional_analysis/motif_enrichment.log",
        shell:
            """
            mkdir -p $(dirname {output.tsv})
            mkdir -p $(dirname {log})

            /cluster/home/t128737uhn/miniconda3/envs/R4.5.3/bin/Rscript {params.script} \
                --dmr-bed {input.dmr_bed} \
                --fasta {input.fasta} \
                --jaspar-db {params.jaspar_db} \
                --output-tsv {output.tsv} \
                --output-png {output.png} \
                --output-pdf {output.pdf} \
                --output-html {output.html} \
                --n-background {params.n_background} \
                --pvalue-cutoff {params.pvalue_cutoff} \
                --threads {params.threads} \
                2>&1 | tee {log}
            """


# =============================================================================
# Convert_dss_dmr_tsv: Convert DSS DMR TSV columns for annotation pipeline
#
# DSS callDMR outputs lowercase column names (chr, start, end, diff.Methy)
# but annotate_dmrs.R expects uppercase (CHROM, START, END, DELTA).
# =============================================================================
rule Convert_dss_dmr_tsv:
    input:
        tsv=DSS_BASE + "/dmr_results.tsv",
    output:
        tsv=FUNC_BASE + "/dss_dmr_results.converted.tsv",
    params:
        script=os.path.join(SCRIPTS_DIR, "convert_dss_dmr_tsv.py"),
    log:
        config["directory"]["output"] + "/logs/functional_analysis/convert_dss_dmr_tsv.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/bin/python {params.script} \
            --input-tsv {input.tsv} \
            --output-tsv {output.tsv} \
            2>&1 | tee {log}
        """


# =============================================================================
# Annotate_dss_dmrs: Annotate DSS DMRs with custom T2T-CHM13 annotations
#
# Uses the same annotate_dmrs.R pipeline but with DSS-called DMRs instead
# of methbat significant DMRs. Includes GO/KEGG enrichment for associated genes.
# =============================================================================
rule Annotate_dss_dmrs:
    input:
        dmr_bed=DSS_BASE + "/dmr_results.bed",
        dmr_tsv=FUNC_BASE + "/dss_dmr_results.converted.tsv",
        genome_fasta=config["genome"]["fasta"],
        gff3=config["genome"]["gff3"],
        cgi_bed=config["genome"]["cgi"],
        enhancer_bed=config["genome"]["enhancer"],
        eglink_bed=config["genome"]["eglink"],
    output:
        tsv=FUNC_BASE + "/dss_dmr_annotations.tsv",
        png=FUNC_BASE + "/dss_dmr_annotation_distribution.png",
        pdf=FUNC_BASE + "/dss_dmr_annotation_distribution.pdf",
        gene_assoc=FUNC_BASE + "/dss_dmr_gene_associations.tsv",
        hyper_go=FUNC_BASE + "/dss_hyper_dmr_go_enrichment.tsv",
        hypo_go=FUNC_BASE + "/dss_hypo_dmr_go_enrichment.tsv",
        hyper_go_png=FUNC_BASE + "/dss_hyper_dmr_go_enrichment.png",
        hypo_go_png=FUNC_BASE + "/dss_hypo_dmr_go_enrichment.png",
        hyper_kegg=FUNC_BASE + "/dss_hyper_dmr_kegg_enrichment.tsv",
        hypo_kegg=FUNC_BASE + "/dss_hypo_dmr_kegg_enrichment.tsv",
        hyper_kegg_png=FUNC_BASE + "/dss_hyper_dmr_kegg_enrichment.png",
        hypo_kegg_png=FUNC_BASE + "/dss_hypo_dmr_kegg_enrichment.png",
    params:
        script=os.path.join(SCRIPTS_DIR, "annotate_dmrs.R"),
    log:
        config["directory"]["output"] + "/logs/functional_analysis/annotate_dss_dmrs.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/envs/R4.5.3/bin/Rscript {params.script} \
            --dmr-bed {input.dmr_bed} \
            --dmr-tsv {input.dmr_tsv} \
            --genome-fasta {input.genome_fasta} \
            --gff3 {input.gff3} \
            --cgi-bed {input.cgi_bed} \
            --enhancer-bed {input.enhancer_bed} \
            --eglink-bed {input.eglink_bed} \
            --output-tsv {output.tsv} \
            --output-png {output.png} \
            --output-pdf {output.pdf} \
            --gene-assoc-tsv {output.gene_assoc} \
            --hyper-go-tsv {output.hyper_go} \
            --hypo-go-tsv {output.hypo_go} \
            --hyper-go-png {output.hyper_go_png} \
            --hypo-go-png {output.hypo_go_png} \
            --hyper-kegg-tsv {output.hyper_kegg} \
            --hypo-kegg-tsv {output.hypo_kegg} \
            --hyper-kegg-png {output.hyper_kegg_png} \
            --hypo-kegg-png {output.hypo_kegg_png} \
            2>&1 | tee {log}
        """


# =============================================================================
# Motif_enrichment_dss_dmrs: TFBS motif enrichment in DSS-called DMRs
#
# Uses JASPAR2024 and TFBSTools to scan for transcription factor binding
# site motifs in DSS DMR sequences and tests for enrichment against
# dinucleotide-shuffled background sequences.
# =============================================================================
rule Motif_enrichment_dss_dmrs:
    input:
        dmr_bed=DSS_BASE + "/dmr_results.bed",
        fasta=config["genome"]["fasta"],
    output:
        tsv=FUNC_BASE + "/dss_dmr_motif_enrichment.tsv",
        png=FUNC_BASE + "/dss_dmr_motif_enrichment.png",
        pdf=FUNC_BASE + "/dss_dmr_motif_enrichment.pdf",
        html=FUNC_BASE + "/dss_dmr_motif_enrichment.html",
    params:
        script=os.path.join(SCRIPTS_DIR, "motif_enrichment.R"),
        n_background=10,
        pvalue_cutoff=0.05,
        threads=16,
        jaspar_db=config.get("genome", {}).get("jaspar_db", ""),
    log:
        config["directory"]["output"] + "/logs/functional_analysis/motif_enrichment_dss_dmrs.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {log})

        /cluster/home/t128737uhn/miniconda3/envs/R4.5.3/bin/Rscript {params.script} \
            --dmr-bed {input.dmr_bed} \
            --fasta {input.fasta} \
            --jaspar-db {params.jaspar_db} \
            --output-tsv {output.tsv} \
            --output-png {output.png} \
            --output-pdf {output.pdf} \
            --output-html {output.html} \
            --n-background {params.n_background} \
            --pvalue-cutoff {params.pvalue_cutoff} \
            --threads {params.threads} \
            2>&1 | tee {log}
        """