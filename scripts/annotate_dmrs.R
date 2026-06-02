#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(annotatr)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Rsamtools)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(enrichplot)
})

args <- commandArgs(trailingOnly = TRUE)
params <- list()
for (i in seq(1, length(args), by = 2)) {
  params[[sub("^--", "", args[i])]] <- args[i + 1]
}

dmr_bed       <- params[["dmr-bed"]]
dmr_tsv       <- params[["dmr-tsv"]]
genome_fasta  <- params[["genome-fasta"]]
gff3_file     <- params[["gff3"]]
cgi_bed       <- params[["cgi-bed"]]
enhancer_bed  <- params[["enhancer-bed"]]
eglink_bed    <- params[["eglink-bed"]]
output_tsv    <- params[["output-tsv"]]
output_png    <- params[["output-png"]]
output_pdf    <- params[["output-pdf"]]
gene_assoc_tsv<- params[["gene-assoc-tsv"]]
hyper_go_tsv  <- params[["hyper-go-tsv"]]
hypo_go_tsv   <- params[["hypo-go-tsv"]]
hyper_go_png  <- params[["hyper-go-png"]]
hypo_go_png   <- params[["hypo-go-png"]]
hyper_kegg_tsv<- params[["hyper-kegg-tsv"]]
hypo_kegg_tsv <- params[["hypo-kegg-tsv"]]
hyper_kegg_png<- params[["hyper-kegg-png"]]
hypo_kegg_png <- params[["hypo-kegg-png"]]

# for debugging
if(0){
  dmr_bed <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/visualization/significant_dmrs.bed"
  dmr_tsv <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/visualization/significant_dmrs.tsv"
  gff3_file <- "/cluster/projects/hardinggroup/Shuye/reference_genome/human/T2T_CHM13/Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3"
  cgi_bed <- "/cluster/projects/hardinggroup/Shuye/resource/chm13v2.0_CGI.bed"
  enhancer_bed <- "/cluster/projects/hardinggroup/Shuye/resource/ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T_6col.bed"
  eglink_bed <- "/cluster/projects/hardinggroup/Shuye/resource/ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T.bed"
  output_tsv <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/significant_dmrs_annotated.tsv"
  output_pdf <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/significant_dmrs_annotated.pdf"
  output_png <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/significant_dmrs_annotated.png"
  gene_assoc_tsv <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/dmr_gene_associations.tsv"
  hyper_go_tsv <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/hyper_dmr_go_enrichment.tsv"
  hypo_go_tsv <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/hypo_dmr_go_enrichment.tsv"
  hyper_go_png <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/hyper_dmr_go_enrichment.png"
  hypo_go_png <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/hypo_dmr_go_enrichment.png"
  hyper_kegg_tsv <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/hyper_dmr_kegg_enrichment.tsv"
  hypo_kegg_tsv <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/hypo_dmr_kegg_enrichment.tsv"
  hyper_kegg_png <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/hyper_dmr_kegg_enrichment.png"
  hypo_kegg_png <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/hypo_dmr_kegg_enrichment.png"
}


# =========================================================================
# 1. Load DMRs with delta values from TSV
# =========================================================================
cat("Loading DMRs:", dmr_bed, "\n")
dmr_raw <- import(dmr_bed, format = "bed")

cat("Loading DMR delta values from:", dmr_tsv, "\n")
dmr_tsv_df <- read.table(dmr_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dmr_tsv_gr <- makeGRangesFromDataFrame(dmr_tsv_df,
  seqnames.field = "CHROM", start.field = "START", end.field = "END", keep.extra.columns = TRUE)
hits <- findOverlaps(dmr_raw, dmr_tsv_gr)
dmr_raw$delta <- NA_real_
dmr_raw$delta[queryHits(hits)] <- dmr_tsv_gr$DELTA[subjectHits(hits)]
dmr_raw$hyper <- dmr_raw$delta > 0
cat(sprintf("Loaded %d DMRs (%d hyper, %d hypo)\n",
            length(dmr_raw), sum(dmr_raw$hyper, na.rm = TRUE), sum(!dmr_raw$hyper, na.rm = TRUE)))

dmrs <- dmr_raw
dmrs$id <- paste0(seqnames(dmrs), ":", start(dmrs), "-", end(dmrs))
dmrs$type <- "dmr"

# =========================================================================
# 2. Build annotations from GFF3
# =========================================================================
cat("Building annotations from GFF3:", gff3_file, "\n")
gff_gr <- import(gff3_file, format = "gff3")

# Build transcript -> gene mapping from GFF3 parent chain
# exon/CDS/UTR Parent=transcript:ENST... -> mRNA ID=transcript:ENST... -> mRNA Parent=gene:ENSG...
gff_type <- as.character(gff_gr$type)
cat(sprintf("  GFF3 feature types: %s\n", paste(unique(gff_type), collapse = ", ")))
mrna_gr <- gff_gr[gff_type == "mRNA"]
cat(sprintf("  mRNA entries: %d\n", length(mrna_gr)))
tx_to_gene <- setNames(as.character(mrna_gr$Parent), as.character(mrna_gr$ID))
cat(sprintf("  transcript->gene map: %d entries\n", length(tx_to_gene)))

resolve_gene <- function(parent_ids) {
  parent_ids <- as.character(parent_ids)
  gene_ids <- tx_to_gene[parent_ids]
  gene_ids[is.na(gene_ids)] <- parent_ids[is.na(gene_ids)]
  gene_ids
}

make_gr <- function(gr, type_name) {
  if (length(gr) == 0) return(NULL)
  gr$id <- resolve_gene(gr$Parent)
  gr$type <- type_name
  gr
}

# Genes: use gene_id attribute directly, keep Name for gene symbol lookup
genes_gr <- gff_gr[gff_type == "gene"]
genes_gr$id <- as.character(genes_gr$gene_id)
genes_gr$type <- "T2T_genes"
gene_name_map <- setNames(as.character(genes_gr$Name), as.character(genes_gr$gene_id))

harmonize_seqlevels <- function(gr, ref_gr) {
  if (is.null(gr) || length(gr) == 0 || is.null(ref_gr) || length(ref_gr) == 0) {
    return(gr)
  }
  ref_style <- tryCatch(GenomeInfoDb::seqlevelsStyle(ref_gr)[1], error = function(e) NA_character_)
  if (!is.na(ref_style) && nzchar(ref_style)) {
    gr <- tryCatch({
      GenomeInfoDb::seqlevelsStyle(gr) <- ref_style
      gr
    }, error = function(e) gr)
  }
  common_seqlevels <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(ref_gr))
  if (length(common_seqlevels) > 0) {
    gr <- GenomeInfoDb::keepSeqlevels(gr, common_seqlevels, pruning.mode = "coarse")
  }
  gr
}

genes_gr <- harmonize_seqlevels(genes_gr, dmr_raw)

exons_gr    <- make_gr(gff_gr[gff_type == "exon"], "T2T_exons")
cds_gr      <- make_gr(gff_gr[gff_type == "CDS"], "T2T_cds")
five_gr     <- make_gr(gff_gr[gff_type == "five_prime_UTR"], "T2T_5UTRs")
three_gr    <- make_gr(gff_gr[gff_type == "three_prime_UTR"], "T2T_3UTRs")

introns_gr <- NULL
if (!is.null(genes_gr) && !is.null(exons_gr)) {
  introns_gr <- GenomicRanges::setdiff(genes_gr, exons_gr, ignore.strand = FALSE)
  if (length(introns_gr) > 0) {
    hits <- findOverlaps(introns_gr, genes_gr, ignore.strand = FALSE, select = "first")
    introns_gr$id <- genes_gr$id[hits]
    introns_gr$type <- "T2T_introns"
  }
}

promoters_gr <- NULL
if (!is.null(genes_gr)) {
  promoters_gr <- promoters(genes_gr, upstream = 2000, downstream = 0)
  promoters_gr$id <- promoters_gr$id
  promoters_gr$type <- "T2T_promoters"
  promoters_gr <- harmonize_seqlevels(promoters_gr, dmr_raw)
}

annotation_list <- list()
add_annot <- function(gr, name) {
  if (!is.null(gr) && length(gr) > 0) {
    annotation_list[[name]] <<- GRanges(
      seqnames = seqnames(gr),
      ranges   = ranges(gr),
      strand   = strand(gr),
      id       = as.character(gr$id),
      type     = as.character(gr$type)
    )
    cat(sprintf("  %s: %d\n", name, length(gr)))
  }
}

add_annot(genes_gr, "T2T_genes")
add_annot(exons_gr, "T2T_exons")
add_annot(introns_gr, "T2T_introns")
add_annot(cds_gr, "T2T_cds")
add_annot(five_gr, "T2T_5UTRs")
add_annot(three_gr, "T2T_3UTRs")
add_annot(promoters_gr, "T2T_promoters")

cat("Loading enhancers:", enhancer_bed, "\n")
enhancer_gr <- import(enhancer_bed, format = "bed")
if (length(enhancer_gr) > 0) {
  enhancer_gr$id <- paste0(seqnames(enhancer_gr), ":", start(enhancer_gr), "-", end(enhancer_gr))
  enhancer_gr$type <- "T2T_enhancer"
  annotation_list[["T2T_enhancer"]] <- enhancer_gr
  cat(sprintf("  T2T_enhancer: %d\n", length(enhancer_gr)))
}

cat("Loading CGI:", cgi_bed, "\n")
cgi_df <- read.table(cgi_bed, sep = "\t", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
cgi_df <- cgi_df[, 1:3]
colnames(cgi_df) <- c("chrom", "start", "end")
cgi_df <- cgi_df[!is.na(cgi_df$start) & !is.na(cgi_df$end), ]
cgi_gr <- makeGRangesFromDataFrame(cgi_df, starts.in.df.are.0based = TRUE)
if (length(cgi_gr) > 0) {
  cgi_gr$id <- paste0(seqnames(cgi_gr), ":", start(cgi_gr), "-", end(cgi_gr))
  cgi_gr$type <- "T2T_cgi"
  annotation_list[["T2T_cgi"]] <- cgi_gr
  cat(sprintf("  T2T_cgi: %d\n", length(cgi_gr)))
}

all_annotations <- unlist(GRangesList(annotation_list))
cat(sprintf("Total annotations: %d\n", length(all_annotations)))

# =========================================================================
# 3. Annotate DMRs
# =========================================================================
cat("Annotating DMRs...\n")
dmrs_annotated <- annotate_regions(dmrs, all_annotations, ignore.strand = TRUE, quiet = FALSE)

cat("Summarizing...\n")

# Calculate total genome size from genome FASTA
# Use total genome length as background for enrichment testing
genome_total_width <- NULL
if (!is.null(genome_fasta) && file.exists(genome_fasta)) {
  tryCatch({
    fa_idx <- indexFa(genome_fasta)
    genome_total_width <- sum(as.numeric(seqlengths(fa_idx)))
    cat(sprintf("Genome total length from FASTA: %.0f bp\n", genome_total_width))
  }, error = function(e) {
    cat("Warning: could not read genome FASTA index, using annotation union as fallback\n")
  })
}
if (is.null(genome_total_width) || genome_total_width == 0) {
  # Fallback: use union of all annotation ranges
  genome_total_width <- sum(as.numeric(width(all_annotations)))
  cat(sprintf("Using annotation union as background: %.0f bp\n", genome_total_width))
}
ann_types <- unique(all_annotations$type)

ann_summary <- do.call(rbind, lapply(ann_types, function(t) {
  sub_ann <- all_annotations[all_annotations$type == t]
  hits <- findOverlaps(dmrs, sub_ann, ignore.strand = TRUE)
  n_dmr_obs <- length(unique(queryHits(hits)))

  # Expected DMRs based on proportion of annotated genome
  ann_width <- sum(as.numeric(width(sub_ann)))
  ann_prop <- ann_width / genome_total_width
  n_dmr_exp <- length(dmrs) * ann_prop

  # Binomial test for enrichment/depletion
  # H0: probability of DMR in this annotation = ann_prop
  # Using genome as background, test if observed > expected
  p_val <- NA_real_
  if (n_dmr_obs > 0 && ann_prop > 0) {
    # Two-sided binomial test
    bt <- binom.test(n_dmr_obs, length(dmrs), p = ann_prop, alternative = "two.sided")
    p_val <- bt$p.value
  } else if (n_dmr_obs == 0 && ann_prop > 0) {
    # If no DMRs observed but annotation exists, test for depletion
    bt <- binom.test(0, length(dmrs), p = ann_prop, alternative = "less")
    p_val <- bt$p.value
  }

  data.frame(
    annotation_type = t,
    n_annotations   = length(sub_ann),
    ann_width_bp    = ann_width,
    ann_prop        = round(ann_prop * 100, 4),
    n_dmrs_obs      = n_dmr_obs,
    n_dmrs_exp      = round(n_dmr_exp, 2),
    pct_dmrs        = round(100 * n_dmr_obs / length(dmrs), 1),
    fold_enrich     = if (n_dmr_exp > 0) round((n_dmr_obs+1) / (n_dmr_exp+1), 2) else NA_real_,
    p_value         = p_val,
    stringsAsFactors = FALSE
  )
}))

# Multiple testing correction
ann_summary$padj <- p.adjust(ann_summary$p_value, method = "BH")

# Significance and direction
ann_summary$significance <- ifelse(
  ann_summary$padj < 0.001, "***",
  ifelse(ann_summary$padj < 0.01, "**",
  ifelse(ann_summary$padj < 0.05, "*", "ns"))
)
ann_summary$direction <- ifelse(
  ann_summary$n_dmrs_obs > ann_summary$n_dmrs_exp, "enriched",
  ifelse(ann_summary$n_dmrs_obs < ann_summary$n_dmrs_exp, "depleted", "no_change")
)

ann_summary <- ann_summary[order(-ann_summary$n_dmrs_obs), ]
rownames(ann_summary) <- NULL

cat("DMR summary:\n")
print(ann_summary)

cat("Writing:", output_tsv, "\n")
dmr_gr_out <- granges(dmrs_annotated)
dmr_mcols <- mcols(dmrs_annotated)
if ("annot" %in% names(dmr_mcols)) {
  annot_col <- dmr_mcols$annot
  if (inherits(annot_col, "GRanges")) {
    annot_types <- mcols(annot_col)$type
    dmr_mcols$annot <- if (!is.null(annot_types)) as.character(annot_types) else NA_character_
  } else {
    dmr_mcols$annot <- as.character(annot_col)
  }
}
df_out <- data.frame(
  seqnames = as.character(seqnames(dmr_gr_out)),
  start    = start(dmr_gr_out),
  end      = end(dmr_gr_out),
  width    = width(dmr_gr_out),
  strand   = as.character(strand(dmr_gr_out)),
  as.data.frame(dmr_mcols),
  row.names = NULL,
  check.names = FALSE
)
write.table(df_out, output_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# =========================================================================
# 4. DMR-Gene association via promoters and enhancers
# =========================================================================
cat("Associating DMRs with genes...\n")

promoter_genes <- data.frame(
  dmr_idx   = integer(),
  gene_id   = character(),
  gene_name = character(),
  source    = character(),
  stringsAsFactors = FALSE
)
if (!is.null(promoters_gr) && length(promoters_gr) > 0) {
  promoter_hits <- findOverlaps(dmr_raw, promoters_gr, ignore.strand = TRUE)
  if (length(promoter_hits) > 0) {
    promoter_gene_ids <- promoters_gr$id[subjectHits(promoter_hits)]
    promoter_genes <- data.frame(
      dmr_idx   = queryHits(promoter_hits),
      gene_id   = promoter_gene_ids,
      gene_name = unname(gene_name_map[promoter_gene_ids]),
      source    = "promoter",
      stringsAsFactors = FALSE
    )
  }
}

eglink_genes <- data.frame()
if (!is.null(eglink_bed) && file.exists(eglink_bed)) {
  cat("Loading enhancer-gene links:", eglink_bed, "\n")
  eglink_df <- read.table(
    eglink_bed,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if ("#chr" %in% names(eglink_df)) {
    names(eglink_df)[names(eglink_df) == "#chr"] <- "chrom"
  }
  required_cols <- c("chrom", "start", "end", "TargetGene", "TargetGeneEnsemblID")
  if (all(required_cols %in% names(eglink_df))) {
    eglink_df <- eglink_df[!is.na(eglink_df$chrom) & !is.na(eglink_df$start) & !is.na(eglink_df$end), required_cols, drop = FALSE]
    if (nrow(eglink_df) > 0) {
      eglink_df$start <- as.numeric(eglink_df$start)
      eglink_df$end <- as.numeric(eglink_df$end)
      eglink_df <- eglink_df[!is.na(eglink_df$start) & !is.na(eglink_df$end), , drop = FALSE]
      if (nrow(eglink_df) > 0) {
        eglink_gr <- makeGRangesFromDataFrame(
          eglink_df,
          seqnames.field = "chrom",
          start.field = "start",
          end.field = "end",
          starts.in.df.are.0based = TRUE,
          keep.extra.columns = TRUE
        )
        eglink_hits <- findOverlaps(dmr_raw, eglink_gr, ignore.strand = TRUE)
        if (length(eglink_hits) > 0) {
          eglink_genes <- data.frame(
            dmr_idx   = queryHits(eglink_hits),
            gene_id   = eglink_gr$TargetGeneEnsemblID[subjectHits(eglink_hits)],
            gene_name = eglink_gr$TargetGene[subjectHits(eglink_hits)],
            source    = "enhancer",
            stringsAsFactors = FALSE
          )
          cat(sprintf("  Enhancer-gene links: %d\n", length(eglink_hits)))
        }
      }
    }
  } else {
    cat("  Skipping enhancer-gene links: required columns not found\n")
  }
}

all_gene_links <- rbind(promoter_genes, eglink_genes)
all_gene_links <- all_gene_links[!duplicated(all_gene_links[, c("dmr_idx", "gene_name")]), ]

dmr_names <- paste0(seqnames(dmr_raw), ":", start(dmr_raw), "-", end(dmr_raw))

gene_assoc <- data.frame(
  dmr_chr    = seqnames(dmr_raw)[all_gene_links$dmr_idx],
  dmr_start  = start(dmr_raw)[all_gene_links$dmr_idx],
  dmr_end    = end(dmr_raw)[all_gene_links$dmr_idx],
  dmr_name   = dmr_names[all_gene_links$dmr_idx],
  dmr_delta  = dmr_raw$delta[all_gene_links$dmr_idx],
  dmr_type   = ifelse(dmr_raw$hyper[all_gene_links$dmr_idx], "hyper", "hypo"),
  gene_id    = all_gene_links$gene_id,
  gene_name  = all_gene_links$gene_name,
  link_source = all_gene_links$source,
  stringsAsFactors = FALSE
)

cat(sprintf("DMR-gene associations: %d rows\n", nrow(gene_assoc)))
cat("Writing:", gene_assoc_tsv, "\n")
write.table(gene_assoc, gene_assoc_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# =========================================================================
# 5. Functional enrichment with clusterProfiler
# =========================================================================
run_enrichment <- function(gene_names, title_prefix, go_tsv, go_png, kegg_tsv, kegg_png) {
  # Ensure output directory exists
  dir.create(dirname(go_tsv), recursive = TRUE, showWarnings = FALSE)

  if (length(gene_names) == 0) {
    cat(sprintf("  %s: no genes, writing empty outputs\n", title_prefix))
    write.table(data.frame(), go_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data.frame(), kegg_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    # Create empty placeholder plots
    p_empty <- ggplot(data.frame(x = 1, y = 1, label = "No enrichment data"),
                      aes(x, y, label = label)) +
      geom_text(size = 6) + theme_void()
    ggsave(go_png, p_empty, width = 10, height = 8, dpi = 150)
    ggsave(kegg_png, p_empty, width = 10, height = 8, dpi = 150)
    return(invisible())
  }

  entrez <- tryCatch({
    bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)

  if (is.null(entrez) || nrow(entrez) == 0) {
    cat(sprintf("  %s: gene ID conversion failed, writing empty outputs\n", title_prefix))
    write.table(data.frame(), go_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data.frame(), kegg_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    p_empty <- ggplot(data.frame(x = 1, y = 1, label = "Gene ID conversion failed"),
                      aes(x, y, label = label)) +
      geom_text(size = 6) + theme_void()
    ggsave(go_png, p_empty, width = 10, height = 8, dpi = 150)
    ggsave(kegg_png, p_empty, width = 10, height = 8, dpi = 150)
    return(invisible())
  }

  entrez_ids <- unique(entrez$ENTREZID)
  cat(sprintf("  %s: %d genes -> %d Entrez IDs\n", title_prefix, length(gene_names), length(entrez_ids)))

  # GO enrichment
  ego <- tryCatch({
    enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)

  if (!is.null(ego) && nrow(ego) > 0) {
    write.table(as.data.frame(ego), go_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    p <- dotplot(ego, showCategory = 20, title = paste(title_prefix, "- GO BP Enrichment"))
    ggsave(go_png, p, width = 10, height = 8, dpi = 150)
    cat(sprintf("  GO terms: %d\n", nrow(ego)))
  } else {
    cat(sprintf("  %s: no GO terms found, writing empty output\n", title_prefix))
    write.table(data.frame(), go_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    p_empty <- ggplot(data.frame(x = 1, y = 1, label = "No GO enrichment"),
                      aes(x, y, label = label)) +
      geom_text(size = 6) + theme_void()
    ggsave(go_png, p_empty, width = 10, height = 8, dpi = 150)
  }

  # KEGG enrichment
  ekegg <- tryCatch({
    enrichKEGG(gene = entrez_ids, organism = "hsa", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)

  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    write.table(as.data.frame(ekegg), kegg_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    p <- dotplot(ekegg, showCategory = 20, title = paste(title_prefix, "- KEGG Enrichment"))
    ggsave(kegg_png, p, width = 10, height = 8, dpi = 150)
    cat(sprintf("  KEGG pathways: %d\n", nrow(ekegg)))
  } else {
    cat(sprintf("  %s: no KEGG pathways found, writing empty output\n", title_prefix))
    write.table(data.frame(), kegg_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    p_empty <- ggplot(data.frame(x = 1, y = 1, label = "No KEGG enrichment"),
                      aes(x, y, label = label)) +
      geom_text(size = 6) + theme_void()
    ggsave(kegg_png, p_empty, width = 10, height = 8, dpi = 150)
  }
}

cat("\nFunctional enrichment for hyper-methylated DMR genes...\n")
hyper_genes <- unique(gene_assoc$gene_name[gene_assoc$dmr_type == "hyper"])
hyper_genes <- hyper_genes[!is.na(hyper_genes) & hyper_genes != ""]
run_enrichment(hyper_genes, "Hyper-DMR", hyper_go_tsv, hyper_go_png,
               hyper_kegg_tsv, hyper_kegg_png)

cat("\nFunctional enrichment for hypo-methylated DMR genes...\n")
hypo_genes <- unique(gene_assoc$gene_name[gene_assoc$dmr_type == "hypo"])
hypo_genes <- hypo_genes[!is.na(hypo_genes) & hypo_genes != ""]
run_enrichment(hypo_genes, "Hypo-DMR", hypo_go_tsv, hypo_go_png,
               hypo_kegg_tsv, hypo_kegg_png)

# =========================================================================
# 6. Annotation distribution plot
# =========================================================================
cat("\nPlotting annotation distribution...\n")

# Prepare plot data with log2 fold enrichment
plot_data <- ann_summary[, c("annotation_type", "n_dmrs_obs", "fold_enrich", "padj", "significance")]
plot_data$log2_fold_enrich <- log2(plot_data$fold_enrich)
plot_data$log2_fold_enrich[is.infinite(plot_data$log2_fold_enrich) | is.na(plot_data$log2_fold_enrich)] <- 0

# Create label: count + p-value
plot_data$label <- sprintf("%s", ifelse(plot_data$significance == "ns", "", plot_data$significance))

# Color by enrichment direction
plot_data$direction <- ifelse(plot_data$log2_fold_enrich > 0, "enriched",
                              ifelse(plot_data$log2_fold_enrich < 0, "depleted", "no_change"))

p <- ggplot(plot_data, aes(x = reorder(annotation_type, log2_fold_enrich), y = log2_fold_enrich, fill = direction)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = label), hjust = ifelse(plot_data$log2_fold_enrich >= 0, -0.1, 1.1),
            size = 3, nudge_x = 0) +
  scale_fill_manual(values = c("enriched" = "firebrick", "depleted" = "steelblue", "no_change" = "gray60"),
                    name = "Enrichment") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  coord_flip() +
  labs(title = "DMR Annotation Enrichment/Depletion",
       subtitle = paste0("Total DMRs: ", length(dmrs), " | Background: genome-wide"),
       x = "Annotation Type", y = expression(log[2]~Fold~Enrichment)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))

ggsave(output_png, p, width = 10, height = 6, dpi = 150)
ggsave(output_pdf, p, width = 10, height = 6)

cat("Done.\n")

