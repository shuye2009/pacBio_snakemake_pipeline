#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(TFBSTools)
  library(JASPAR2024)
  library(Biostrings)
  library(GenomicRanges)
  library(rtracklayer)
  library(Rsamtools)
  library(ggplot2)
  library(parallel)
  library(RSQLite)
  library(ggseqlogo)
})

args <- commandArgs(trailingOnly = TRUE)
params <- list()
for (i in seq(1, length(args), by = 2)) {
  params[[sub("^--", "", args[i])]] <- args[i + 1]
}

dmr_bed     <- params[["dmr-bed"]]
fasta_file  <- params[["fasta"]]
jaspar_db   <- params[["jaspar-db"]]
output_tsv  <- params[["output-tsv"]]
output_png  <- params[["output-png"]]
output_pdf  <- params[["output-pdf"]]
output_html <- params[["output-html"]]
n_bg        <- as.integer(params[["n-background"]])
pval_cutoff <- as.numeric(params[["pvalue-cutoff"]])
threads     <- as.integer(params[["threads"]])

if(0){
  dmr_bed <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/visualization/significant_dmrs.bed"
  fasta_file <- "/cluster/projects/hardinggroup/Shuye/reference_genome/human/T2T_CHM13/bwameth_index/chm13v2.0_maskedMY_lambda_pUC19c.fa"
  jaspar_db <- "/cluster/home/t128737uhn/.cache/R/BiocFileCache/1e2926a17d297_JASPAR2024.sqlite"
  output_tsv <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/dmr_motif_enrichment.tsv"
  output_png <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/dmr_motif_enrichment.png"
  output_pdf <- "/cluster/projects/hardinggroup/Shuye/pacbio/output/functional_analysis/dmr_motif_enrichment.pdf"
  n_bg <- 10
  pval_cutoff <- 0.05
  threads <- 8
}

cat("Loading DMRs:", dmr_bed, "\n")
dmr_bed_df <- read.table(dmr_bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
dmr_bed_df[, 2] <- as.numeric(dmr_bed_df[, 2])
dmr_bed_df[, 3] <- as.numeric(dmr_bed_df[, 3])
dmr_gr <- makeGRangesFromDataFrame(dmr_bed_df, seqnames.field = "V1", start.field = "V2", end.field = "V3")
dmr_gr <- reduce(dmr_gr)
cat(sprintf("Loaded %d DMRs (%d after merging overlaps)\n",
            nrow(dmr_bed_df), length(dmr_gr)))

cat("Extracting sequences from:", fasta_file, "\n")
fa_file <- FaFile(fasta_file)
indexFa(fasta_file)
dmr_seqs <- getSeq(fa_file, dmr_gr)
names(dmr_seqs) <- paste0(seqnames(dmr_gr), ":", start(dmr_gr), "-", end(dmr_gr))
cat(sprintf("Extracted %d sequences, total %d bp\n", length(dmr_seqs), sum(width(dmr_seqs))))
if (sum(width(dmr_seqs)) > 50e6) {
  cat(sprintf("Warning: large total DMR sequence (%d Mb) may require significant memory.\n",
              round(sum(width(dmr_seqs)) / 1e6)))
}

cat("Loading JASPAR2024 vertebrate motifs...\n")
opts <- list(species = 9606, collection = "CORE", all_versions = FALSE)

# Try cached DB first (from prior JASPAR2024() call on internet-connected machine),
# then fall back to package extdata, or use user-supplied path.
# Run `JASPAR2024()` on login node to cache the DB: Conda activate R4.5.3 / R/4.5.3; R/ JASPAR2024() 
# The cached DB will be in "/cluster/home/t128737uhn/.cache/R/BiocFileCache/1e2926a17d297_JASPAR2024.sqlite"
if (is.null(jaspar_db) || !file.exists(jaspar_db)) {
  cache_dir <- tools::R_user_dir("BiocFileCache", "cache")
  db_file <- Sys.glob(file.path(cache_dir, "*_JASPAR2024.sqlite"))
  db_file <- db_file[1] # ensure we have one and only one DB file in the cache, if multiple exist, only keep the latest one.
  
} else {
  db_file <- jaspar_db
}
if (!file.exists(db_file)) {
  stop("JASPAR2024.sqlite not found. Download it from https://jaspar2022.genereg.net/download/database/JASPAR2024.sqlite")
}

con <- dbConnect(SQLite(), db_file)
pfm_list <- getMatrixSet(con, opts)
dbDisconnect(con)
cat(sprintf("Loaded %d motifs\n", length(pfm_list)))

cat("Converting to PWM...\n")
pwm_list <- toPWM(pfm_list, type = "log2probratio")

cat("Scanning motifs in DMR sequences...\n")
scan_motifs <- function(pwm, seqs, min_score = "90%") {
  sites <- searchSeq(pwm, seqs, seqname = "DMRs", strand = "*", min.score = min_score)
  nrow(as.data.frame(sites))
}

n_motif_cores <- min(threads, 8)
dmr_hits <- mclapply(pwm_list, scan_motifs, seqs = dmr_seqs, mc.cores = n_motif_cores)
dmr_hits <- unlist(dmr_hits)
names(dmr_hits) <- names(pwm_list)
cat(sprintf("DMR motif hits: min=%d, median=%.1f, max=%d, total=%d\n",
            min(dmr_hits), median(dmr_hits), max(dmr_hits), sum(dmr_hits)))

cat("Generating background sequences...\n")
n_bg_cores <- min(threads, n_bg, 4)
bg_hits_matrix <- do.call(cbind, mclapply(1:n_bg, function(i) {
  bg_seqs <- DNAStringSet(lapply(dmr_seqs, function(s) {
    chars <- strsplit(as.character(s), "")[[1]]
    if (length(chars) < 2) return(s)
    dinucs <- paste0(chars[-length(chars)], chars[-1])
    dinucs_shuffled <- sample(dinucs)
    result <- chars[1]
    for (j in seq_along(dinucs_shuffled)) {
      result <- paste0(result, substr(dinucs_shuffled[j], 2, 2))
    }
    DNAString(result)
  }))
  names(bg_seqs) <- names(dmr_seqs)
  bg_hits <- unlist(lapply(pwm_list, scan_motifs, seqs = bg_seqs))
  rm(bg_seqs); gc()
  bg_hits
}, mc.cores = n_bg_cores))

cat("Computing enrichment...\n")
total_dmr_bp <- sum(width(dmr_seqs))
total_bg_bp  <- total_dmr_bp

enrichment <- do.call(rbind, lapply(names(pwm_list), function(motif_name) {
  dmr_count <- dmr_hits[motif_name]
  bg_counts <- bg_hits_matrix[motif_name, ]
  bg_mean   <- mean(bg_counts)
  bg_sd     <- sd(bg_counts)

  if (bg_sd == 0) bg_sd <- 1
  if (bg_mean == 0) bg_mean <- 1

  z_score <- (dmr_count - bg_mean) / bg_sd
  p_value <- 2 * pnorm(-abs(z_score))

  fold_change <- if (bg_mean > 0) dmr_count / bg_mean else Inf

  data.frame(
    motif_id    = motif_name,
    tf_name     = name(pfm_list[[motif_name]]),
    dmr_hits    = dmr_count,
    bg_mean     = round(bg_mean, 2),
    bg_sd       = round(bg_sd, 2),
    fold_change = round(fold_change, 3),
    z_score     = round(z_score, 3),
    p_value     = signif(p_value, 4),
    stringsAsFactors = FALSE
  )
}))

enrichment <- enrichment[order(enrichment$p_value), ]
enrichment$adj_pvalue <- signif(p.adjust(enrichment$p_value, method = "BH"), 4)

cat(sprintf("Significant motifs (adj.p < %s): %d\n", pval_cutoff,
            sum(enrichment$adj_pvalue < pval_cutoff)))

cat("Writing:", output_tsv, "\n")
dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(enrichment, output_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Plotting...\n")
sig <- enrichment[enrichment$adj_pvalue < pval_cutoff, ]
if (nrow(sig) > 30) sig <- sig[1:30, ]

if (nrow(sig) > 0) {
  sig$tf_name <- factor(sig$tf_name, levels = rev(sig$tf_name))
  sig$mlog10p <- -log10(sig$adj_pvalue)
  sig$log2fc <- log2(sig$fold_change)

  p <- ggplot(sig, aes(x = log2fc, y = tf_name, size = mlog10p, color = log2fc)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(low = "blue", mid = "grey80", high = "red",
                          midpoint = 0, name = "log2(FC)") +
    scale_size_continuous(range = c(2, 8), name = "-log10(adj.P)") +
    labs(
      title = "TFBS Motif Enrichment in DMRs",
      subtitle = paste0("JASPAR2024 | Top ", nrow(sig), " significant motifs"),
      x = "log2(Fold Change vs Background)", y = ""
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
} else {
  p <- ggplot() +
    annotate("text", x = 1, y = 1, label = "No significantly enriched motifs found") +
    theme_void()
}

ggsave(output_png, p, width = 10, height = max(6, nrow(sig) * 0.25), dpi = 150)
ggsave(output_pdf, p, width = 10, height = max(6, nrow(sig) * 0.25))

cat("Generating HTML report...\n")
logo_dir <- file.path(dirname(output_html), "logos")
dir.create(logo_dir, recursive = TRUE, showWarnings = FALSE)

# Build HTML report
html <- c(
  '<!DOCTYPE html>',
  '<html><head>',
  '<meta charset="UTF-8">',
  '<title>DMR Motif Enrichment Report</title>',
  '<style>',
  'body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;',
  '       max-width: 1100px; margin: 40px auto; padding: 0 20px; color: #333; background: #fafafa; }',
  'h1 { border-bottom: 2px solid #2c3e50; padding-bottom: 10px; color: #2c3e50; }',
  'h2 { color: #2c3e50; margin-top: 30px; }',
  '.summary { background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 20px; margin: 20px 0; }',
  '.summary table { border-collapse: collapse; width: 100%; }',
  '.summary td { padding: 4px 12px; }',
  '.summary td:first-child { font-weight: 600; width: 200px; }',
  'table.motifs { width: 100%; border-collapse: collapse; margin: 20px 0; background: #fff;',
  '               border: 1px solid #ddd; border-radius: 6px; overflow: hidden; }',
  'table.motifs th { background: #2c3e50; color: #fff; padding: 10px 12px; text-align: left; font-size: 13px; }',
  'table.motifs td { padding: 8px 12px; border-bottom: 1px solid #eee; font-size: 13px; }',
  'table.motifs tr:hover { background: #f5f8fa; }',
  '.enriched { color: #c0392b; font-weight: 600; }',
  '.depleted { color: #2980b9; font-weight: 600; }',
  '.motif-card { background: #fff; border: 1px solid #ddd; border-radius: 6px;',
  '              padding: 16px; margin: 16px 0; display: flex; align-items: center; gap: 20px; }',
  '.motif-card .logo { flex-shrink: 0; }',
  '.motif-card .logo img { max-height: 80px; }',
  '.motif-card .info { flex: 1; }',
  '.motif-card .info h3 { margin: 0 0 6px 0; color: #2c3e50; }',
  '.motif-card .info .stats { font-size: 13px; color: #666; }',
  '.motif-card .info .stats span { margin-right: 16px; }',
  '</style>',
  '</head><body>',
  paste0('<h1>DMR Motif Enrichment Report</h1>'),
  '<div class="summary">',
  '<h2>Analysis Summary</h2>',
  '<table>',
  paste0('<tr><td>Total DMRs</td><td>', length(dmr_gr), '</td></tr>'),
  paste0('<tr><td>Total DMR bp</td><td>', format(sum(width(dmr_seqs)), big.mark=","), '</td></tr>'),
  paste0('<tr><td>Motifs scanned</td><td>', length(pfm_list), ' (JASPAR2024 CORE vertebrate)</td></tr>'),
  paste0('<tr><td>Background iterations</td><td>', n_bg, '</td></tr>'),
  paste0('<tr><td>Significance cutoff</td><td>adj.P &lt; ', pval_cutoff, '</td></tr>'),
  paste0('<tr><td>Significant motifs</td><td>', sum(enrichment$adj_pvalue < pval_cutoff), '</td></tr>'),
  '</table>',
  '</div>'
)

# Top 30 table
sig <- enrichment[enrichment$adj_pvalue < pval_cutoff, ]
if (nrow(sig) > 30) sig <- sig[1:30, ]

if (nrow(sig) > 0) {
  html <- c(html,
    '<h2>Top ', nrow(sig), ' Significant Motifs</h2>',
    '<table class="motifs">',
    '<tr><th>Rank</th><th>TF Name</th><th>Motif ID</th><th>DMR Hits</th><th>BG Mean</th>',
    '<th>log2(FC)</th><th>Z-score</th><th>adj.P</th></tr>'
  )
  for (i in 1:nrow(sig)) {
    fc_class <- if (sig$fold_change[i] >= 1) "enriched" else "depleted"
    html <- c(html, paste0(
      '<tr>',
      '<td>', i, '</td>',
      '<td><strong>', sig$tf_name[i], '</strong></td>',
      '<td><code>', sig$motif_id[i], '</code></td>',
      '<td>', sig$dmr_hits[i], '</td>',
      '<td>', sig$bg_mean[i], '</td>',
      '<td class="', fc_class, '">', round(log2(sig$fold_change[i]), 3), '</td>',
      '<td>', sig$z_score[i], '</td>',
      '<td>', sig$adj_pvalue[i], '</td>',
      '</tr>'
    ))
  }
  html <- c(html, '</table>')

  # Sequence logos
  html <- c(html, '<h2>Sequence Logos</h2>')
  for (i in 1:nrow(sig)) {
    motif_id <- sig$motif_id[i]
    tf_name <- sig$tf_name[i]
    logo_file <- file.path(logo_dir, paste0(motif_id, ".png"))
    logo_rel <- paste0("logos/", motif_id, ".png")

    pfm <- as.matrix(pfm_list[[motif_id]])
    p <- ggseqlogo(pfm, method = "bits") +
      labs(title = paste0(tf_name, " (", motif_id, ")")) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(size = 11, face = "bold"))
    ggsave(logo_file, p, width = 6, height = 1.8, dpi = 150)

    fc_str <- sprintf("log2(FC)=%.3f", log2(sig$fold_change[i]))
    pval_str <- sprintf("adj.P=%.2e", sig$adj_pvalue[i])
    fc_class <- if (sig$fold_change[i] >= 1) "enriched" else "depleted"
    html <- c(html, paste0(
      '<div class="motif-card">',
      '<div class="logo"><img src="', logo_rel, '" alt="', tf_name, ' logo"></div>',
      '<div class="info">',
      '<h3>', tf_name, ' <code>', motif_id, '</code></h3>',
      '<div class="stats">',
      '<span>DMR hits: <strong>', sig$dmr_hits[i], '</strong></span>',
      '<span>BG mean: <strong>', sig$bg_mean[i], '</strong></span>',
      '<span class="', fc_class, '">', fc_str, '</span>',
      '<span>', pval_str, '</span>',
      '</div>',
      '</div>',
      '</div>'
    ))
  }
} else {
  html <- c(html, '<p>No significantly enriched motifs found.</p>')
}

html <- c(html, '</body></html>')
writeLines(html, output_html)
cat("HTML report written:", output_html, "\n")

cat("Done.\n")
