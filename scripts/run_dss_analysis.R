#!/usr/bin/env Rscript
#
# DSS differential methylation analysis: DML test + DMR calling
#
# Input:  pb-CpG-tools combined.bed.gz files (one per sample)
# Output: DML and DMR tables (TSV + BED)

suppressPackageStartupMessages({
  library(DSS)
  library(data.table)
  library(argparse)
})


load_sample <- function(bed_path, sample_name) {
  # pb-CpG-tools columns: chrom, begin, end, mod_score, type, cov,
  #   est_mod_count, est_unmod_count, discretized_mod_score
  dt <- fread(
    cmd = sprintf("zcat %s | grep -v '^#'", bed_path),
    header = FALSE,
    select = c(1, 2, 6, 7),
    col.names = c("chr", "pos", "N", "X"),
    showProgress = FALSE
  )
  dt[, pos := pos + 1L]          # BED 0-based -> DSS 1-based
  dt <- dt[N >= 1]
  dt[X > N, X := N]
  setnames(dt, old = c("N", "X"),
           new = c(sprintf("N_%s", sample_name),
                   sprintf("X_%s", sample_name)))
  setkey(dt, chr, pos)
  dt[]
}


merge_samples <- function(sample_list) {
  Reduce(function(x, y) merge(x, y, by = c("chr", "pos"), all = TRUE),
         sample_list)
}


build_bsseq <- function(merged, sample_names) {
  N_cols <- paste0("N_", sample_names)
  X_cols <- paste0("X_", sample_names)
  for (col in c(N_cols, X_cols)) {
    set(merged, which(is.na(merged[[col]])), col, 0L)
  }
  N_mat <- as.matrix(merged[, ..N_cols])
  X_mat <- as.matrix(merged[, ..X_cols])
  storage.mode(N_mat) <- "integer"
  storage.mode(X_mat) <- "integer"
  colnames(N_mat) <- sample_names
  colnames(X_mat) <- sample_names
  BSseq(M = X_mat, Cov = N_mat, chr = merged$chr, pos = merged$pos,
        sampleNames = sample_names)
}


main <- function() {
  parser <- ArgumentParser(description = "DSS DML/DMR analysis")
  parser$add_argument("--beds",            required = TRUE)
  parser$add_argument("--samples",         required = TRUE)
  parser$add_argument("--groups",          required = TRUE)
  parser$add_argument("--output-dml-tsv",  required = TRUE)
  parser$add_argument("--output-dml-bed",  required = TRUE)
  parser$add_argument("--output-dmr-tsv",  required = TRUE)
  parser$add_argument("--output-dmr-bed",  required = TRUE)
  parser$add_argument("--delta",       type = "double",  default = 0.1)
  parser$add_argument("--pvalue",      type = "double",  default = 1e-5)
  parser$add_argument("--fdr",         type = "double",  default = 0.05)
  parser$add_argument("--min-cpg",     type = "integer", default = 3)
  parser$add_argument("--min-length",  type = "integer", default = 50)
  parser$add_argument("--merge-dist",  type = "integer", default = 100)
  parser$add_argument("--smoothing",   type = "logical", default = TRUE)
  parser$add_argument("--smoothing-span", type = "integer", default = 500)
  parser$add_argument("--threads",     type = "integer", default = 1)
  args <- parser$parse_args()

  beds    <- strsplit(args$beds,    ",")[[1]]
  samples <- strsplit(args$samples, ",")[[1]]
  groups  <- strsplit(args$groups,  ",")[[1]]

  stopifnot(length(beds) == length(samples),
            length(samples) == length(groups))

  cat(sprintf("Loading %d samples...\n", length(samples)))

  sample_list <- mapply(
    function(b, s) load_sample(b, s),
    beds, samples, SIMPLIFY = FALSE, USE.NAMES = FALSE
  )

  cat("Merging CpG sites across samples...\n")
  merged <- merge_samples(sample_list)

  cat(sprintf("Merged %d CpG sites.\n", nrow(merged)))

  cat("Building BSseq object...\n")
  bs <- build_bsseq(merged, samples)

  cat(sprintf("Running DMLtest (smoothing=%s, span=%d, ncores=%d)...\n",
              args$smoothing, args$smoothing_span, args$threads))
  dml_test <- DMLtest(
    bs,
    group1 = samples[groups == "case"],
    group2 = samples[groups == "control"],
    smoothing = args$smoothing,
    smoothing.span = args$smoothing_span,
    ncores = args$threads
  )

  cat("Calling DMLs...\n")
  dml <- callDML(
    dml_test,
    delta = args$delta,
    p.threshold = args$pvalue
  )
  cat(sprintf("Called %d DMLs.\n", nrow(dml)))

  cat("Calling DMRs...\n")
  dmr <- callDMR(
    dml_test,
    delta = args$delta,
    p.threshold = args$pvalue,
    minlen = args$min_length,
    minCG = args$min_cpg,
    dis.merge = args$merge_dist,
    pct.sig = 0.5
  )
  cat(sprintf("Called %d DMRs.\n", nrow(dmr)))

  # --- Write DML outputs ---
  dml_out <- as.data.table(dml)
  setorder(dml_out, chr, pos)
  dml_out[, pos := as.integer(pos)]
  fwrite(dml_out, args$output_dml_tsv, sep = "\t")

  # DML BED: 0-based, score = -log10(pval) scaled to 0-1000
  dml_bed <- copy(dml_out)
  dml_bed[, `:=`(
    start  = pos - 1L,
    end    = pos,
    name   = sprintf("%s:%d", chr, pos),
    score  = as.integer(pmin(pmax(-log10(pval) * 100, 0), 1000)),
    strand = "."
  )]
  dml_bed[, `:=`(start = as.integer(start), end = as.integer(end))]
  setcolorder(dml_bed, c("chr", "start", "end", "name", "score", "strand"))
  fwrite(dml_bed, args$output_dml_bed, sep = "\t", col.names = FALSE)

  # --- Write DMR outputs ---
  dmr_out <- as.data.table(dmr)
  if (nrow(dmr_out) > 0) {
    setorder(dmr_out, chr, start)
    dmr_out[, `:=`(start = as.integer(start), end = as.integer(end))]
    fwrite(dmr_out, args$output_dmr_tsv, sep = "\t")

    dmr_bed <- dmr_out[, .(
      chr, start = start - 1L, end,
      name   = sprintf("%s:%d-%d", chr, start, end),
      score  = as.integer(pmin(pmax(-log10(areaStat) * 100, 0), 1000)),
      strand = "."
    )]
    dmr_bed[, `:=`(start = as.integer(start), end = as.integer(end))]
    fwrite(dmr_bed, args$output_dmr_bed, sep = "\t", col.names = FALSE)
  } else {
    fwrite(dmr_out, args$output_dmr_tsv, sep = "\t")
    file.create(args$output_dmr_bed)
  }

  cat("Done.\n")
}

if (sys.nframe() == 0) main()
