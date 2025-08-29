# ================================
# DESeq2 LRT Time-course (Robust)
# ================================
# Input: a counts table where the first column is gene_id,
#        the header row contains time labels (e.g., 0,0,0,1,1,1,...),
#        with three replicates per time point.
#
# Usage (CLI, optional):
#   Rscript deseq2_lrt_timecourse.R counts.tsv out_dir 3 10 FALSE 3
#       counts.tsv   => path to counts file (TSV/CSV)
#       out_dir      => output directory (default: "deseq2_lrt_out")
#       expected_reps=> expected # replicates per time (default: 3)
#       min_row_sum  => gene prefilter threshold (default: 10)
#       use_spline   => TRUE/FALSE; if TRUE and time is numeric, uses spline model (default: FALSE)
#       spline_df    => degrees of freedom for spline (default: 3)
#
# Key outputs:
#   - DESeq2_LRT_results_by_padj.csv
#   - DESeq2_LRT_results_by_geneid.csv
#   - DESeq2_normalized_counts.tsv
#   - PCA_time.png, DispersionPlot.png
#   - colData_used.tsv, size_factors.tsv, sessionInfo.txt

# ----------------
# 0) Dependencies
# ----------------
options(warn = 1)
pkgs <- c("data.table","DESeq2","stringr","ggplot2","splines","matrixStats")
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
	  stop(paste0("Install required packages: ", paste(missing, collapse = ", ")))
}
library(data.table)
library(DESeq2)
library(stringr)
library(ggplot2)
library(splines)
library(matrixStats)

# -----------------------------
# 1) Parameters / CLI arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
counts_file    <- if (length(args) >= 1) args[1] else "counts.tsv"
output_dir     <- if (length(args) >= 2) args[2] else "deseq2_lrt_out"
expected_reps  <- if (length(args) >= 3) as.integer(args[3]) else 3
min_row_sum    <- if (length(args) >= 4) as.integer(args[4]) else 10
use_spline     <- if (length(args) >= 5) tolower(args[5]) %in% c("true","t","1","yes","y") else FALSE
spline_df      <- if (length(args) >= 6) as.integer(args[6]) else 3

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

cat("Parameters:\n",
	    "  counts_file   :", counts_file, "\n",
		    "  output_dir    :", output_dir, "\n",
			    "  expected_reps :", expected_reps, "\n",
				    "  min_row_sum   :", min_row_sum, "\n",
					    "  use_spline    :", use_spline, "\n",
						    "  spline_df     :", spline_df, "\n\n", sep = "")

# -------------------------
# 2) Read & sanitize counts
# -------------------------
# Preserve duplicate headers if times repeat (e.g., 0,0,0)
counts_dt <- fread(counts_file, check.names = FALSE)
if (ncol(counts_dt) < 4) stop("Need at least 1 gene_id column + 3 samples.")

# Enforce first column name and clean rows
setnames(counts_dt, 1, "gene_id")
counts_dt <- counts_dt[!is.na(gene_id) & gene_id != ""]
counts_dt <- unique(counts_dt, by = "gene_id")

# Sort genes alphabetically
setorder(counts_dt, gene_id)

# Build integer count matrix
gene_ids  <- counts_dt$gene_id
count_mat <- as.matrix(counts_dt[, -1, with = FALSE])
mode(count_mat) <- "numeric"
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- gene_ids

# ------------------------------------
# 3) Infer time and replicate metadata
# ------------------------------------
original_headers <- colnames(count_mat)

# Extract time tokens from headers:
#  - bare numbers (e.g., "0", "1.5")
#  - number with 'h' (e.g., "24h")
#  - prefixed "t" or "time" (e.g., "T0", "time48")
time_from_numeric_only <- stringr::str_match(original_headers, "^(\\d+(?:\\.\\d+)?)$")[,2]
time_from_timeh        <- stringr::str_match(original_headers, "(?i)(?:^|[_-])(\\d+(?:\\.\\d+)?)\\s*h?$")[,2]
time_from_t_prefix     <- stringr::str_match(original_headers, "(?i)(?:^|[_-])(?:t|time)\\s*(\\d+(?:\\.\\d+)?)")[,2]

time_char <- ifelse(!is.na(time_from_numeric_only), time_from_numeric_only,
					             ifelse(!is.na(time_from_timeh),        time_from_timeh,
										             ifelse(!is.na(time_from_t_prefix),     time_from_t_prefix,
															                    original_headers)))

# Extract replicate numbers if encoded, else assign sequentially within each time
rep_from_token <- stringr::str_match(original_headers, "(?i)(?:^|[_-])(rep|r)\\s*(\\d+)")[,3]
replicate_num  <- suppressWarnings(as.integer(rep_from_token))
if (any(is.na(replicate_num))) {
	  replicate_num <- ifelse(
							      is.na(replicate_num),
								      ave(seq_along(time_char), time_char, FUN = seq_along),
									      replicate_num
								    )
}

# Validate replicate count per time
tab_reps <- table(time_char)
if (any(tab_reps != expected_reps)) {
	  stop(paste0(
				      "Each time point must have exactly ", expected_reps, " replicates.\n",
					      "Observed per-time replicate counts: ",
						      paste(names(tab_reps), tab_reps, sep="=", collapse=", ")
					    ))
}

# Determine numeric vs. discrete time and order samples
time_numeric    <- suppressWarnings(as.numeric(time_char))
time_is_numeric <- all(!is.na(time_numeric))
if (time_is_numeric) {
	  time_levels <- sort(unique(time_numeric))
  time_factor <- factor(time_char, levels = as.character(time_levels))
    ord <- order(time_numeric, replicate_num)
} else {
	  time_levels <- sort(unique(time_char))
  time_factor <- factor(time_char, levels = time_levels)
    ord <- order(as.character(time_factor), replicate_num)
}

# Rename columns to unique sample IDs and reorder by time → replicate
sample_id <- paste0("t", time_char, "_rep", replicate_num)
count_mat <- count_mat[, ord, drop = FALSE]
time_char <- time_char[ord]
time_factor <- factor(time_char, levels = levels(time_factor))
time_numeric <- if (time_is_numeric) time_numeric[ord] else rep(NA_real_, length(ord))
replicate_num <- replicate_num[ord]
sample_id <- sample_id[ord]
colnames(count_mat) <- sample_id

# Final colData as base data.frame (no tibble rowname warnings)
coldata <- data.frame(
					    sample_id    = sample_id,
						  time_factor  = time_factor,
						    time_numeric = time_numeric,
							  replicate    = factor(replicate_num),
							    stringsAsFactors = FALSE
						)
rownames(coldata) <- coldata$sample_id

# Save colData used
fwrite(coldata, file.path(output_dir, "colData_used.tsv"), sep = "\t")

# -----------------------------
# 4) Prefilter low-count genes
# -----------------------------
keep <- rowSums(count_mat) >= min_row_sum
count_mat <- count_mat[keep, , drop = FALSE]
cat("Genes retained after prefilter (rowSums >=", min_row_sum, "):", nrow(count_mat), "\n")

# ---------------------------------------
# 5) Build DESeqDataSet and run LRT robust
# ---------------------------------------
if (use_spline && time_is_numeric) {
	  design_formula <- ~ ns(time_numeric, df = spline_df)
  reduced_formula <- ~ 1
    cat("Design: spline on numeric time with df =", spline_df, "\n")
} else {
	  design_formula <- ~ time_factor
  reduced_formula <- ~ 1
    cat("Design: discrete time (factor)\n")
}

dds <- DESeqDataSetFromMatrix(
							    countData = count_mat,
								  colData   = coldata,
								    design    = design_formula
								)

# Basic sanity checks
libsizes <- colSums(counts(dds))
if (any(libsizes == 0)) stop("One or more samples have library size 0.")
if (any(duplicated(colnames(dds)))) stop("Duplicate sample names exist after processing.")

# Helper: attempt multiple strategies, then gene-wise fallback
run_deseq_lrt_with_fallback <- function(dds, reduced) {
	  strategies <- expand.grid(
								    sfType  = c("poscounts","ratio"),
									    fitType = c("parametric","local","mean"),
										    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
									  )
  for (i in seq_len(nrow(strategies))) {
	      sfType  <- strategies$sfType[i]
      fitType <- strategies$fitType[i]
	      message(sprintf("Trying DESeq LRT with sfType='%s', fitType='%s'...", sfType, fitType))
	      ans <- try(
					       DESeq(dds, test = "LRT", reduced = reduced,
								             sfType = sfType, fitType = fitType, quiet = TRUE),
					       silent = TRUE
						       )
		      if (!inherits(ans, "try-error")) {
				        message(sprintf("Succeeded with sfType='%s', fitType='%s'.", sfType, fitType))
			        return(ans)
					    }
		    }
    message("Falling back to gene-wise dispersion estimates.")
    dds2 <- estimateSizeFactors(dds, type = "poscounts")
	  dds2 <- estimateDispersionsGeneEst(dds2, quiet = TRUE)
	  dispersions(dds2) <- mcols(dds2)$dispGeneEst
	    dds2 <- nbinomLRT(dds2, reduced = reduced)
	    dds2
}

dds <- run_deseq_lrt_with_fallback(dds, reduced = reduced_formula)

# ---------------
# 6) Save results
# ---------------
res <- results(dds)  # LRT: tests any time effect
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Sorted outputs
ord_p <- order(res_df$padj, res_df$pvalue, na.last = TRUE)
ord_g <- order(res_df$gene_id)
res_by_padj   <- res_df[ord_p, ]
res_by_geneid <- res_df[ord_g, ]

fwrite(res_by_padj,   file.path(output_dir, "DESeq2_LRT_results_by_padj.csv"))
fwrite(res_by_geneid, file.path(output_dir, "DESeq2_LRT_results_by_geneid.csv"))

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
norm_counts <- cbind(gene_id = rownames(norm_counts), as.data.frame(norm_counts))
fwrite(norm_counts, file.path(output_dir, "DESeq2_normalized_counts.tsv"), sep = "\t")

# Size factors
sf <- sizeFactors(dds)
fwrite(data.frame(sample_id = names(sf), size_factor = as.numeric(sf)),
	          file.path(output_dir, "size_factors.tsv"), sep = "\t")

# -------------
# 7) QC plots
# -------------
vsd <- vst(dds, blind = TRUE)
p_pca <- plotPCA(vsd, intgroup = "time_factor") + ggtitle("PCA: colored by time")
ggsave(file.path(output_dir, "PCA_time.png"), p_pca, width = 7, height = 6, dpi = 300)

png(file.path(output_dir, "DispersionPlot.png"), width = 1200, height = 900, res = 150)
plotDispEsts(dds)
dev.off()

# ------------------------
# 8) Diagnostics (optional)
# ------------------------
cor_mat <- cor(counts(dds, normalized = FALSE), method = "pearson")
fwrite(as.data.frame(cor_mat), file.path(output_dir, "sample_correlation.tsv"), sep = "\t")

gene_var <- matrixStats::rowVars(as.matrix(counts(dds, normalized = TRUE)))
summary_var <- summary(gene_var)
capture.output(summary_var, file = file.path(output_dir, "gene_variance_summary.txt"))

# Reproducibility
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

cat("\nDone. Outputs in: ", normalizePath(output_dir), "\n",
	    "  - DESeq2_LRT_results_by_padj.csv\n",
		    "  - DESeq2_LRT_results_by_geneid.csv\n",
			    "  - DESeq2_normalized_counts.tsv\n",
				    "  - PCA_time.png, DispersionPlot.png\n",
					    "  - colData_used.tsv, size_factors.tsv, sample_correlation.tsv\n",
						    "  - gene_variance_summary.txt, sessionInfo.txt\n", sep = "")

