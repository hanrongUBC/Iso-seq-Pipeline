# Replicate-aware differential poly(A) site usage on the re-filtered LAPA
# count matrix (count_matrix.tsv from lapa_refilter.py).
# Usage: Rscript run_dexseq_apa.R count_matrix.tsv samples.tsv out_dir [threads] [control]
# Each treatment is compared with the control using all its replicates.
# Significance = cluster FDR (DEXSeq padj) < 0.05 and |delta PAU| >= 0.1;
# gene-level FDR comes from perGeneQValue. Only clusters assigned to exactly
# one annotated gene are tested; intergenic clusters are left out.
args <- commandArgs(TRUE)
counts_file <- args[1]; samples_file <- args[2]; out <- args[3]
threads <- if (length(args) >= 4) as.integer(args[4]) else 1L
control <- if (length(args) >= 5) args[5] else "Control"
suppressPackageStartupMessages({library(DEXSeq); library(data.table)})
dir.create(out, recursive = TRUE, showWarnings = FALSE)
bp <- if (threads > 1) BiocParallel::MulticoreParam(workers = threads) else BiocParallel::SerialParam()

x <- fread(counts_file, colClasses = list(character = c("cluster_id", "gene_id", "feature")))
s <- fread(samples_file)
stopifnot(all(s$sample %in% names(x)), control %in% s$group)
x <- x[n_genes == 1]

all_results <- list()
gene_results <- list()
for (treatment in setdiff(unique(s$group), control)) {
  ss <- s[group %in% c(control, treatment)]
  cts <- as.matrix(x[, ss$sample, with = FALSE])
  rownames(cts) <- x$cluster_id
  keep <- rowSums(cts) > 0
  a <- x[keep]; cts <- cts[keep, , drop = FALSE]
  multi <- a[, .N, by = gene_id][N >= 2, gene_id]
  keep <- a$gene_id %in% multi
  a <- a[keep]; cts <- cts[keep, , drop = FALSE]
  colnames(cts) <- make.names(colnames(cts))
  cd <- data.frame(condition = factor(ss$group, levels = c(control, treatment)), row.names = colnames(cts))
  dxd <- DEXSeqDataSet(countData = cts, sampleData = cd, design = ~sample + exon + condition:exon,
                       featureID = a$cluster_id, groupID = a$gene_id)
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, BPPARAM = bp)
  dxd <- testForDEU(dxd, BPPARAM = bp)
  dxr <- DEXSeqResults(dxd)

  # Poly(A) site usage (PAU) per sample = cluster count / gene total; group means.
  totals <- rowsum(cts, a$gene_id)[a$gene_id, , drop = FALSE]
  pau <- ifelse(totals > 0, cts / totals, NA)
  is_t <- cd$condition == treatment
  res <- data.table(cluster_id = as.character(dxr$featureID), gene_id = as.character(dxr$groupID),
                    pvalue = dxr$pvalue, padj = dxr$padj)
  pau_dt <- data.table(cluster_id = rownames(cts),
                       PAU_control = rowMeans(pau[, !is_t, drop = FALSE], na.rm = TRUE),
                       PAU_treatment = rowMeans(pau[, is_t, drop = FALSE], na.rm = TRUE))
  res <- merge(res, pau_dt, by = "cluster_id")
  res <- merge(res, a[, .(cluster_id, chrom, strand, polyA_site, last_base_1based, feature)], by = "cluster_id")
  res[, delta_PAU := PAU_treatment - PAU_control]
  res[, contrast := paste0(treatment, "_vs_", control)]
  res[, significant := !is.na(padj) & padj < 0.05 & abs(delta_PAU) >= 0.1]
  q <- perGeneQValue(dxr)
  genes <- data.table(gene_id = names(q), gene_qvalue = as.numeric(q), contrast = paste0(treatment, "_vs_", control))
  genes <- merge(genes, res[, .(n_clusters = .N, n_significant_clusters = sum(significant)), by = gene_id], by = "gene_id")
  genes[, significant_gene := gene_qvalue < 0.05 & n_significant_clusters > 0]
  all_results[[treatment]] <- res
  gene_results[[treatment]] <- genes
  saveRDS(dxd, file.path(out, paste0(treatment, "_vs_", control, "_DEXSeq.rds")))
}

res <- rbindlist(all_results)
# BH over every test of every contrast, as a check on contrast-wise FDR.
res[!is.na(pvalue), padj_all_contrasts := p.adjust(pvalue, method = "BH")]
fwrite(res[order(contrast, padj)], file.path(out, "cluster_results.tsv"), sep = "\t")
genes <- rbindlist(gene_results)
fwrite(genes[order(contrast, gene_qvalue)], file.path(out, "gene_results.tsv"), sep = "\t")
summary <- genes[, .(genes_tested = .N, significant_genes = sum(significant_gene)), by = contrast]
summary <- merge(summary, res[, .(clusters_tested = .N, significant_clusters = sum(significant)), by = contrast], by = "contrast")
fwrite(summary, file.path(out, "summary.tsv"), sep = "\t")
print(summary)
capture.output(sessionInfo(), file = file.path(out, "R_session.txt"))
