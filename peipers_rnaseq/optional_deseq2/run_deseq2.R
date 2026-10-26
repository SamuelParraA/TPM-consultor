#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(readr)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Uso: run_deseq2.R QUANT_DIR SAMPLE_METADATA TX2GENE OUTPUT_DIR")
quant_dir <- normalizePath(args[[1]], mustWork = TRUE)
metadata_path <- normalizePath(args[[2]], mustWork = TRUE)
tx2gene_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read_tsv(metadata_path, show_col_types = FALSE)
required <- c("SampleID", "genotype", "treatment", "replicate")
if (!all(required %in% colnames(metadata))) stop("Faltan columnas de metadata")
if (nrow(metadata) != 18 || anyDuplicated(metadata$SampleID)) stop("Se requieren 18 SampleID únicos")
metadata$genotype <- factor(metadata$genotype, levels = c("WT", "epi", "ACCD"))
metadata$treatment <- factor(metadata$treatment, levels = c("Control", "Salt"))
if (anyNA(metadata$genotype) || anyNA(metadata$treatment)) stop("Factores genotype/treatment inválidos")
rownames(metadata) <- metadata$SampleID

files <- file.path(quant_dir, metadata$SampleID, "quant.sf")
names(files) <- metadata$SampleID
if (!all(file.exists(files))) stop("Faltan quant.sf: ", paste(names(files)[!file.exists(files)], collapse = ", "))
tx2gene <- read_tsv(tx2gene_path, show_col_types = FALSE)
if (!all(c("TXNAME", "GENEID") %in% colnames(tx2gene))) stop("tx2gene.tsv inválido")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("TXNAME", "GENEID")])
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ genotype * treatment)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
message("Genes retenidos después del filtro: ", sum(keep))
dds <- DESeq(dds)

result_names <- resultsNames(dds)
writeLines(result_names, file.path(output_dir, "results_names.txt"))
main_treatment <- grep("^treatment_Salt_vs_Control$", result_names, value = TRUE)
epi_main <- grep("^genotype_epi_vs_WT$", result_names, value = TRUE)
accd_main <- grep("^genotype_ACCD_vs_WT$", result_names, value = TRUE)
epi_interaction <- grep("genotypeepi\\.treatmentSalt", result_names, value = TRUE)
accd_interaction <- grep("genotypeACCD\\.treatmentSalt", result_names, value = TRUE)
if (length(c(main_treatment, epi_main, accd_main, epi_interaction, accd_interaction)) != 5) {
  stop("Coeficientes inesperados: ", paste(result_names, collapse = ", "))
}

write_result <- function(result, label) {
  table <- as.data.frame(result)
  table$GeneID <- rownames(table)
  table <- table[, c("GeneID", setdiff(colnames(table), "GeneID"))]
  table <- table[order(is.na(table$padj), table$padj, table$pvalue), ]
  write_tsv(table, file.path(output_dir, paste0(label, ".tsv")), na = "NA")
}

write_result(results(dds, name = main_treatment, alpha = 0.05), "WT_Salt_vs_Control")
write_result(results(dds, contrast = list(c(main_treatment, epi_interaction)), alpha = 0.05), "epi_Salt_vs_Control")
write_result(results(dds, contrast = list(c(main_treatment, accd_interaction)), alpha = 0.05), "ACCD_Salt_vs_Control")
write_result(results(dds, name = epi_main, alpha = 0.05), "epi_vs_WT_Control")
write_result(results(dds, name = accd_main, alpha = 0.05), "ACCD_vs_WT_Control")
write_result(results(dds, contrast = list(c(epi_main, epi_interaction)), alpha = 0.05), "epi_vs_WT_Salt")
write_result(results(dds, contrast = list(c(accd_main, accd_interaction)), alpha = 0.05), "ACCD_vs_WT_Salt")
write_result(results(dds, name = epi_interaction, alpha = 0.05), "interaction_epi_vs_WT")
write_result(results(dds, name = accd_interaction, alpha = 0.05), "interaction_ACCD_vs_WT")

normalized <- as.data.frame(counts(dds, normalized = TRUE))
normalized$GeneID <- rownames(normalized)
write_tsv(normalized[, c("GeneID", metadata$SampleID)], file.path(output_dir, "normalized_counts.tsv"))

vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = c("genotype", "treatment"), returnData = TRUE)
percent_var <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color = genotype, shape = treatment, label = name)) +
  geom_point(size = 4) +
  geom_text(nudge_y = 1.5, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "%")) +
  ylab(paste0("PC2: ", percent_var[2], "%")) +
  theme_bw()
ggsave(file.path(output_dir, "PCA_samples.png"), p, width = 8, height = 6, dpi = 180)

distances <- as.matrix(dist(t(assay(vsd))))
annotation <- as.data.frame(metadata[, c("genotype", "treatment", "replicate")])
png(file.path(output_dir, "sample_distances.png"), width = 1600, height = 1400, res = 180)
pheatmap(distances, annotation_col = annotation, annotation_row = annotation,
         main = "Distancia entre bibliotecas - VST")
dev.off()

saveRDS(dds, file.path(output_dir, "deseq2_dataset.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("OK: resultados DESeq2 escritos en ", output_dir)
