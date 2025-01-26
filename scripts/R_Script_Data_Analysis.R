# ==============================================================================
# Title: RNA-Seq Analysis Workflow with DESeq2 and Functional Enrichment
# Author: Isra Boughanmi
# Date: 26/01/2025
# Description:
#   This script performs RNA-Seq data analysis for samples at different cycle 
#   days (LH+5, LH+8, and LH+11). It includes the following steps:
#     1. Data import and preprocessing of metadata and count files.
#     2. Differential expression analysis using DESeq2.
#     3. Exploratory data analysis with PCA plots, correlation heatmaps, 
#      gene count distributions etc.
#     4. Pairwise differential gene expression (DEG) analysis with MA and 
#        volcano plots.
#     5. Functional enrichment analysis with directionality visualisation.
#     6. Batch effect detection using surrogate variable analysis (SVA).

# Requirements:
#   - R version >= 4.0.0
#   - Bioconducor, DESeq2, ggplot2, pheatmap, RColorBrewer, GenomicFeatures, clusterProfiler
#   - org.Hs.eg.db, sva, dplyr, readxl, ggrepel

# Notes:
#   - Ensure the input files (metadata and count files) are in the specified 
#     directories.
# ==============================================================================

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(GenomicFeatures)# For GTF-based annotations
library(clusterProfiler) # For enrichment analysis
library(org.Hs.eg.db) # GO term mapping (necessary for functional enrichment)
library(sva)  # Surrogate variable analysis for batch effects
library(dplyr)# For data manipulation
library(readxl) # Reading Excel metadata
library(ggrepel)

set.seed(123) # For reproducibility

# ----------------------------
# Data Import and Preprocessing
# ----------------------------
# Setting working directory
setwd("/Users/u5677580/Desktop/")

# Loading metadata
sample_info <- read_excel("/Users/u5677580/Downloads/Assignment_metadata.xls", skip = 10)
colnames(sample_info) <- c("Sample_name", "Title", "Source_name", "Organism", 
                           "Cycle_day", "Live_births", "Pregnancy_losses", "Age", 
                           "BMI", "Molecule", "Description", "Processed_data_file", 
                           "Raw_file_1", "Raw_file_2")
sample_info_clean <- sample_info[grep("^Day_", sample_info$Sample_name), ]
sample_info_clean$Base_name <- gsub("_1\\.fastq", "", basename(sample_info_clean$Raw_file_1))

# Loading count data
count_data_path <- "/Users/u5677580/Desktop/gene_counts"
count_files <- list.files(count_data_path, pattern = "*.txt", full.names = TRUE)
count_data_list <- lapply(count_files, function(file) {
  data <- read.delim(file, row.names = 1, header = FALSE, check.names = FALSE)
  colname <- gsub("_counts\\.txt", "", basename(file))
  colnames(data) <- colname
  data[, 1] <- as.numeric(as.character(data[, 1]))
  return(data)
})
count_data <- do.call(cbind, count_data_list)

# Matching metadata with count data
colnames(count_data) <- sample_info_clean$Base_name
sample_info_clean <- sample_info_clean[match(colnames(count_data), sample_info_clean$Base_name), ]
sample_info_clean$Cycle_day <- factor(sample_info_clean$Cycle_day, levels = c("LH+5", "LH+8", "LH+11"))

# ---------------------------
# DESeq2 Workflow
# ---------------------------
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info_clean, design = ~ Cycle_day)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)
rld <- rlog(dds, blind = FALSE)

# ---------------------------
# Exploratory Data Analysis
# ---------------------------
# PCA Plot
plotPCA(rld, intgroup = "Cycle_day") + ggtitle("PCA Plot of Samples")

# Computing sample-to-sample correlation matrix
sampleCor <- cor(assay(rld), method = "pearson")  # Pearson correlation

# Converting the correlation matrix to a heatmap-compatible format
sampleDistMatrix <- as.matrix(sampleCor)  # Directly use the correlation matrix
colors <- colorRampPalette((brewer.pal(9, "Blues")))(255)

# Plotting heatmap: sample-to-sample correlation
pheatmap(sampleDistMatrix,
         clustering_distance_rows = as.dist(1 - sampleCor),  # Clustering based on dissimilarity
         clustering_distance_cols = as.dist(1 - sampleCor),
         col = colors,
         display_numbers = TRUE,  
         main = "Sample-to-Sample Correlation Heatmap")

# Gene Count Distribution
boxplot(log2(counts(dds) + 1), main = "Gene Count Distribution", 
        ylab = "Log2 Counts", col = "lightblue", las = 2)

# ---------------------------
# Pairwise DEG Analysis
# ---------------------------
# LH+8 vs LH+5
res_LH8_vs_LH5 <- results(dds, contrast = c("Cycle_day", "LH+8", "LH+5"), alpha = 0.05)
write.csv(as.data.frame(res_LH8_vs_LH5[order(res_LH8_vs_LH5$padj), ]), file = "DESeq2_results_LH8_vs_LH5.csv")
plotMA(res_LH8_vs_LH5, main = "MA Plot: LH+8 vs LH+5 (padj < 0.05)", ylim = c(-5, 5))

# LH+11 vs LH+5
res_LH11_vs_LH5 <- results(dds, contrast = c("Cycle_day", "LH+11", "LH+5"), alpha = 0.05)
write.csv(as.data.frame(res_LH11_vs_LH5[order(res_LH11_vs_LH5$padj), ]), file = "DESeq2_results_LH11_vs_LH5.csv")
plotMA(res_LH11_vs_LH5, main = "MA Plot: LH+11 vs LH+5 (padj < 0.05)", ylim = c(-5, 5))

# LH+11 vs LH+8
res_LH11_vs_LH8 <- results(dds, contrast = c("Cycle_day", "LH+11", "LH+8"), alpha = 0.05)
write.csv(as.data.frame(res_LH11_vs_LH8[order(res_LH11_vs_LH8$padj), ]), file = "DESeq2_results_LH11_vs_LH8.csv")
plotMA(res_LH11_vs_LH8, main = "MA Plot: LH+11 vs LH+8 (padj < 0.05)", ylim = c(-5, 5))

# ---------------------------
# More Filtering for significant genes
# ---------------------------
# For LH+8 vs LH+5
sig_genes_LH8_vs_LH5 <- res_LH8_vs_LH5[res_LH8_vs_LH5$padj < 0.05 & !is.na(res_LH8_vs_LH5$padj), ]
write.csv(as.data.frame(sig_genes_LH8_vs_LH5), file = "Significant_DEGs_LH8_vs_LH5.csv")

# For LH+11 vs LH+5
sig_genes_LH11_vs_LH5 <- res_LH11_vs_LH5[res_LH11_vs_LH5$padj < 0.05 & !is.na(res_LH11_vs_LH5$padj), ]
write.csv(as.data.frame(sig_genes_LH11_vs_LH5), file = "Significant_DEGs_LH11_vs_LH5.csv")

# For LH+11 vs LH+8
sig_genes_LH11_vs_LH8 <- res_LH11_vs_LH8[res_LH11_vs_LH8$padj < 0.05 & !is.na(res_LH11_vs_LH8$padj), ]
write.csv(as.data.frame(sig_genes_LH11_vs_LH8), file = "Significant_DEGs_LH11_vs_LH8.csv")
# ---------------------------
# Volcano plot for DEG Analysis
# ---------------------------
create_volcano_plot <- function(res, title) {
  res_df <- as.data.frame(res)
  res_df$Significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                                "Significant", "Not Significant")
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_text_repel(data = subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1), 
                    aes(label = rownames(subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1))),
                    size = 3, max.overlaps = 10) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", color = "Significance") +
    theme(legend.position = "top")
}

# Volcano plot for LH+8 vs LH+5
volcano_LH8_vs_LH5 <- create_volcano_plot(res_LH8_vs_LH5, "Volcano Plot: LH+8 vs LH+5")
print(volcano_LH8_vs_LH5)

# Volcano plot for LH+11 vs LH+5
volcano_LH11_vs_LH5 <- create_volcano_plot(res_LH11_vs_LH5, "Volcano Plot: LH+11 vs LH+5")
print(volcano_LH11_vs_LH5)

# Volcano plot for LH+11 vs LH+8
volcano_LH11_vs_LH8 <- create_volcano_plot(res_LH11_vs_LH8, "Volcano Plot: LH+11 vs LH+8")
print(volcano_LH11_vs_LH8)

# ---------------------------
# XY Scatter Plot with Log2 Fold change comparison
# ---------------------------
create_xy_plot <- function(res1, res2, title) {
  res1_df <- as.data.frame(res1)
  res2_df <- as.data.frame(res2)
  common_genes <- intersect(rownames(res1_df), rownames(res2_df))
  res1_df <- res1_df[common_genes, ]
  res2_df <- res2_df[common_genes, ]
  
  xy_df <- data.frame(Log2FC1 = res1_df$log2FoldChange, Log2FC2 = res2_df$log2FoldChange)
  
  correlation <- cor(xy_df$Log2FC1, xy_df$Log2FC2, use = "complete.obs")
  
  ggplot(xy_df, aes(x = Log2FC1, y = Log2FC2)) +
    geom_point(alpha = 0.8, size = 1.5, color = "blue") +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    annotate("text", x = max(xy_df$Log2FC1, na.rm = TRUE) * 0.7, 
             y = min(xy_df$Log2FC2, na.rm = TRUE) * 0.7, 
             label = paste("Correlation:", round(correlation, 2)),
             size = 4, color = "black", hjust = 0) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change (Condition 1)", 
         y = "Log2 Fold Change (Condition 2)") +
    theme(legend.position = "none")
}

# XY scatter plot for LH+8 vs LH+5 and LH+11 vs LH+5
xy_LH8_vs_LH5_LH11_vs_LH5 <- create_xy_plot(res_LH8_vs_LH5, res_LH11_vs_LH5, 
                                            "XY Plot: LH+8 vs LH+5 vs LH+11 vs LH+5")
print(xy_LH8_vs_LH5_LH11_vs_LH5)

# XY scatter plot for LH+11 vs LH+5 and LH+11 vs LH+8
xy_LH11_vs_LH5_LH11_vs_LH8 <- create_xy_plot(res_LH11_vs_LH5, res_LH11_vs_LH8, 
                                             "XY Plot: LH+11 vs LH+5 vs LH+11 vs LH+8")
print(xy_LH11_vs_LH5_LH11_vs_LH8)


# ---------------------------
# Functional Enrichment Analysis with visualisation of directionality for LH+8 vs LH+5
# ---------------------------

# Explanation: GTF files provide gene annotations, but functional annotations like GO terms require databases such as org.Hs.eg.db.
# org.Hs.eg.db contains curated mappings between genes (e.g., Entrez IDs) and their biological functions.

# Extract significant genes for enrichment analysis (LH+8 vs LH+5)
sig_genes <- rownames(res_LH8_vs_LH5[which(res_LH8_vs_LH5$padj < 0.05), ])

# Convert gene symbols to Entrez IDs using org.Hs.eg.db
entrez_ids <- mapIds(org.Hs.eg.db, keys = sig_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Perform GO Enrichment Analysis
ego_LH8_vs_LH5 <- enrichGO(gene = na.omit(entrez_ids), OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", readable = TRUE)

# Calculating directionality (Upregulated/Downregulated/Mixed) for each GO term
ego_LH8_vs_LH5@result$Directionality <- sapply(ego_LH8_vs_LH5@result$Description, function(go_term) {
  go_genes <- ego_LH8_vs_LH5@result[ego_LH8_vs_LH5@result$Description == go_term, "geneID"]
  go_genes_list <- strsplit(go_genes, "/")[[1]]
  go_genes_de <- res_LH8_vs_LH5[rownames(res_LH8_vs_LH5) %in% go_genes_list, ]
  upregulated <- sum(go_genes_de$log2FoldChange > 0, na.rm = TRUE)
  downregulated <- sum(go_genes_de$log2FoldChange < 0, na.rm = TRUE)
  if (upregulated > downregulated) return("Upregulated")
  if (downregulated > upregulated) return("Downregulated")
  return("Mixed")
})

# Colours based on directionality
bar_colors <- ifelse(ego_LH8_vs_LH5@result$Directionality == "Upregulated", "lightblue",
                     ifelse(ego_LH8_vs_LH5@result$Directionality == "Downregulated", "pink","grey"))


# Margin adjusted for labels
par(mar = c(10, 4, 4, 2))  

# Creating the bar plot without axis labels
bp <- barplot(ego_LH8_vs_LH5@result$Count[1:10], 
              names.arg = NA,  # Suppress x-axis labels
              col = bar_colors[1:10],
              main = "Top 10 GO Biological Processes with Directionality (LH+8 vs LH+5)",
              ylab = "Gene Count",
              cex.names = 0.8)

text(x = bp, y = par("usr")[3] - 0.5,  # Position labels just below the axis
     labels = ego_LH8_vs_LH5@result$Description[1:10], 
     srt = 45,  # Set rotation angle to 45 degrees
     adj = 1,   # Right-align the labels
     xpd = TRUE, # Allow text to plot outside plot region
     cex = 0.8)  # Adjust label size

# Legend
legend("topright", legend = c("Upregulated", "Downregulated"),
       fill = c("lightblue", "pink"), cex = 0.8)

# ---------------------------
# Functional Enrichment Analysis with Directionality Visualization for LH+5 vs LH+11
# ---------------------------

# Extracting significant genes for enrichment analysis (LH+5 vs LH+11)
sig_genes_LH5_vs_LH11 <- rownames(res_LH11_vs_LH5[which(res_LH11_vs_LH5$padj < 0.05), ])

# Converting gene symbols to Entrez IDs
entrez_ids_LH5_vs_LH11 <- mapIds(org.Hs.eg.db, keys = sig_genes_LH5_vs_LH11, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

#  GO Enrichment Analysis
ego_LH5_vs_LH11 <- enrichGO(gene = na.omit(entrez_ids_LH5_vs_LH11), OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", readable = TRUE)

# Calculating directionality as well
ego_LH5_vs_LH11@result$Directionality <- sapply(ego_LH5_vs_LH11@result$Description, function(go_term) {
  go_genes <- ego_LH5_vs_LH11@result[ego_LH5_vs_LH11@result$Description == go_term, "geneID"]
  go_genes_list <- strsplit(go_genes, "/")[[1]]
  go_genes_de <- res_LH11_vs_LH5[rownames(res_LH11_vs_LH5) %in% go_genes_list, ]
  upregulated <- sum(go_genes_de$log2FoldChange > 0, na.rm = TRUE)
  downregulated <- sum(go_genes_de$log2FoldChange < 0, na.rm = TRUE)
  if (upregulated > downregulated) return("Upregulated")
  if (downregulated > upregulated) return("Downregulated")
  return("Mixed")
})

bar_colors <- ifelse(ego_LH8_vs_LH5@result$Directionality == "Upregulated", "lightblue",
                     ifelse(ego_LH8_vs_LH5@result$Directionality == "Downregulated", "pink", "grey"))

par(mar = c(10, 6, 4, 2))  # Adjust bottom margin for angled labels

bp <- barplot(ego_LH5_vs_LH11@result$Count[1:10], 
              names.arg = NA,  # Suppress x-axis labels
              col = bar_colors[1:10],
              main = "Top 10 GO Biological Processes with Directionality (LH+5 vs LH+11)",
              ylab = "Gene Count",
              cex.names = 0.8)
text(x = bp, y = par("usr")[3] - 0.5,  # Position labels just below the axis
     labels = ego_LH5_vs_LH11@result$Description[1:10], 
     srt = 45,  # Set rotation angle to 45 degrees
     adj = 1,   # Right-align the labels
     xpd = TRUE, # Allow text to plot outside plot region
     cex = 0.6)  # Adjust label size
legend("topright", legend = c("Upregulated", "Downregulated", "Mixed"),
       fill = c("lightblue", "pink", "grey"), cex = 0.8)

# ---------------------------
# Batch Effect Detection 
# ---------------------------
mod <- model.matrix(~ Cycle_day, colData(dds))
sv <- sva(assay(dds), mod)

