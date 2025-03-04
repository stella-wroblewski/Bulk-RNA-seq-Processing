# Bulk RNA-seq Data Analysis Pipeline
# Author: Your Name
# Date: YYYY-MM-DD

###############################################################################
# This script outlines a generalized workflow for bulk RNA-seq data processing.
# The pipeline includes:
#   1) Data loading and cleaning
#   2) DESeq2-based differential expression analysis
#   3) PCA visualization (overall and condition/time-specific)
#   4) Volcano plots
#   5) Heatmap of key genes
#   6) KEGG pathway enrichment
#
# Users should adjust:
#   - File paths, conditions, and sample grouping
#   - Desired thresholds for significance
#   - Gene lists for specialized heatmaps, etc.
#
# This script is intended to showcase a robust analysis pipeline for publication.
###############################################################################

############################################
# Step 0: Load Required Libraries
############################################

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(openxlsx)
library(FactoMineR)
library(factoextra)
library(DESeq2)
library(tibble)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

############################################
# Step 1: Define Input File Paths
############################################

# Replace these with your actual data file paths. Each file should contain:
#   - A 'gene' column
#   - Columns representing samples (e.g., expression counts)
file_1 <- "Path/To/DataFile_1.xlsx"
file_2 <- "Path/To/DataFile_2.xlsx"
file_3 <- "Path/To/DataFile_3.xlsx"

# Map your files to condition/time labels:
files <- list(
  "TimePoint_1" = file_1,
  "TimePoint_2" = file_2,
  "TimePoint_3" = file_3
  # Add more as needed
)

############################################
# Step 2: Data Loading Function
############################################

# This function:
#   - Reads data from an Excel file
#   - Removes any unneeded statistical columns
#   - Pivots into a long format with:
#       gene, sample, expression, time/condition, group
#   - Adjust the group detection logic as needed

load_data <- function(file, label) {
  cat("Loading data for", label, "...\n")
  data <- read.xlsx(file)
  cat("Data dimensions for", label, ":", dim(data), "\n")

  # Remove columns matching certain patterns if they exist
  data <- data %>% select(-matches("wald\\.test|lrt\\.test"), everything())

  # Pivot data into a long format
  data_long <- data %>%
    pivot_longer(
      cols = -gene,
      names_to = "sample",
      values_to = "expression"
    ) %>%
    mutate(
      condition = label,
      # Modify logic to define groups (e.g., treated vs control) if relevant
      group = ifelse(grepl("treated", sample, ignore.case = TRUE), "Treated", "Control")
    )

  return(data_long)
}

############################################
# Step 3: Load & Combine All Data
############################################

combined_data_long <- bind_rows(
  lapply(names(files), function(x) load_data(files[[x]], x))
)

# Remove duplicate (gene, sample) entries by taking the mean expression
combined_data_long <- combined_data_long %>%
  group_by(gene, sample) %>%
  summarise(expression = mean(expression), .groups = "drop") %>%
  left_join(
    combined_data_long %>% distinct(gene, condition, sample, group),
    by = c("gene", "sample")
  )

############################################
# Step 4: Prepare Counts & Metadata for DESeq2
############################################

cat("Preparing data for DESeq2...\n")

# Reshape to wide format for count matrix
count_matrix <- combined_data_long %>%
  select(gene, sample, expression) %>%
  pivot_wider(names_from = sample, values_from = expression, values_fill = 0)

# Set row names to the 'gene' column, remove 'gene'
count_matrix <- as.data.frame(count_matrix)
rownames(count_matrix) <- count_matrix$gene
count_matrix$gene <- NULL

# Round expression to integers (DESeq2 expects integer counts)
count_matrix <- round(as.matrix(count_matrix))

# Build metadata (sample annotations)
metadata <- combined_data_long %>%
  select(sample, group, condition) %>%
  distinct() %>%
  arrange(match(sample, colnames(count_matrix)))

# Ensure alignment of metadata samples and count matrix columns
if(!all(metadata$sample == colnames(count_matrix))) {
  stop("Sample names in metadata and count matrix do not match.")
}

# Use sample names as row names
rownames(metadata) <- metadata$sample

# Convert condition to factor if needed
metadata$condition <- factor(metadata$condition, levels = names(files))

############################################
# Step 5: Create DESeqDataSet & Run DESeq2
############################################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  # Example design: condition + group (modify as needed)
  design = ~ condition + group
)

cat("Running DESeq2...\n")
dds <- DESeq(dds)

# Extract normalized counts for downstream analysis
cat("Extracting normalized counts...\n")
vsd <- vst(dds, blind = FALSE)
normalized_counts <- assay(vsd)

############################################
# Step 6: PCA (Overall)
############################################

cat("Performing overall PCA...\n")

pca_data <- t(normalized_counts)
pca_data <- as.data.frame(pca_data)
pca_data$sample <- rownames(pca_data)

# Merge with metadata
pca_data <- left_join(pca_data, metadata, by = "sample")

# PCA with FactoMineR
pca_result <- PCA(pca_data %>% select(-sample, -group, -condition), graph = FALSE)

# Extract PCA results
pca_scores <- as.data.frame(pca_result$ind$coord)
pca_scores$sample <- rownames(pca_scores)
pca_scores <- left_join(pca_scores, metadata, by = "sample")

# Variance explained
variance <- pca_result$eig[, 2]
variance <- variance / sum(variance) * 100

# Plot PCA
pca_plot <- ggplot(pca_scores, aes(
  x = Dim.1,
  y = Dim.2,
  color = group,
  shape = condition
)) +
  geom_point(size = 5) +
  labs(
    title = "PCA of Transcriptomes",
    x = paste0("PC1: ", round(variance[1], 2), "% variance"),
    y = paste0("PC2: ", round(variance[2], 2), "% variance")
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

print(pca_plot)
ggsave("pca_plot_overall.png", plot = pca_plot, width = 10, height = 8)

############################################
# Step 7: Condition-Specific PCA (Optional)
############################################

cat("Performing condition/time-specific PCA...\n")
unique_conditions <- levels(metadata$condition)

for (cond in unique_conditions) {
  cat(paste("  - PCA for", cond, "\n"))

  # Subset normalized counts to this condition's samples
  selected_samples <- metadata %>% filter(condition == cond) %>% pull(sample)
  data_subset <- normalized_counts[, selected_samples]

  pca_sub <- t(data_subset)
  pca_sub <- as.data.frame(pca_sub)
  pca_sub$sample <- rownames(pca_sub)
  pca_sub <- left_join(pca_sub, metadata, by = "sample")

  pca_result_sub <- PCA(pca_sub %>% select(-sample, -group, -condition), graph = FALSE)
  pca_scores_sub <- as.data.frame(pca_result_sub$ind$coord)
  pca_scores_sub$sample <- rownames(pca_scores_sub)
  pca_scores_sub <- left_join(pca_scores_sub, metadata, by = "sample")

  variance_sub <- pca_result_sub$eig[, 2]
  variance_sub <- variance_sub / sum(variance_sub) * 100

  pca_plot_sub <- ggplot(pca_scores_sub, aes(
    x = Dim.1,
    y = Dim.2,
    color = group
  )) +
    geom_point(size = 5) +
    labs(
      title = paste0("PCA - ", cond),
      x = paste0("PC1: ", round(variance_sub[1], 2), "% variance"),
      y = paste0("PC2: ", round(variance_sub[2], 2), "% variance")
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )

  print(pca_plot_sub)
  ggsave(paste0("pca_plot_", cond, ".png"), plot = pca_plot_sub, width = 10, height = 8)
}

############################################
# Step 8: Differential Expression Analysis
############################################

cat("Running differential expression analysis...\n")

differential_expression_results <- list()

for (cond in unique_conditions) {
  cat(paste("  - DE analysis for", cond, "\n"))

  dds_cond <- dds[, dds$condition == cond]

  # Check if at least two groups are present in this subset
  if (length(unique(dds_cond$group)) < 2) {
    warning(paste("Not enough groups at", cond))
    next
  }

  # Re-level the group factor to ensure "Control" is the reference (adjust if needed)
  dds_cond$group <- droplevels(dds_cond$group)
  dds_cond$group <- relevel(dds_cond$group, ref = "Control")

  # Redefine design if needed
  design(dds_cond) <- ~ group

  dds_cond <- DESeq(dds_cond)

  # Contrast: Treated vs Control (modify as needed)
  res <- results(dds_cond, contrast = c("group", "Treated", "Control"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$condition <- cond

  differential_expression_results[[cond]] <- res_df

  # Volcano Plot
  padj_cutoff <- 0.05
  logfc_cutoff <- 1

  res_df <- res_df %>%
    mutate(Significance = case_when(
      padj < padj_cutoff & log2FoldChange >= logfc_cutoff ~ "Upregulated",
      padj < padj_cutoff & log2FoldChange <= -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))

  res_df$Significance <- factor(
    res_df$Significance,
    levels = c("Upregulated", "Downregulated", "Not Significant")
  )

  volcano <- ggplot(res_df, aes(
    x = log2FoldChange,
    y = -log10(pvalue),
    color = Significance
  )) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")
    ) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(
      title = paste0("Volcano Plot - ", cond),
      x = "Log2 Fold Change",
      y = "-Log10(p-value)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )

  print(volcano)
  ggsave(paste0("volcano_plot_", cond, ".png"), plot = volcano, width = 8, height = 6)
}

# Combine DE results into one data frame
all_de_results <- bind_rows(differential_expression_results)

# Reorder columns if desired
all_de_results <- all_de_results %>%
  select(condition, gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# Save final DE results
output_file <- "differential_expression_results_DESeq2.xlsx"
cat("Saving DE results to", output_file, "\n")
write.xlsx(all_de_results, output_file, rowNames = FALSE)

############################################
# Step 9: Heatmap of Selected Genes (Optional)
############################################

# Example code to create a heatmap of (Treated - Control) differences
# for selected genes of interest. Adjust gene list, color scale, etc.

genes_of_interest <- c("GENE1", "GENE2") # Modify as needed

# Convert normalized counts to data frame
counts_df <- as.data.frame(normalized_counts)
counts_df$gene <- rownames(counts_df)

# Merge counts with metadata
counts_long <- counts_df %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  left_join(metadata, by = "sample")

# Filter genes
counts_long_filtered <- counts_long %>% filter(gene %in% genes_of_interest)

# Average expression by (gene, condition, group)
counts_summarized <- counts_long_filtered %>%
  group_by(gene, condition, group) %>%
  summarise(mean_expression = mean(expression), .groups = "drop")

# Pivot so columns = group
counts_wide <- counts_summarized %>%
  pivot_wider(names_from = group, values_from = mean_expression)

# Compute difference (Treated - Control)
counts_wide <- counts_wide %>%
  mutate(difference = Treated - Control)

# Pivot so rows = genes, columns = condition
heatmap_data <- counts_wide %>%
  select(gene, condition, difference) %>%
  pivot_wider(names_from = condition, values_from = difference)

# Create matrix
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- heatmap_data$gene

# Define color scale centered at 0 (white)
max_diff <- max(abs(heatmap_matrix), na.rm = TRUE)
my_breaks <- seq(-max_diff, max_diff, length.out = 101)
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)

pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = my_colors,
  breaks = my_breaks,
  main = "Treated - Control Differences"
)

############################################
# Step 10: KEGG Pathway Enrichment (Optional)
############################################

# Example pipeline for KEGG enrichment analysis using clusterProfiler.
# Adjust p-value cutoff, organism DB, etc.

kegg_results_list <- list()

for (cond in unique(all_de_results$condition)) {
  cat(paste("Performing KEGG for", cond, "...\n"))

  subset_de <- all_de_results %>% filter(condition == cond)

  # Filter significant genes
  sig_genes <- subset_de %>% filter(padj < 0.05) %>% drop_na(log2FoldChange)

  if (nrow(sig_genes) > 0) {
    gene_map <- bitr(sig_genes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    gene_map <- gene_map[!duplicated(gene_map$ENTREZID), ]
    gene_map <- gene_map %>% left_join(sig_genes, by = c("SYMBOL" = "gene"))

    fold_changes <- gene_map$log2FoldChange
    names(fold_changes) <- gene_map$ENTREZID

    kegg_result <- enrichKEGG(
      gene          = names(fold_changes),
      organism      = "hsa",  # Human
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )

    if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
      kegg_df <- kegg_result@result %>%
        mutate(
          GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
          NumGenes  = sapply(strsplit(geneID, "/"), length)
        )
      kegg_results_list[[cond]] <- kegg_df
    } else {
      cat("No enriched KEGG pathways at", cond, "\n")
    }
  } else {
    cat("No significant genes at", cond, "\n")
  }
}

# Save KEGG results
if (length(kegg_results_list) > 0) {
  write.xlsx(kegg_results_list, file = "KEGG_Pathway_Analysis.xlsx", rowNames = FALSE)
}

# Example Dot Plot for top KEGG pathways
if (length(kegg_results_list) > 0) {
  # Determine global color scale for -log10 p-values
  pval_range <- range(
    unlist(lapply(kegg_results_list, function(df) if (!is.null(df)) -log10(df$p.adjust))),
    na.rm = TRUE
  )

  for (cond in names(kegg_results_list)) {
    if (!is.null(kegg_results_list[[cond]])) {
      df_kegg <- kegg_results_list[[cond]] %>% head(10)

      dotplot <- ggplot(df_kegg, aes(
        x = GeneRatio,
        y = reorder(Description, -p.adjust),
        size = NumGenes,
        color = -log10(p.adjust)
      )) +
        geom_point(alpha = 0.8) +
        scale_size_continuous(range = c(3, 12), name = "Gene Count") +
        scale_color_gradient(low = "blue", high = "red", limits = pval_range) +
        labs(
          title = paste0("KEGG - ", cond),
          x = "Gene Ratio",
          y = "Pathway",
          color = "-Log10(p-value)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid.major = element_line(color = "gray80", size = 0.3, linetype = "solid"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14, face = "bold")
        )

      print(dotplot)
      ggsave(filename = paste0("kegg_dotplot_", cond, ".png"), plot = dotplot, width = 10, height = 8)
    }
  }
}

cat("All analyses completed!\n")
