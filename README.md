# Detailed README for Bulk RNA-seq Data Analysis Pipeline

## Overview
This pipeline provides a structured approach to analyzing bulk RNA sequencing (RNA-seq) data, from raw expression counts to functional annotation using KEGG pathways. The workflow is organized into several steps, ensuring a clear and reproducible process for identifying differentially expressed genes and higher-level biological insights.

## Pipeline Steps

### 1. Data Loading and Cleaning
- **Input**: Excel files containing gene expression data for each time point or condition.
  - Each file contains a `gene` column and additional columns representing samples.
- **Actions**:
  1. Remove unwanted statistical columns if they exist (e.g., Wald or LRT test columns from previous analyses).
  2. Transform the data from a wide format to a long format using `tidyr::pivot_longer`.
  3. Tag each sample with a `condition` label (e.g., time point or treatment condition) and a `group` label (e.g., Treated vs Control).
- **Output**: A combined data frame (`combined_data_long`) in long format, containing expression data for all genes across all samples and conditions.

### 2. Data Integration
- **Process**:
  - Concatenate data from multiple conditions using `dplyr::bind_rows`.
  - Remove duplicate `(gene, sample)` entries by taking the mean of the expression values.
  - Prepare a final merged data set.

### 3. DESeq2 Data Preparation
- **Process**:
  1. Convert the long-format data into a wide format to produce a raw count matrix.
  2. Round expression values to integers (required by DESeq2).
  3. Build a metadata (colData) data frame containing `sample`, `group`, and `condition`.
  4. Ensure the order of samples in the metadata matches the columns of the count matrix.

### 4. Differential Expression Analysis (DESeq2)
- **Process**:
  1. Create a `DESeqDataSet` from the count matrix and the metadata (`colData`).
  2. Specify the design formula (e.g., `~ condition + group`).
  3. Run the DESeq pipeline (`DESeq(dds)`).
  4. Extract normalized counts using the variance-stabilizing transformation (VST) for downstream analyses (`vst(dds, blind=FALSE)`).

### 5. Principal Component Analysis (PCA)
- **Overall PCA**: Analyze the entire data set across all conditions/time points.
- **Condition-Specific PCA**: Conduct PCA separately for each condition/time point.
- **Tools**: `FactoMineR` and `factoextra` are used for PCA, while `ggplot2` visualizes the PCA results.
- **Output**: PCA plots showing the variation among samples and highlighting differences between groups.

### 6. Volcano Plots
- **Purpose**: Visualize differentially expressed genes (DEGs) for each condition/time point.
- **Process**:
  - For each subset (condition), run DESeq2.
  - Define the contrast of interest (e.g., `Treated` vs `Control`).
  - Mark significance based on thresholds (e.g., `padj < 0.05` and `|log2FoldChange| > 1`).
  - Plot gene distribution using `ggplot2` with color-coded points for upregulated, downregulated, or not significant genes.

### 7. Heatmap of Key Genes
- **Purpose**: Provide a focused view of gene expression changes for a set of key or significant genes.
- **Process**:
  1. Extract normalized counts.
  2. Filter for genes of interest (e.g., biomarkers or top DEGs).
  3. Average expression per condition/group.
  4. Create a matrix of expression differences (e.g., Treated – Control).
  5. Use `pheatmap` or similar tools to generate a heatmap.

### 8. KEGG Pathway Enrichment
- **Purpose**: Explore the biological pathways enriched among the DEGs.
- **Process**:
  1. Identify significant genes (e.g., `padj < 0.05`).
  2. Convert gene symbols to Entrez IDs using `clusterProfiler::bitr`.
  3. Run KEGG enrichment analysis (`enrichKEGG`).
  4. Visualize top enriched pathways using a dot plot.
- **Output**: List of enriched pathways with p-values and gene ratios.

## How to Use the Pipeline
1. **Install Required R Packages**: Make sure to install all listed packages (`DESeq2`, `ggplot2`, `pheatmap`, `FactoMineR`, etc.).
2. **Modify Input File Paths**: Update `file_1`, `file_2`, `file_3`, etc., to match your local paths. Add or remove files as needed, and update the named list.
3. **Adjust Conditions/Groups**:
   - Adapt the `group = ifelse(...)` logic to reflect your study design. For instance, you may have multiple treatments or multiple control groups.
4. **Tweak Statistical Thresholds**:
   - The script uses `padj < 0.05` and `|log2FoldChange| > 1` as default significance thresholds. Adjust these as appropriate for your research.
5. **Run the Script**:
   - Source the script in an R session or run it line by line.
6. **Examine Outputs**:
   - PCA plots, volcano plots, DE results, heatmaps, and KEGG results (if performed) will be saved to your working directory. The final DE results are also saved to an Excel file (`differential_expression_results_DESeq2.xlsx`).

## Customization Tips
1. **Experimental Design**:
   - If you have more complex designs (e.g., `~ batch + condition + group`), update the `design` accordingly.
2. **Partial Analyses**:
   - If you only need the DE analysis for one specific condition, subset your data and skip the rest.
3. **Additional Plots**:
   - Consider adding hierarchical clustering, MDS plots, or other data visualization steps for quality control.
4. **Gene Ontology (GO) Enrichment**:
   - The approach for GO enrichment is similar to KEGG, using `enrichGO` from `clusterProfiler`.
5. **Automate Repetitive Tasks**:
   - You can wrap PCA, volcano plot generation, and heatmap creation in functions to automatically handle multiple conditions.

## Software & Environment
- R version: `>= 4.0`
- Operating System: Linux/Mac/Windows
- Key dependencies:
  - `BiocManager` for installing Bioconductor packages
  - `DESeq2` for differential expression
  - `FactoMineR` and `factoextra` for PCA
  - `pheatmap` for heatmaps
  - `clusterProfiler` and `org.Hs.eg.db` for KEGG analysis
  - `openxlsx` for Excel import/export

## Troubleshooting
1. **Mismatch in Sample Names**:
   - If the metadata sample names do not match the columns of the count matrix, the script will throw an error. Double-check your column names and sample labels.
2. **Insufficient Replicates**:
   - DESeq2 and other statistical methods require replicates. If you have fewer than 2 replicates per group, the analysis may fail or yield unreliable results.
3. **No Significant Genes**:
   - This can occur if the data has low statistical power or if differences between groups are subtle. Consider adjusting thresholds or verifying that the data is correct.
4. **Pathway Enrichment Returned Null**:
   - If there aren’t enough significant genes, or if gene mapping fails, enrichment may return no significant pathways.

## References
- **DESeq2**: [Love et al., Genome Biology, 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
- **clusterProfiler**: [Yu et al., OMICS: A Journal of Integrative Biology, 2012](https://doi.org/10.1089/omi.2011.0118)
- **FactoMineR**: [Lê et al., Journal of Statistical Software, 2008](https://www.jstatsoft.org/index.php/jss/article/view/v025i01)

## Contact
- For questions or troubleshooting, please contact: swroblewski@tulane.edu

