# Load necessary libraries
library(consensusDE)
library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(magrittr)
library(edgeR)
library(RUVSeq)
library(Biobase)
library(DESeq2)
library(pheatmap)
library(limma)
library(patchwork)
library(matrixStats)  


# Load data
summarised_de <- readRDS("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/summarised_de_iso.rds")
sample_table <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/sample_table.tsv")

# Library size
library_sizes <- assay(summarised_de) %>% 
    colSums()

summary_stats <- summary(library_sizes)
capture.output(summary_stats, file = "library_sizes.tsv")

# Initialize log
filter_log <- NULL
log_filter_counts <- function(se, step_label, log_df = NULL) {
  tibble(step = step_label, n_samples = ncol(se), n_transcripts = nrow(se)) %>%
    bind_rows(log_df, .)
}

# STEP 0: Raw data
filter_log <- log_filter_counts(summarised_de, "Initial (raw summarised_de)", filter_log)

# STEP 1: consensusDE filtering
summarised_de_filter <- buildSummarized(
  summarized = summarised_de, 
  sample_table = sample_table,
  filter = TRUE, 
  output_log = "/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/consensusDE_log"
)
filter_log <- log_filter_counts(summarised_de_filter, "After buildSummarized (filter = TRUE)", filter_log)

# STEP 2: Remove samples with missing group
summarised_de_filter <- summarised_de_filter[, !is.na(colData(summarised_de_filter)$group)]
filter_log <- log_filter_counts(summarised_de_filter, "Removed samples with NA group", filter_log)

# STEP 3: Remove all-zero transcripts
summarised_de_filter <- summarised_de_filter[rowSums(assay(summarised_de_filter)) > 0, ]
filter_log <- log_filter_counts(summarised_de_filter, "Removed transcripts with all-zero counts", filter_log)

# Library size - post filtering and pre CPM
library_sizes <- assay(summarised_de_filter) %>% 
    colSums()

summary_stats <- summary(library_sizes)
capture.output(summary_stats, file = "library_sizes_filter_pre_cpm.tsv")

# STEP 4: CPM > X in ≥ X samples - only has an impact if you adjust the sample number
# otherwise, buildSummarized will already clean up the data
keep <- assay(summarised_de_filter) %>%
  DGEList() %>%
  cpm() %>%
  { rowSums(. > 1) >= 10 }
summarised_de_filter <- summarised_de_filter[keep, ]
filter_log <- log_filter_counts(summarised_de_filter, "Filtered by CPM > 1 in ≥3 samples", filter_log)

# Library size - post CPM
library_sizes <- assay(summarised_de_filter) %>% 
    colSums()

summary_stats <- summary(library_sizes)
capture.output(summary_stats, file = "library_sizes_filter_post_cpm.tsv")

# STEP 5: Remove high-count transcripts (mean > 50k)
transcript_means <- rowMeans(assay(summarised_de_filter))
summarised_de_filter <- summarised_de_filter[transcript_means < 5e4, ]
filter_log <- log_filter_counts(summarised_de_filter, "Removed transcripts with mean > 50,000", filter_log)

# Save high-count summary
tibble(
  threshold = c(">50k", ">100k"),
  n_transcripts = c(sum(transcript_means > 5e4), sum(transcript_means > 1e5))
) %>% write_tsv("qc_high_expression_transcripts.tsv")

# STEP 6: Create SeqExpressionSet for RUV
# Extract colData and convert to data frame

# Combine IDC and metadata
idc_df <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/idc_mapping/pk_ref/scaden_mean_proportions.txt") %>% 
  pivot_longer(cols = -sample_id, names_to = "celltype", values_to = "proportion") %>% 
  group_by(sample_id) %>%
  slice_max(order_by = proportion, n = 1) %>%
  ungroup() %>% 
  select(-proportion) %>%
  left_join(
      read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/idc_mapping/pk_ref/scaden_mean_proportions.txt")
  ) %>%
  mutate(sample_id = str_replace(sample_id, "-", "_")) %>%
  mutate(schizont = schizont * 100,
    ring = ring * 100,
    trophozoite = trophozoite * 100)

meta_idc_combined <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/full_meta_for_de.tsv") %>%
  dplyr::rename(sample = file) %>%
  mutate(sample_id = str_remove(sample, ".txt")) %>%
  left_join(
    idc_df,
    by = "sample_id"
  )
    
# Combine with de object
coldata <- colData(summarised_de_filter) %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(
    meta_idc_combined, 
    by = "sample") %>%
  column_to_rownames("sample") %>%
  DataFrame() 

colData(summarised_de_filter) <- coldata

pheno <- AnnotatedDataFrame(as(colData(summarised_de_filter), "data.frame"))
ruv_set <- newSeqExpressionSet(
  counts = as.matrix(assay(summarised_de_filter)),
  phenoData = pheno
)

# STEP 7: Compute residuals from a first-pass GLM regression
design <- model.matrix(~batch + group, data = pData(ruv_set))
dge <- DGEList(counts = counts(ruv_set), group = pData(ruv_set)$group)
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
residuals_matrix <- residuals(fit, type = "deviance")

# STEP 8: Apply RUVr to estimate unwanted variation factors
control_genes <- rownames(ruv_set)  # Using all genes as controls
k <- 1  # Number of unwanted factors to estimate; adjust based on your data
ruv_set <- RUVr(ruv_set, cIdx = control_genes, k = k, residuals = residuals_matrix)

# STEP 9: Incorporate RUV factors into colData
colData(summarised_de_filter)$W_1 <- pData(ruv_set)$W_1

# STEP 10: Remove rows with any NA in counts
summarised_de_filter <- summarised_de_filter[rowSums(is.na(assay(summarised_de_filter))) == 0, ]
filter_log <- log_filter_counts(summarised_de_filter, "Removed transcripts with any NA values", filter_log)

# STEP 11: Export log
write_tsv(filter_log, "filtering_summary_log.tsv")

# Save as GDS for easy access (step above can take a while)
#saveRDS(summarised_de_filter, file = "summarized_de_filter.rds")
#summarised_de_filter <- readRDS("summarized_de_filter.rds")

#####################################################


# Plot distribution of counts across transcripts

# Convert assay to long format
count_distribution <- assay(summarised_de_filter) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count")

# Optional: remove 0s for cleaner plot (keeps > 0 counts only)
count_distribution_filtered <- count_distribution %>%
  filter(count > 0)

# Plot transcript expression distribution
transcript_means <- assay(summarised_de_filter) %>%
  rowMeans()

# Calculate quantiles
low_cutoff <- quantile(transcript_means, 0.05)
high_cutoff <- quantile(transcript_means, 0.99)

# Plot with overlays
export_plot <- ggplot(data.frame(mean_count = transcript_means), aes(x = mean_count)) +
  geom_histogram(bins = 100, fill = "grey70") +
  geom_vline(xintercept = low_cutoff, color = "#440154FF", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = high_cutoff, color = "#55C667FF", linetype = "dashed", linewidth = 1) +
  scale_x_log10() +
  labs(
    title = "Transcript Expression Distribution",
    x = "Mean Counts (log10)",
    y = "Transcript Count"
  )

ggsave("transcript_expression_distribution_with_cutoffs.png", plot = export_plot, dpi = 300)



#####################################################

# DE analysis
# Build DESeq2 object
# Relevel the group factor in colData
colData(summarised_de_filter)$group <- relevel(colData(summarised_de_filter)$group, ref = "Um")

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(assay(summarised_de_filter)),
  colData = colData(summarised_de_filter),
  design = ~ batch + W_1 + group
)

# Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)

# Now this line reflects "Sm vs Um"
res <- results(dds, name = "group_Sm_vs_Um")


# Save p-value distribution plot
export_plot <- ggplot(as.data.frame(res), aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "grey30") +
  labs(title = "DESeq2 P-value Distribution", x = "p-value", y = "Frequency")

ggsave("pval_hist_deseq2.png", plot = export_plot, dpi = 300)


# PCA
plot_pca_pairs <- function(pca_scores,
                           color_var,
                           tool_label = "deseq2",
                           filename_prefix = "pca",
                           log_transform = FALSE,
                           color_scheme = NULL,
                           pcs = c(1, 2, 3)) {

  # If log transformation is needed
  if (log_transform) {
    color_col <- paste0("log10_", color_var)
    pca_scores[[color_col]] <- log10(pca_scores[[color_var]])
    color_var <- color_col
  }

  # Loop over all PC combinations (e.g., PC1 vs PC2, PC1 vs PC3)
  for (i in seq_along(pcs)) {
    for (j in seq_along(pcs)) {
      if (i < j) {
        pc_x <- pcs[i]
        pc_y <- pcs[j]
        pcx_name <- paste0("PC", pc_x)
        pcy_name <- paste0("PC", pc_y)

        x_lab <- paste0("PC", pc_x, ": ", round(attr(pca_scores, "var_expl")[pc_x] * 100, 2), "% variance")
        y_lab <- paste0("PC", pc_y, ": ", round(attr(pca_scores, "var_expl")[pc_y] * 100, 2), "% variance")

        export_plot <- ggplot(pca_scores, aes_string(x = pcx_name, y = pcy_name)) +
          geom_point(aes_string(colour = color_var)) +
          xlab(x_lab) +
          ylab(y_lab)

        # Add color scale
        if (is.numeric(pca_scores[[color_var]])) {
          export_plot <- export_plot + scale_color_viridis_c()
        } else if (!is.null(color_scheme)) {
          export_plot <- export_plot + scale_color_manual(values = color_scheme)
        }

        # Save the file
        filename <- sprintf("%s_%s_%s_PC%s_PC%s.png", filename_prefix, tool_label, gsub("log10_", "", color_var), pc_x, pc_y)
        ggsave(filename, plot = export_plot, dpi = 300)
      }
    }
  }
}

# Create PCA object and attach variance info
vsd <- vst(dds, blind = TRUE)
pca_obj <- prcomp(t(assay(vsd)))
var_expl <- summary(pca_obj)$importance[2, ]

pca_scores <- data.frame(
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    PC3 = pca_obj$x[, 3],
    sample = rownames(pca_obj$x),
    group = colData(vsd)$group
    ) %>%
    mutate(sample = str_remove(sample, "_DK.*|_HALF.*|.txt"),
        sample = str_remove(sample, "ZS_")) %>%
        left_join(
          meta_idc_combined %>%
            mutate(sample = sampleid)
        )

attr(pca_scores, "var_expl") <- var_expl

plot_pca_pairs(pca_scores, color_var = "group", tool_label = "deseq2", color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "cluster", tool_label = "deseq2", color_scheme = c("Mf" = "#440154FF", "Mn" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "parasitemia", tool_label = "deseq2", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "age", tool_label = "deseq2", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "batch", tool_label = "deseq2", color_scheme = c("batch_1" = "#440154FF", "batch_2" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "sex", tool_label = "deseq2", color_scheme = c("F" = "#440154FF", "M" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "deseq2", color_scheme = c("trophozoite" = "#440154FF", "schizont" = "#55C667FF", "ring" = "#FDE725"))
plot_pca_pairs(pca_scores, color_var = "schizont", tool_label = "deseq2", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "ring", tool_label = "deseq2", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "trophozoite", tool_label = "deseq2", log_transform = TRUE)


# Pf IDC

idc_df <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/idc_mapping/pf_ref/scaden_mean_proportions.txt") %>% 
  pivot_longer(cols = -sample_id, names_to = "celltype", values_to = "proportion") %>% 
  group_by(sample_id) %>%
  slice_max(order_by = proportion, n = 1) %>%
  ungroup() %>% 
  select(-proportion) %>%
  mutate(sample_id = str_replace(sample_id, "-", "_")) 

pca_scores <- pca_scores %>% 
  select(-c(celltype, schizont, ring, trophozoite)) %>%
  left_join(
    idc_df
  )

plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "deseq2_pf", 
  color_scheme = c("trophozoite" = "#440154FF", 
                   "schizont" = "#55C667FF", 
                   "ring" = "black", 
                   "gametocyte_female" = "#3B528BFF", 
                     "gametocyte_developing" = "#21908CFF", 
                     "gametocyte_male" = "#FDE725FF"))


# Volcano plot
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(sig = padj < 0.001 & abs(log2FoldChange) > 1)

export_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey70", "#f57d15")) +
  labs(title = "Volcano Plot (DESeq2)", x = "log2 Fold Change", y = "-log10 p-value")

ggsave("volcano_deseq2.png", plot = export_plot, dpi = 300)

# Heatmap of top 30 DE genes
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(30) %>%
  pull(gene)

export_plot <- pheatmap(assay(vsd)[top_genes, ],
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        annotation_col = as.data.frame(colData(vsd)[, "group", drop = FALSE]),
                        show_rownames = TRUE,
                        show_colnames = FALSE)

# Save heatmap as file
ggsave("heatmap_deseq2_top30.png", plot = export_plot[[4]], dpi = 300)

res_deseq2_full <- as.data.frame(res) %>%
  rownames_to_column("gene")

write_tsv(res_deseq2_full, "deseq2_results_full.tsv")



# PCA on top 1,000 most variable genes, NOT accounting for batch effect
vsd <- vst(dds, blind = TRUE)
top_variable_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 1000)
vsd_top_var <- vsd[top_variable_genes, ]
pca_results <- prcomp(t(assay(vsd_top_var)))
var_expl <- summary(pca_results)$importance[2, ]


# Prepare PCA scores for plotting
pca_scores <- as.data.frame(pca_results$x) %>%
  rownames_to_column("file") %>%
  left_join(as.data.frame(colData(vsd_top_var)), by = "file")

attr(pca_scores, "var_expl") <- var_expl

plot_pca_pairs(pca_scores, color_var = "group", tool_label = "deseq2_top_1000", color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "cluster", tool_label = "deseq2_top_1000", color_scheme = c("Mf" = "#440154FF", "Mn" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "batch", tool_label = "deseq2_top_1000", color_scheme = c("batch_1" = "#440154FF", "batch_2" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "deseq2_top_1000", color_scheme = c("trophozoite" = "#440154FF", "schizont" = "#55C667FF", "ring" = "#FDE725"))
plot_pca_pairs(pca_scores, color_var = "schizont", tool_label = "deseq2_top_1000", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "ring", tool_label = "deseq2_top_1000", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "trophozoite", tool_label = "deseq2_top_1000", log_transform = TRUE)


# Pf IDC

idc_df <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/idc_mapping/pf_ref/scaden_mean_proportions.txt") %>% 
  pivot_longer(cols = -sample_id, names_to = "celltype", values_to = "proportion") %>% 
  group_by(sample_id) %>%
  slice_max(order_by = proportion, n = 1) %>%
  ungroup() %>% 
  select(-proportion) %>%
  mutate(sample_id = str_replace(sample_id, "-", "_")) 

pca_scores <- pca_scores %>% 
  select(-c(celltype, schizont, ring, trophozoite)) %>%
  left_join(
    idc_df
  )

plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "deseq2_pf_1000", 
  color_scheme = c("trophozoite" = "#440154FF", 
                   "schizont" = "#55C667FF", 
                   "ring" = "black", 
                   "gametocyte_female" = "#3B528BFF", 
                     "gametocyte_developing" = "#21908CFF", 
                     "gametocyte_male" = "#FDE725FF"))

# PCA on top 1,000 most variable genes whilst accounting for batch effect

# Apply variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
assay(vsd) <- removeBatchEffect(assay(vsd), 
  batch = vsd$batch, 
  covariates = vsd$W_1)
top_variable_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 1000)
vsd_top_var <- vsd[top_variable_genes, ]
pca_results <- prcomp(t(assay(vsd_top_var)))
var_expl <- summary(pca_results)$importance[2, ]

# Prepare PCA scores for plotting
pca_scores <- as.data.frame(pca_results$x) %>%
  rownames_to_column("file") %>%
  left_join(as.data.frame(colData(vsd_top_var)), by = "file")

attr(pca_scores, "var_expl") <- var_expl

plot_pca_pairs(pca_scores, color_var = "group", tool_label = "deseq2_top_1000_batch_corrected", color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "cluster", tool_label = "deseq2_top_1000_batch_corrected", color_scheme = c("Mf" = "#440154FF", "Mn" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "batch", tool_label = "deseq2_top_1000_batch_corrected", color_scheme = c("batch_1" = "#440154FF", "batch_2" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "deseq2_top_1000_batch_corrected", color_scheme = c("trophozoite" = "#440154FF", "schizont" = "#55C667FF", "ring" = "#FDE725"))

# Pf IDC

idc_df <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/idc_mapping/pf_ref/scaden_mean_proportions.txt") %>% 
  pivot_longer(cols = -sample_id, names_to = "celltype", values_to = "proportion") %>% 
  group_by(sample_id) %>%
  slice_max(order_by = proportion, n = 1) %>%
  ungroup() %>% 
  select(-proportion) %>%
  mutate(sample_id = str_replace(sample_id, "-", "_")) 

pca_scores <- pca_scores %>% 
  select(-c(celltype, schizont, ring, trophozoite)) %>%
  left_join(
    idc_df
  )

plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "deseq2_pf_1000_batch_corrected", 
  color_scheme = c("trophozoite" = "#440154FF", 
                   "schizont" = "#55C667FF", 
                   "ring" = "black", 
                   "gametocyte_female" = "#3B528BFF", 
                     "gametocyte_developing" = "#21908CFF", 
                     "gametocyte_male" = "#FDE725FF"))

# PCA accounting for batch on all genes

# Perform PCA
pca_results <- prcomp(t(assay(vsd)))

# Prepare PCA scores for plotting
pca_scores <- as.data.frame(pca_results$x) %>%
  rownames_to_column("file") %>%
  left_join(as.data.frame(colData(vsd_top_var)), by = "file")

plot_pca_pairs(pca_scores, color_var = "group", tool_label = "deseq2_batch_corrected_group", color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))



# PCA on significant genes

sig_genes_deseq2 <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(gene)

vsd_sig <- vsd[sig_genes_deseq2, ]
pca_results <- prcomp(t(assay(vsd_top_var)))
var_expl <- summary(pca_results)$importance[2, ]


# Prepare PCA scores for plotting
pca_scores <- as.data.frame(pca_results$x) %>%
  rownames_to_column("file") %>%
  left_join(as.data.frame(colData(vsd_top_var)), by = "file")

attr(pca_scores, "var_expl") <- var_expl

plot_pca_pairs(pca_scores, color_var = "group", tool_label = "deseq2_sig", color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "cluster", tool_label = "deseq2_sig", color_scheme = c("Mf" = "#440154FF", "Mn" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "batch", tool_label = "deseq2_sig", color_scheme = c("batch_1" = "#440154FF", "batch_2" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "deseq2_sig", color_scheme = c("trophozoite" = "#440154FF", "schizont" = "#55C667FF", "ring" = "#FDE725"))


# Pf IDC

idc_df <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/idc_mapping/pf_ref/scaden_mean_proportions.txt") %>% 
  pivot_longer(cols = -sample_id, names_to = "celltype", values_to = "proportion") %>% 
  group_by(sample_id) %>%
  slice_max(order_by = proportion, n = 1) %>%
  ungroup() %>% 
  select(-proportion) %>%
  mutate(sample_id = str_replace(sample_id, "-", "_")) 

pca_scores <- pca_scores %>% 
  select(-c(celltype, schizont, ring, trophozoite)) %>%
  left_join(
    idc_df
  )

plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "deseq2_pf_sig", 
  color_scheme = c("trophozoite" = "#440154FF", 
                   "schizont" = "#55C667FF", 
                   "ring" = "black", 
                   "gametocyte_female" = "#3B528BFF", 
                     "gametocyte_developing" = "#21908CFF", 
                     "gametocyte_male" = "#FDE725FF"))


# Dispersion plots - should follow red trend line
png("deseq2_dispersion_plot.png", width = 1800, height = 1500, res = 300)
plotDispEsts(dds)
dev.off()

# Cooks distance - influential samples
png("deseq2_cooks_distance.png", width = 1800, height = 1500, res = 300)
plotCounts(dds, gene = which.max(assays(dds)[["cooks"]][,1]), intgroup = "group")
dev.off()

# MA Plot
res_df_ma <- as.data.frame(res)
res_df_ma$neg_log_pval <- -log10(res_df_ma$pvalue) # Add a column for -log10(p-value) for coloring

# Create the MA plot
ma_plot_deseq2 <- ggplot(res_df_ma, aes(x = log10(baseMean), y = log2FoldChange, color = neg_log_pval)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_c(option = "viridis") +
  labs(
    x = "Log10 Mean Expression",
    y = "Log2 Fold Change",
    color = "-Log10(P-Value)"
  )

# Save the plot to a file
ggsave(filename = "deseq2_MA_Plot.png", plot = ma_plot_deseq2, width = 10, dpi = 300)



#####################################################



# edgeR

# Setup
counts <- round(assay(summarised_de_filter))
meta <- colData(summarised_de_filter)

# Extract RUV factors (e.g., W_1, W_2 if you used k = 2)
W_vars <- grep("^W_", colnames(meta), value = TRUE)
design_formula <- as.formula(paste("~ batch +", paste(c(W_vars, "group"), collapse = " + ")))
design <- model.matrix(design_formula, data = meta)

# Create DGEList
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# Estimate dispersions
dge <- estimateDisp(dge, design)

# Fit GLM
fit <- glmFit(dge, design)

# Contrast (adjust based on your actual group levels)
# If your levels are group_Control and group_Treated, check colnames(design)
contrast_name <- grep("^group", colnames(design), value = TRUE)
lrt <- glmLRT(fit, coef = contrast_name)

# Extract results
res_edgeR <- topTags(lrt, n = Inf)$table

# Save p-value distribution plot
export_plot <- ggplot(res_edgeR, aes(x = PValue)) +
  geom_histogram(bins = 50, fill = "grey30") +
  labs(title = "edgeR P-value Distribution", x = "p-value", y = "Frequency")

ggsave("pval_hist_edgeR.png", plot = export_plot, dpi = 300)

# PCA from log-CPM
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
pca_obj <- prcomp(t(logcpm))
var_expl <- summary(pca_obj)$importance[2, ]

# Join metadata
pca_scores <- data.frame(
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    PC3 = pca_obj$x[, 3],
    sample = colnames(logcpm),
    group = meta$group,
    W_1 = meta$W_1
    ) %>%
    mutate(sample = str_remove(sample, "_DK.*|_HALF.*|.txt"),
        sample = str_remove(sample, "ZS_")) %>%
        left_join(
          meta_idc_combined %>%
            mutate(sample = sampleid)
        )

attr(pca_scores, "var_expl") <- var_expl

plot_pca_pairs(pca_scores, color_var = "group", tool_label = "edger", color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "cluster", tool_label = "edger", color_scheme = c("Mf" = "#440154FF", "Mn" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "parasitemia", tool_label = "edger", color_scheme = "viridis_continuous")
plot_pca_pairs(pca_scores, color_var = "age", tool_label = "edger", color_scheme = "viridis_continuous")
plot_pca_pairs(pca_scores, color_var = "batch", tool_label = "edger", color_scheme = c("batch_1" = "#440154FF", "batch_2" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "sex", tool_label = "edger", color_scheme = c("F" = "#440154FF", "M" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "edger", color_scheme = c("trophozoite" = "#440154FF", "schizont" = "#55C667FF", "ring" = "#FDE725"))
plot_pca_pairs(pca_scores, color_var = "schizont", tool_label = "edger", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "ring", tool_label = "edger", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "trophozoite", tool_label = "edger", log_transform = TRUE)

# Volcano plot
res_edgeR_vol <- res_edgeR %>%
  rownames_to_column("gene") %>%
  mutate(sig = FDR < 0.05 & abs(logFC) > 1)

export_plot <- ggplot(res_edgeR_vol, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey70", "#f57d15")) +
  labs(title = "Volcano Plot (edgeR)", x = "log2 Fold Change", y = "-log10 p-value")

ggsave("volcano_edgeR.png", plot = export_plot, dpi = 300)

# Heatmap of top 30 DE genes
top_genes <- res_edgeR_vol %>%
  arrange(FDR) %>%
  head(30) %>%
  pull(gene)

# Set up annotation with proper rownames
annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- colnames(logcpm)

# Generate heatmap and assign to export_plot
export_plot <- pheatmap(logcpm[top_genes, ],
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        annotation_col = annotation_col,
                        show_rownames = TRUE,
                        show_colnames = FALSE)

ggsave("heatmap_edgeR_top30.png", plot = export_plot[[4]], dpi = 300)

res_edgeR_full <- topTags(lrt, n = Inf)$table %>%
  rownames_to_column("gene")

write_tsv(res_edgeR_full, "edgeR_results_full.tsv")

# PCA on significant genes

# Subset significant genes
sig_genes_edger <- res_edgeR %>%
  rownames_to_column("gene") %>%
  filter(!is.na(FDR), FDR < 0.001, abs(logFC) > 1) %>%
  pull(gene)

logcpm_sig <- logcpm[sig_genes_edger, ]

# PCA
pca_data <- prcomp(t(logcpm_sig))
pca_scores <- data.frame(PC1 = pca_data$x[, 1],
                         PC2 = pca_data$x[, 2],
                         sample = colnames(logcpm_sig),
                         group = meta$group)

export_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group)) +
  labs(
    title = "PCA on Significant DE Genes (edgeR)",
    x = paste0("PC1: ", round(summary(pca_data)$importance[2, 1] * 100, 2), "% variance"),
    y = paste0("PC2: ", round(summary(pca_data)$importance[2, 2] * 100, 2), "% variance")
  ) +
  scale_color_manual(values = c("Um" = "#440154FF", "Sm" = "#55C667FF"))

ggsave("pca_edger_sig_genes.png", plot = export_plot, dpi = 300)


# PCA on top 1,000 most variable genes

# Calculate variance for each gene
gene_variances <- rowVars(logcpm)

# Select the top 1,000 most variable genes
top_variable_genes <- order(gene_variances, decreasing = TRUE)[1:1000]

# Subset the log-CPM matrix to include only these genes
logcpm_top_var <- logcpm[top_variable_genes, ]

# Perform PCA
pca_results <- prcomp(t(logcpm_top_var))

pca_scores <- data.frame(
  PC1 = pca_results$x[, 1],
  PC2 = pca_results$x[, 2],
  sample = colnames(logcpm_top_var),
    group = meta$group,
    W_1 = meta$W_1
)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(
    title = "PCA on Top 1,000 Most Variable Genes (edgeR)",
    x = paste0("PC1: ", round(summary(pca_results)$importance[2, 1] * 100, 2), "% variance"),
    y = paste0("PC2: ", round(summary(pca_results)$importance[2, 2] * 100, 2), "% variance")
  ) +
  scale_color_manual(values = c("Um" = "#440154FF", "Sm" = "#55C667FF"))

ggsave("pca_edger_top_1000_variable_genes.png", plot = pca_plot, dpi = 300)

# Biological coefficient of variation (BCV) plot - “clean” curve suggests stable dispersion estimate
png("edger_bcv_plot.png", width = 1800, height = 1500, res = 300)
plotBCV(dge)
dev.off()


# MA Plot
res_df_ma <- as.data.frame(res_edgeR)
res_df_ma$neg_log_pval <- -log10(res_df_ma$PValue)

ma_plot_edger <- ggplot(res_df_ma, aes(x = logCPM, y = logFC, color = neg_log_pval)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_c(option = "viridis") +
  labs(
    x = "Log CPM",
    y = "Log Fold Change",
    color = "-Log10(P-Value)"
  )

ggsave(filename = "edger_MA_Plot.png", plot = ma_plot_edger, width = 10, dpi = 300)



#####################################################



# voom-limma

# Setup
counts <- round(assay(summarised_de_filter))
meta <- colData(summarised_de_filter)

# Extract RUV factors (e.g., W_1, W_2 if k = 2)
W_vars <- grep("^W_", colnames(meta), value = TRUE)

# Ensure 'group' is a factor
meta$group <- factor(meta$group)

# Design matrix including RUV factors and group
design_formula <- as.formula(paste("~ batch +", paste(c(W_vars, "group"), collapse = " + ")))
design <- model.matrix(design_formula, data = meta)

# Create DGEList object
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# Apply voom transformation
v <- voom(dge, design, plot = TRUE)

# Save the mean-variance trend plot
ggsave("voom_mean_variance_plot.png", dpi = 300)

# Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract results
voom_res <- topTable(fit, coef = "groupSm", number = Inf, sort.by = "P")

# Save p-value histogram
export_plot <- ggplot(voom_res, aes(x = P.Value)) +
  geom_histogram(bins = 50, fill = "grey30") +
  labs(title = "Voom P-value Distribution", x = "p-value", y = "Frequency")

ggsave("pval_hist_voom.png", plot = export_plot, dpi = 300)

# PCA from voom logCPM
pca_obj <- prcomp(t(v$E))
var_expl <- summary(pca_obj)$importance[2, ]

# Create PCA scores dataframe
pca_scores <- data.frame(
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    PC3 = pca_obj$x[, 3],
    sample = colnames(v$E),
    group = meta$group,
    W_1 = meta$W_1
    ) %>%
    mutate(sample = str_remove(sample, "_DK.*|_HALF.*|.txt"),
        sample = str_remove(sample, "ZS_")) %>%
        left_join(
          meta_idc_combined %>%
            mutate(sample = sampleid)
        )


attr(pca_scores, "var_expl") <- var_expl

plot_pca_pairs(pca_scores, color_var = "group", tool_label = "voom", color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "cluster", tool_label = "voom", color_scheme = c("Mf" = "#440154FF", "Mn" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "parasitemia", tool_label = "voom", color_scheme = "viridis_continuous")
plot_pca_pairs(pca_scores, color_var = "age", tool_label = "voom", color_scheme = "viridis_continuous")
plot_pca_pairs(pca_scores, color_var = "batch", tool_label = "voom", color_scheme = c("batch_1" = "#440154FF", "batch_2" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "sex", tool_label = "voom", color_scheme = c("F" = "#440154FF", "M" = "#55C667FF"))
plot_pca_pairs(pca_scores, color_var = "celltype", tool_label = "voom", color_scheme = c("trophozoite" = "#440154FF", "schizont" = "#55C667FF", "ring" = "#FDE725"))
plot_pca_pairs(pca_scores, color_var = "schizont", tool_label = "voom", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "ring", tool_label = "voom", log_transform = TRUE)
plot_pca_pairs(pca_scores, color_var = "trophozoite", tool_label = "voom", log_transform = TRUE)

# Volcano plot
voom_res_df <- voom_res %>%
  rownames_to_column("gene") %>%
  mutate(sig = adj.P.Val < 0.05 & abs(logFC) > 1)

export_plot <- ggplot(voom_res_df, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey70", "#f57d15")) +
  labs(title = "Volcano Plot (voom)", x = "log2 Fold Change", y = "-log10 p-value")

ggsave("volcano_voom.png", plot = export_plot, dpi = 300)

# Heatmap of top 30 DE genes
top_genes_voom <- voom_res_df %>%
  arrange(adj.P.Val) %>%
  head(30) %>%
  pull(gene)

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- colnames(v$E)

export_plot <- pheatmap(v$E[top_genes_voom, ],
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        annotation_col = annotation_col,
                        show_rownames = TRUE,
                        show_colnames = FALSE)

png("heatmap_voom_top30.png", width = 2000, height = 2000, res = 300)
grid::grid.newpage()
grid::grid.draw(export_plot$gtable)
dev.off()

voom_res_full <- voom_res %>%
  rownames_to_column("gene")

write_tsv(voom_res_full, "voom_results_full.tsv")

# PCA on significant genes

# Subset significant genes
sig_genes_voom <- voom_res_df %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) > 1) %>%
  pull(gene)

voom_sig <- v$E[sig_genes_voom, ]

# PCA
pca_data <- prcomp(t(voom_sig))
pca_scores <- data.frame(PC1 = pca_data$x[, 1],
                         PC2 = pca_data$x[, 2],
                         sample = colnames(voom_sig),
                         group = meta$group)

export_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group)) +
  labs(
    title = "PCA on Significant DE Genes (voom)",
    x = paste0("PC1: ", round(summary(pca_data)$importance[2, 1] * 100, 2), "% variance"),
    y = paste0("PC2: ", round(summary(pca_data)$importance[2, 2] * 100, 2), "% variance")
  ) +
  scale_color_manual(values = c("Um" = "#440154FF", "Sm" = "#55C667FF"))

ggsave("pca_voom_sig_genes.png", plot = export_plot, dpi = 300)


# PCA on top 1,000 most variable genes

# Calculate variance for each gene
gene_variances <- rowVars(v$E)

# Select the top 1,000 most variable genes
top_variable_genes <- order(gene_variances, decreasing = TRUE)[1:1000]

# Subset the expression matrix to include only these genes
voom_top_var <- v$E[top_variable_genes, ]

# Perform PCA
pca_results <- prcomp(t(voom_top_var))

pca_scores <- data.frame(
  PC1 = pca_results$x[, 1],
  PC2 = pca_results$x[, 2],
  sample = colnames(voom_top_var),
  group = meta$group
)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(
    title = "PCA on Top 1,000 Most Variable Genes (voom)",
    x = paste0("PC1: ", round(summary(pca_results)$importance[2, 1] * 100, 2), "% variance"),
    y = paste0("PC2: ", round(summary(pca_results)$importance[2, 2] * 100, 2), "% variance")
  ) +
  scale_color_manual(values = c("Um" = "#440154FF", "Sm" = "#55C667FF"))

ggsave("pca_voom_top_1000_variable_genes.png", plot = pca_plot, dpi = 300)

# Residual Standard Deviation Plot
design_simple <- model.matrix(~ meta$group)
v_simple <- voom(dge, design_simple, plot = TRUE)
ggsave("voom_mean_variance_plot.png", dpi = 300)

saveRDS(v$weights, file = "voom_weights.rds")

png("voom_mean_variance_trend.png", width = 1800, height = 1500, res = 300)
v <- voom(dge, design, plot = TRUE)
dev.off()

# MA Plot
voom_res$neg_log_pval <- -log10(voom_res$P.Value)
ma_plot_voom <- ggplot(voom_res, aes(x = AveExpr, y = logFC, color = neg_log_pval)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_c(option = "viridis") +
  labs(
    x = "Average Expression",
    y = "Log Fold Change",
    color = "-Log10(P-Value)"
  )

# Save the plot to a file
ggsave(filename = "voom_MA_Plot.png", plot = ma_plot_voom, width = 10, dpi = 300)


# Batch-corrected PCA on top 1,000 most variable genes (voom)

# Apply removeBatchEffect to voom logCPM
voom_corrected <- removeBatchEffect(v$E, batch = meta$batch, covariates = meta$W_1)

# Identify top 1,000 most variable genes after correction
top_variable_genes_voom <- order(rowVars(voom_corrected), decreasing = TRUE)[1:1000]
voom_top_var_corrected <- voom_corrected[top_variable_genes_voom, ]

# Perform PCA
pca_results_voom <- prcomp(t(voom_top_var_corrected))
var_expl_voom <- summary(pca_results_voom)$importance[2, ]

# Prepare PCA scores
meta_df <- as.data.frame(meta) %>%
  rownames_to_column("sample")

pca_scores_voom <- data.frame(pca_results_voom$x) %>%
  rownames_to_column("sample") %>%
  left_join(meta_df, by = "sample")

# Add IDC classification
pca_scores_voom <- pca_scores_voom %>%
  select(-c(celltype, schizont, ring, trophozoite)) %>%
  left_join(idc_df, by = c("sample" = "sample_id"))

attr(pca_scores_voom, "var_expl") <- var_expl_voom

# Plot
plot_pca_pairs(pca_scores_voom, color_var = "group", tool_label = "voom_top_1000_batch_corrected", color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))
plot_pca_pairs(pca_scores_voom, color_var = "cluster", tool_label = "voom_top_1000_batch_corrected", color_scheme = c("Mf" = "#440154FF", "Mn" = "#55C667FF"))
plot_pca_pairs(pca_scores_voom, color_var = "batch", tool_label = "voom_top_1000_batch_corrected", color_scheme = c("batch_1" = "#440154FF", "batch_2" = "#55C667FF"))
plot_pca_pairs(pca_scores_voom, color_var = "celltype", tool_label = "voom_top_1000_batch_corrected", color_scheme = c("trophozoite" = "#440154FF", "schizont" = "#55C667FF", "ring" = "#FDE725"))


# Pf IDC for voom

idc_df <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/idc_mapping/pf_ref/scaden_mean_proportions.txt") %>% 
  pivot_longer(cols = -sample_id, names_to = "celltype", values_to = "proportion") %>% 
  group_by(sample_id) %>%
  slice_max(order_by = proportion, n = 1) %>%
  ungroup() %>% 
  select(-proportion) %>%
  mutate(sample_id = str_replace(sample_id, "-", "_")) 

pca_scores_voom <- pca_scores_voom %>% 
  select(-any_of(c("celltype", "schizont", "ring", "trophozoite"))) %>%
  left_join(idc_df, by = c("sample" = "sample_id"))

plot_pca_pairs(pca_scores_voom, color_var = "celltype", tool_label = "voom_pf_1000_batch_corrected", 
  color_scheme = c("trophozoite" = "#440154FF", 
                   "schizont" = "#55C667FF", 
                   "ring" = "black", 
                   "gametocyte_female" = "#3B528BFF", 
                   "gametocyte_developing" = "#21908CFF", 
                   "gametocyte_male" = "#FDE725FF"))

# PCA on all genes with batch correction (voom)

# Apply removeBatchEffect across all genes
voom_corrected_all <- removeBatchEffect(v$E, batch = meta$batch, covariates = meta$W_1)

# PCA
pca_results_all <- prcomp(t(voom_corrected_all))
var_expl_all <- summary(pca_results_all)$importance[2, ]

# Prepare PCA scores
pca_scores_all <- data.frame(pca_results_all$x) %>%
  rownames_to_column("sample") %>%
  left_join(meta_df, by = "sample")


attr(pca_scores_all, "var_expl") <- var_expl_all

# Plot
plot_pca_pairs(pca_scores_all, color_var = "group", tool_label = "voom_batch_corrected_group", 
  color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))


# Batch correction on significant genes (voom)
# Subset significant genes
sig_genes_voom <- voom_res_df %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) > 1) %>%
  pull(gene)

# Expression matrix of significant genes
voom_sig <- v$E[sig_genes_voom, ]

# Apply batch correction
voom_sig_corrected <- removeBatchEffect(voom_sig,
  batch = meta$batch,
  covariates = meta$W_1
)

# PCA
pca_results_sig_corrected <- prcomp(t(voom_sig_corrected))
var_expl <- summary(pca_results_sig_corrected)$importance[2, ]

# Prepare scores
pca_scores_sig_corrected <- data.frame(pca_results_sig_corrected$x) %>%
  rownames_to_column("sample") %>%
  left_join(meta_df, by = "sample")

attr(pca_scores_sig_corrected, "var_expl") <- var_expl

# Plot
plot_pca_pairs(pca_scores_sig_corrected,
               color_var = "group",
               tool_label = "voom_sig_batch_corrected",
               color_scheme = c("Um" = "#440154FF", "Sm" = "#55C667FF"))


#####################################################



# Discordance/concordance

# DESeq2 significant genes
sig_deseq2 <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  select(gene, padj, log2FoldChange) %>%
  dplyr::rename(padj_deseq2 = padj, logFC_deseq2 = log2FoldChange)

# edgeR significant genes
sig_edger <- res_edgeR %>%
  rownames_to_column("gene") %>%
  filter(!is.na(FDR), FDR < 0.05, abs(logFC) > 1) %>%
  select(gene, FDR, logFC) %>%
  dplyr::rename(padj_edger = FDR, logFC_edger = logFC)

# voom significant genes
sig_voom <- voom_res_df %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) > 1) %>%
  select(gene, adj.P.Val, logFC) %>%
  dplyr::rename(padj_voom = adj.P.Val, logFC_voom = logFC)

# Save individually
write_tsv(sig_deseq2, "significant_genes_deseq2.tsv")
write_tsv(sig_edger,   "significant_genes_edger.tsv")
write_tsv(sig_voom,    "significant_genes_voom.tsv")

# Overlap (3-way)
overlap_all <- purrr::reduce(list(sig_deseq2, sig_edger, sig_voom), inner_join, by = "gene")
write_tsv(overlap_all, "significant_overlap_all_methods.tsv")

# Pairwise overlaps
overlap_deseq2_edger <- inner_join(sig_deseq2, sig_edger, by = "gene")
overlap_deseq2_voom  <- inner_join(sig_deseq2, sig_voom,  by = "gene")
overlap_edger_voom   <- inner_join(sig_edger,  sig_voom,  by = "gene")

# Summary
overlap_summary <- tibble(
  Sig_DESeq2  = nrow(sig_deseq2),
  Sig_edgeR   = nrow(sig_edger),
  Sig_voom    = nrow(sig_voom),
  Shared_all  = nrow(overlap_all),
  DESeq2_only = nrow(sig_deseq2) - nrow(bind_rows(overlap_deseq2_edger, overlap_deseq2_voom) %>% distinct(gene)),
  edgeR_only  = nrow(sig_edger)  - nrow(bind_rows(overlap_deseq2_edger, overlap_edger_voom) %>% distinct(gene)),
  voom_only   = nrow(sig_voom)   - nrow(bind_rows(overlap_deseq2_voom, overlap_edger_voom) %>% distinct(gene))
)

write_tsv(overlap_summary, "significant_overlap_summary.tsv")
print(overlap_summary)

# LFC vs LFC
# Extract and align LFCs
lfc_deseq2 <- res_df %>%
  select(log2FoldChange) %>%
  rownames_to_column("gene") %>%
  dplyr::rename(LFC_deseq2 = log2FoldChange)

lfc_edger <- res_edgeR %>%
  select(logFC) %>%
  rownames_to_column("gene") %>%
  dplyr::rename(LFC_edgeR = logFC)

lfc_voom <- voom_res_df %>%
  select(logFC) %>%
  rownames_to_column("gene") %>%
  dplyr::rename(LFC_voom = logFC)

# Merge all into one dataframe
lfc_merged <- purrr::reduce(list(lfc_deseq2, lfc_edger, lfc_voom), full_join, by = "gene") %>%
  drop_na()

# Create each plot
plot1 <- ggplot(lfc_merged, aes(x = LFC_deseq2, y = LFC_edgeR)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "DESeq2 vs edgeR", x = "DESeq2 LFC", y = "edgeR LFC")

plot2 <- ggplot(lfc_merged, aes(x = LFC_deseq2, y = LFC_voom)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "DESeq2 vs voom", x = "DESeq2 LFC", y = "voom LFC")

plot3 <- ggplot(lfc_merged, aes(x = LFC_edgeR, y = LFC_voom)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "edgeR vs voom", x = "edgeR LFC", y = "voom LFC")

# Combine plots
combined_plot <- (plot1 | plot2) / plot3 +
  plot_annotation(title = "Log2 Fold Change Concordance Between Methods")

# Save as PNG
ggsave("lfc_concordance_all_methods.png", plot = combined_plot, width = 12, height = 10, dpi = 300)



####################################################

# Post analysis filtering


# Extract the count matrix and sample metadata
counts_matrix <- assay(summarised_de_filter)
sample_metadata <- colData(summarised_de_filter) %>% as.data.frame()

# Check if 'group' column exists in sample_metadata
if (!"group" %in% colnames(sample_metadata)) {
  stop("The 'group' column is not present in the sample metadata.")
}

# Initialize vectors to store counts and proportions
# overlap_all <- read_tsv("significant_overlap_all_methods.tsv")
num_nonzero_Sm <- integer(nrow(overlap_all))
num_nonzero_Um <- integer(nrow(overlap_all))
prop_nonzero_Sm <- numeric(nrow(overlap_all))
prop_nonzero_Um <- numeric(nrow(overlap_all))

# Total number of samples in each group
total_Sm <- sum(sample_metadata$group == "Sm")
total_Um <- sum(sample_metadata$group == "Um")

# Loop through each gene to calculate counts and proportions
for (i in seq_len(nrow(overlap_all))) {
  gene <- overlap_all$gene[i]
  
  if (gene %in% rownames(counts_matrix)) {
    gene_counts <- counts_matrix[gene, ]
    non_zero_samples <- gene_counts > 10
    
    # Count non-zero samples in each group
    num_nonzero_Sm[i] <- sum(non_zero_samples & sample_metadata$group == "Sm")
    num_nonzero_Um[i] <- sum(non_zero_samples & sample_metadata$group == "Um")
    
    # Calculate proportions
    prop_nonzero_Sm[i] <- num_nonzero_Sm[i] / total_Sm
    prop_nonzero_Um[i] <- num_nonzero_Um[i] / total_Um
  } else {
    # Assign NA if gene is not found in the counts matrix
    num_nonzero_Sm[i] <- NA
    num_nonzero_Um[i] <- NA
    prop_nonzero_Sm[i] <- NA
    prop_nonzero_Um[i] <- NA
  }
}

# Add these counts to the overlap_all data frame
overlap_all_sample_counts <- overlap_all %>%
  mutate(num_nonzero_Sm = num_nonzero_Sm,
         num_nonzero_Um = num_nonzero_Um,
         prop_nonzero_Sm = prop_nonzero_Sm,
         prop_nonzero_Um = prop_nonzero_Um) %>%
         arrange(desc(prop_nonzero_Sm), desc(prop_nonzero_Um))

# Save the updated data frame to a TSV file
write_tsv(overlap_all_sample_counts, "significant_overlap_with_nonzero_counts.tsv")

# Filter the genes based on >= 30% in both groups
overlap_all_sample_counts %>% 
  mutate(parent = str_remove(gene, "_i.*")) %>%
  relocate(parent, .after = gene) %>%
  filter(prop_nonzero_Um >= 0.3 & prop_nonzero_Sm >= 0.3) %>%
  write_tsv("significant_overlap_with_nonzero_counts_30_cutoff.tsv")

# Filter the genes based on >= 50% in both groups
overlap_all_sample_counts %>% 
  mutate(parent = str_remove(gene, "_i.*")) %>%
  relocate(parent, .after = gene) %>%
  filter(prop_nonzero_Um >= 0.5 & prop_nonzero_Sm >= 0.5) %>%
  write_tsv("significant_overlap_with_nonzero_counts_50_cutoff.tsv")

# Filter the genes based on >= 70% in both groups
overlap_all_sample_counts %>% 
  mutate(parent = str_remove(gene, "_i.*")) %>%
  relocate(parent, .after = gene) %>%
  filter(prop_nonzero_Um >= 0.7 & prop_nonzero_Sm >= 0.7) %>%
  write_tsv("significant_overlap_with_nonzero_counts_70_cutoff.tsv")

# Isoforms that come up > 1
overlap_all_sample_counts %>% 
  mutate(parent = str_remove(gene, "_i.*")) %>%
  relocate(parent, .after = gene) %>%
  filter(prop_nonzero_Um >= 0.3 & prop_nonzero_Sm >= 0.3) %>%
  group_by(parent) %>%
  filter(n() > 1) %>%  # keep only parents that appear in >1 row
  ungroup() %>%
  write_tsv("significant_overlap_with_nonzero_counts_30_cutoff_multiple_iso.tsv")








#########################################################