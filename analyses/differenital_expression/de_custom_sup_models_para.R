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
summarised_de <- readRDS("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/summarised_de.rds")
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

# Add log1p-transformed parasitemia BEFORE RUV
colData(summarised_de_filter)$parasitemia_log <- log1p(colData(summarised_de_filter)$parasitemia)

pheno <- AnnotatedDataFrame(as(colData(summarised_de_filter), "data.frame"))
ruv_set <- newSeqExpressionSet(
  counts = as.matrix(assay(summarised_de_filter)),
  phenoData = pheno
)

# STEP 7: Compute residuals from a first-pass GLM regression
design <- model.matrix(~batch + parasitemia_log, data = pData(ruv_set))
dge <- DGEList(counts = counts(ruv_set), group = pData(ruv_set)$parasitemia_log)
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

#colData(summarised_de_filter)$parasitemia_log <- log1p(colData(summarised_de_filter)$parasitemia)


#####################################################



# DE analysis
# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(assay(summarised_de_filter)),
  colData = colData(summarised_de_filter),
  design = ~ batch + W_1 + parasitemia_log
)

# Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name = "parasitemia_log")



#####################################################



# edgeR

# Setup
counts <- round(assay(summarised_de_filter))
meta <- colData(summarised_de_filter)

# Extract RUV factors (e.g., W_1, W_2 if you used k = 2)
W_vars <- grep("^W_", colnames(meta), value = TRUE)
design_formula <- as.formula(paste("~ batch +", paste(c(W_vars, "parasitemia_log"), collapse = " + ")))
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
contrast_name <- grep("^parasitemia_log", colnames(design), value = TRUE)
lrt <- glmLRT(fit, coef = contrast_name)

# Extract results
res_edgeR <- topTags(lrt, n = Inf)$table


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
design_formula <- as.formula(paste("~ batch +", paste(c(W_vars, "parasitemia_log"), collapse = " + ")))
design <- model.matrix(design_formula, data = meta)

# Create DGEList object
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# Apply voom transformation
v <- voom(dge, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract results
voom_res <- topTable(fit, coef = "parasitemia_log", number = Inf, sort.by = "P")

#####################################################

# Discordance/concordance

# Estimate IQR range for log-parasitemia
iqr_vals <- quantile(colData(summarised_de_filter)$parasitemia_log, probs = c(0.25, 0.75))
delta_log <- diff(iqr_vals)  # 75th - 25th percentile


# DESeq2 significant genes
sig_deseq2 <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(delta_IQR_LFC = log2FoldChange * delta_log) %>%
  filter(!is.na(padj), padj < 0.05, abs(delta_IQR_LFC) > 0.75) %>%
  select(gene, padj, log2FoldChange) %>%
  dplyr::rename(padj_deseq2 = padj, logFC_deseq2 = log2FoldChange)

# edgeR significant genes
sig_edger <- res_edgeR %>%
  rownames_to_column("gene") %>%
  mutate(delta_IQR_LFC = logFC * delta_log) %>%
  filter(!is.na(FDR), FDR < 0.05, abs(delta_IQR_LFC) > 0.75) %>%
  select(gene, FDR, logFC) %>%
  dplyr::rename(padj_edger = FDR, logFC_edger = logFC)

# voom significant genes
sig_voom <- voom_res %>%
  rownames_to_column("gene") %>% 
  mutate(delta_IQR_LFC = logFC * delta_log) %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(delta_IQR_LFC) > 0.75) %>%
  select(gene, adj.P.Val, logFC) %>%
  dplyr::rename(padj_voom = adj.P.Val, logFC_voom = logFC)

# Overlap (3-way)
overlap_all <- purrr::reduce(list(sig_deseq2, sig_edger, sig_voom), inner_join, by = "gene")
write_tsv(overlap_all, "significant_overlap_all_methods.tsv")


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
para_de <- overlap_all_sample_counts %>% 
  mutate(parent = str_remove(gene, "_i.*")) %>%
  relocate(parent, .after = gene) %>%
  filter(prop_nonzero_Um >= 0.3 & prop_nonzero_Sm >= 0.3)

para_de %>%
  write_tsv("significant_overlap_with_nonzero_counts_30_cutoff_lfc_1.tsv")


# Compare to severity model

severity_de <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/de_covar/um_control/custom_1cpm_10_samples_ruv/significant_overlap_with_nonzero_counts_30_cutoff.tsv")


overlap_full_transcript <- inner_join(para_de, severity_de, by = "gene")

overlap_full_transcript %>%
  write_tsv("parasitemia_severity_overlap_lfc_1.tsv")

# Need to annotate using workflow in consensusDE.Rmd - generates significant_genes_description_combined_brief.tsv
overlap_full_transcript %>% 
  select(gene) %>%
  left_join(
    read_tsv("significant_genes_description_combined_brief_lfc_1.tsv") %>%
      select(gene, description, padj_voom, logFC_voom, description_source)
  ) %>%
  write_tsv("parasitemia_severity_overlap_with_description_lfc_1.tsv")
  