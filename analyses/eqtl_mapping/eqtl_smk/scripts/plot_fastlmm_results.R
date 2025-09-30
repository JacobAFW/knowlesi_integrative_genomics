#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(scales)
library(QCEWAS)
library(ggExtra)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_fastlmm_results.R <result_tsv> <bim_file> [cis|trans] [window_kb]")

file_path <- args[1]
bim_path  <- args[2]
mode      <- ifelse(length(args) >= 3, tolower(args[3]), "cis")
win_kb    <- ifelse(length(args) >= 4, as.numeric(args[4]), 10)
if (!mode %in% c("cis","trans")) stop("Mode must be 'cis' or 'trans'")
win_bp    <- as.integer(win_kb * 1000)

# Parse gene ID (used for filenames)
gene_id_clean <- str_remove(basename(file_path), "\\.tsv$") %>% str_trim() %>% str_squish()
outdir <- dirname(file_path)

# Load eQTL/GWAS results
gwasResults_all <- read_tsv(file_path, show_col_types = FALSE) %>%
  mutate(
    Chr = as.factor(Chr),
    SNP = str_remove(SNP, "ordered_PKNH_|new"),
    SNP = as.factor(SNP)
  )

# Default to all results; optionally subset cis-window
gwasResults <- gwasResults_all
mode_tag <- mode

if (mode == "cis") {
  gene_pos <- read_tsv("/g/data/pq84/malaria/Pk_trancsriptomics/outputs/analyses/eqtl/gene_positions.tsv",
                       show_col_types = FALSE) %>%
    mutate(trinity = str_trim(trinity) %>% str_squish())

  region <- gene_pos %>%
    filter(trinity == gene_id_clean) %>%
    mutate(start = pmax(start - win_bp, 0L),
           end   = end + win_bp)

  if (nrow(region) == 0) stop("Gene not found in gene_positions.tsv")

  chrom <- region$chr %>%
    str_remove("ordered_PKNH_") %>% str_remove("_v2") %>% str_remove("^0+")
  cis_start <- region$start
  cis_end   <- region$end

  gwasResults <- gwasResults_all %>%
    filter(Chr == chrom, ChrPos >= cis_start, ChrPos <= cis_end)
}

# Thresholds based on the analysis set (cis window or all)
n_snps <- nrow(gwasResults)
if (mode == "cis" && n_snps == 0) stop("No SNPs in cis-window for this gene.")

p_cutoff <- 0.05 / n_snps
suggestive_cutoff <- min(1 / n_snps, 0.05)  # guard against >0.05 suggestives

# Manhattan plot function (plots genome-wide for context)
my_manhattan <- function(gwasResults_all, p_cutoff, suggestive_cutoff) {
  gwasResults_all %>%
    mutate(Chr = as.factor(as.numeric(as.character(Chr)))) %>%
    arrange(Chr, ChrPos) %>%
    mutate(SNP = row_number(),
           logP = -log10(PValue),
           size = ifelse(logP > -log10(p_cutoff), rescale(logP, to = c(1.5, 5)), 1)) -> gwas

  chr_labels <- gwas %>% group_by(Chr) %>% summarise(center = mean(SNP), .groups = "drop")

  ggplot(gwas, aes(x = SNP, y = logP, colour = Chr, size = size)) +
    geom_point() +
    geom_hline(yintercept = -log10(p_cutoff), linetype = 2) +
    geom_hline(yintercept = -log10(suggestive_cutoff), linetype = "dotted") +
    scale_colour_viridis_d() +
    scale_size_identity() +
    scale_x_continuous(breaks = chr_labels$center, labels = chr_labels$Chr) +
    labs(x = "Chromosome", y = expression(-log[10]("P-value"))) +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.position = "none")
}

# QQ plot
my_qq_plot <- function(gwasResults_all) {
  pvector <- gwasResults_all %>% filter(PValue > 0 & PValue < 1 & !is.na(PValue)) %>% pull(PValue)
  if (length(pvector) < 2) return(ggplot() + theme_void())
  o <- -log10(sort(pvector))
  e <- -log10(ppoints(length(pvector)))
  tmp <- data.frame(e = e, o = o)
  ggMarginal(
    ggplot(tmp, aes(e, o)) +
      geom_point(color = "#1F968BFF") +
      geom_abline(color = "#440154FF") +
      xlab(expression(Expected ~ -log[10](italic(p)))) +
      ylab(expression(Observed ~ -log[10](italic(p)))) +
      theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14)),
    type = "histogram", bins = 50
  )
}

# Lambda
my_lambda <- function(gwasResults_all) {
  gwasResults_all %>% select(PValue) %>% drop_na() %>% pull() %>% P_lambda()
}

# Hits (based on the analysis set)
sig_hits <- gwasResults %>% filter(PValue < p_cutoff)
sug_hits <- gwasResults %>% filter(PValue < 0.05, PValue < suggestive_cutoff)

has_hits <- any(gwasResults$PValue < p_cutoff)

# Outputs
manhattan_fp <- file.path(outdir, paste0("manhattan_", gene_id_clean, "_", ".png"))
qq_fp       <- file.path(outdir, paste0("qq_", gene_id_clean, "_", ".png"))
lambda_fp   <- file.path(outdir, paste0("lambda_", gene_id_clean, "_", ".txt"))
sig_fp      <- file.path(outdir, paste0("sig_snps_", gene_id_clean, "_", ".txt"))
sug_fp      <- file.path(outdir, paste0("suggestive_snps_", gene_id_clean, "_", ".txt"))

if (has_hits) {
  ggsave(manhattan_fp, my_manhattan(gwasResults_all, p_cutoff, suggestive_cutoff), dpi = 300, width = 20)
  ggsave(qq_fp,        my_qq_plot(gwasResults_all), dpi = 300)
  write.table(my_lambda(gwasResults_all), lambda_fp, col.names = FALSE, sep = "\t", quote = FALSE)
}

write_tsv(sig_hits, sig_fp)
write_tsv(sug_hits, sug_fp)