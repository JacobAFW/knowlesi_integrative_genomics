#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(scales)
library(QCEWAS)
library(ggExtra)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript plot_fastlmm_results.R <result_tsv> <bim_file>")
file_path <- args[1]
bim_path <- args[2]

# P-value threshold
n_snps <- fread(bim_path, header = FALSE) %>% nrow()
p_cutoff <- 0.05 / n_snps

# Manhattan plot
my_manhattan <- function(gwasResults) {
  gwasResults %>%
    mutate(Chr = as.factor(as.numeric(as.character(Chr)))) %>%
    arrange(Chr, ChrPos) %>%
    mutate(SNP = row_number(),
           logP = -log10(PValue),
           size = ifelse(logP > -log10(p_cutoff), rescale(logP, to = c(1.5, 5)), 1)) -> gwas

  chr_labels <- gwas %>%
    group_by(Chr) %>%
    summarise(center = mean(SNP), .groups = "drop")

  ggplot(gwas, aes(x = SNP, y = logP, colour = Chr, size = size)) +
    geom_point() +
    geom_hline(yintercept = -log10(p_cutoff), linetype = 2) +
    geom_hline(yintercept = -log10(1e-4), linetype = "dotted") +
    scale_colour_viridis_d() +
    scale_size_identity() +
    scale_x_continuous(breaks = chr_labels$center, labels = chr_labels$Chr) +
    labs(x = "Chromosome", y = expression(-log[10]("P-value"))) +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.position = "none")
}

# QQ plot
my_qq_plot <- function(gwasResults) {
  pvector <- gwasResults %>%
    filter(PValue > 0 & PValue < 1 & !is.na(PValue)) %>%
    pull(PValue)

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
my_lambda <- function(gwasResults) {
  gwasResults %>% 
    select(PValue) %>% 
    drop_na() %>% 
    pull() %>% 
    P_lambda()
}

# Main function
gwasResults <- read_tsv(file_path) %>%
  mutate(Chr = as.factor(Chr),
         SNP = str_remove(SNP, "ordered_PKNH_|new"),
         SNP = as.factor(SNP))

gene_id <- str_remove(basename(file_path), ".tsv")

outdir <- dirname(file_path)

ggsave(file.path(outdir, paste0("manhattan_", gene_id, ".png")),
       my_manhattan(gwasResults), dpi = 300, width = 20)

ggsave(file.path(outdir, paste0("qq_", gene_id, ".png")),
       my_qq_plot(gwasResults), dpi = 300)

write.table(my_lambda(gwasResults),
            file.path(outdir, paste0("lambda_", gene_id, ".txt")),
            col.names = FALSE, sep = "\t", quote = FALSE)

gwasResults %>%
  filter(PValue < p_cutoff) %>%
  write_tsv(file.path(outdir, paste0("sig_snps_", gene_id, ".txt")))

gwasResults %>%
  filter(PValue < 1e-4) %>%
  write_tsv(file.path(outdir, paste0("suggestive_snps_", gene_id, ".txt")))