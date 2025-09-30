#!/usr/bin/Rscript 
#!/usr/bin/env R 4.3.1

# Packages
library(tidyverse)
library(data.table)
library(scales)
library(QCEWAS)
library(ggExtra)

# p-value threshold
n_snps <- read_tsv("fastLMM_ld.bim", col_names = FALSE) %>% nrow()
p_cutoff <- 0.05 / n_snps

# Manhattan
my_manhattan <- function(gwasResults){
  # Convert Chr to factor and arrange
  gwasResults <- gwasResults %>%
    mutate(Chr = as.factor(as.numeric(as.character(Chr)))) %>% 
    select(-SNP) %>%
    arrange(Chr, ChrPos)

  # Add SNP index
  gwasResults <- gwasResults %>%
    add_column(SNP = 1:nrow(.))

  # Compute chromosome label positions (midpoint of SNP index per Chr)
  chr_labels <- gwasResults %>%
    group_by(Chr) %>%
    summarise(center = mean(SNP), .groups = "drop")

  # Compute -log10(p) and size
  gwasResults <- gwasResults %>%
    mutate(logP = -log10(PValue),
           size = ifelse(logP > -log10(1e-4), scales::rescale(logP, to = c(1.5, 5)), 1))

  # Plot
  ggplot(gwasResults, aes(x = SNP, y = logP, colour = Chr, size = size)) +
    geom_point() + 
    geom_hline(yintercept = -log10(p_cutoff), linetype = 2) +               # Bonferroni line
    geom_hline(yintercept = -log10(1e-4), linetype = "dotted") +           # Suggestive line
    scale_colour_viridis_d() +
    scale_size_identity() +
    scale_x_continuous(
        breaks = chr_labels$center,
        labels = chr_labels$Chr
    ) +
    labs(
        x = "Chromosome",
        y = expression(-log[10]("P-value"))
    ) +
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "none"
    )
}


#QQ Plot
my_qq_plot <- function(gwasResults) {
    pvector <- gwasResults %>%
        filter(PValue > 0 & PValue < 1 & !is.na(PValue) & !is.nan(PValue) & is.finite(PValue)) %>%
        pull(PValue) 

    o = -log10(sort(pvector, decreasing = FALSE))
    e = -log10(ppoints(length(pvector)))

    tmp_plot <- as.data.frame(o) %>%
        cbind(
            as.data.frame(e)
        ) %>% 
        ggplot(aes(x = e, y = o)) +
            geom_point(colour = "#1F968BFF") + 
            xlab(expression(Expected ~ ~-log[10](italic(p)))) +
            ylab(expression(Observed ~ ~-log[10](italic(p)))) +
            geom_abline(colour = "#440154FF")  +
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
    )

    ggMarginal(tmp_plot, type = "histogram", bins = 50)
}

#lambda
my_lambda <- function(gwasResults){
    gwasResults %>% 
        select(PValue) %>% 
        na.omit() %>%
        pull() %>%
        P_lambda() 
}

# Pulling it all together
my_gwas_outputs <- function(file_path){
    gwasResults <- read_tsv(file_path) %>%
        mutate(Chr = as.factor(Chr), 
            SNP = str_remove(SNP, "ordered_PKNH_"),
            SNP = str_remove(SNP, "new"),
            SNP = as.factor(SNP)) 

    man_ggplot <- my_manhattan(gwasResults)
    ggsave(paste0("manhattan_", str_remove(file_path, ".tsv"), ".png"), dpi = 300, width = 20, man_ggplot)

    qq_ggplot <- my_qq_plot(gwasResults)
    ggsave(paste0("qq_", str_remove(file_path, ".tsv"), ".png"), dpi = 300, qq_ggplot)

    my_lambda(gwasResults) %>%
        write.table(paste0("lambda_", str_remove(file_path, ".tsv"), ".txt"), col.names = FALSE, sep = "\t") 

    gwasResults %>% 
        filter(PValue < p_cutoff) %>%
        write.table(paste0("sig_snps_", str_remove(file_path, ".tsv"), ".txt"), col.names = TRUE, sep = "\t") 

    gwasResults %>% 
        filter(PValue < 1e-4) %>%
        write.table(paste0("suggestive_snps_", str_remove(file_path, ".tsv"), ".txt"), col.names = TRUE, row.names = FALSE, sep = "\t")
}

gwasResults_path_list <- list.files(pattern="^testing") 

for (i in gwasResults_path_list){
    my_gwas_outputs(paste0(i))
}
