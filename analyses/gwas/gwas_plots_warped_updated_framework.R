#! /usr/bin/Rscript
#!/usr/bin/env R 4.3.1

# Packages
library(tidyverse)
library(data.table)
library(scales)
library(QCEWAS)
library(janitor)
library(ggExtra)

# Manhattan
my_manhattan <- function(gwasResults, p_cutoff){
  # Convert Chr to factor and arrange
  gwasResults <- gwasResults %>%
    mutate(Chr = as.factor(as.numeric(as.character(Chr)))) %>% 
    #select(-SNP) %>%
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
           size = ifelse(logP > -log10(p_cutoff), scales::rescale(logP, to = c(1.5, 5)), 1))

  # Plot
    ggplot(gwasResults, aes(x = SNP, y = logP, colour = Chr, size = size)) +
    geom_point() + 
    geom_hline(yintercept = -log10(p_cutoff), linetype = 2) + # Bonferroni line
    geom_hline(yintercept = -log10(1e-4), linetype = "dotted") + # Suggestive line
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
    gwasResults <- read.table(file_path, sep = "\t") %>% 
                        row_to_names(1) %>%
                        #slice(-1) %>%
                        mutate(PValue = as.numeric(PValue))  %>%
                        mutate(Chr = as.factor(Chr),
                        ChrPos = as.factor(ChrPos)) 

    # LD-based p-value threshold
    p_cutoff <- read.table("warpedLMM_covar_ld.prune.in") %>%
    nrow()
    p_cutoff <- 0.05/p_cutoff

    # Manhattan
    man_ggplot <- my_manhattan(gwasResults, p_cutoff)
    ggsave(paste0(str_remove(file_path, "gwas_results.feather"), "manhattan_gwas_results_warped_", str_remove(file_path, "\\./") %>% str_remove("/.*"), ".png"), width = 20, dpi = 300, man_ggplot)

    # QQ
    qq_ggplot <- my_qq_plot(gwasResults)
    ggsave(paste0(str_remove(file_path, "gwas_results.feather"), "qq_gwas_results_warped_", str_remove(file_path, "\\./") %>% str_remove("/.*"), ".png"), dpi = 300, qq_ggplot)

    #Lambda
    my_lambda(gwasResults) %>%
        write.table(paste0(str_remove(file_path, "gwas_results.feather"), "lambda_gwas_results_warped_", str_remove(file_path, "\\./") %>% str_remove("/.*"), ".txt"), col.names = FALSE, sep = "\t")
    
    # Sig SNPs
    gwasResults %>% 
        filter(PValue < p_cutoff) %>%
        write.table(paste0(str_remove(file_path, "gwas_results.feather"), "sig_snps_", str_remove(file_path, "\\./") %>% str_remove("/.*"), ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
    # Suggestive SNPs
    gwasResults %>% 
        filter(PValue < 1e-4) %>%
        write.table(paste0(str_remove(file_path, "gwas_results.feather"), "suggestive_snps_", str_remove(file_path, "\\./") %>% str_remove("/.*"), ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
}


gwasResults_path_list <- list.files(pattern = "gwas_results.feather", recursive = TRUE, full.names = TRUE)

for (i in gwasResults_path_list){
    my_gwas_outputs(paste0(i))
}


# Top 10 SNPs from each model

top_10_warped <- data.frame(matrix(ncol = 6, nrow = 0))
names <- c("Chr", "ChrPos", "Dist", "PValue", "warped", "Model")
colnames(top_10_warped) <- names

for (i in gwasResults_path_list){
    top_10_warped <- top_10_warped %>%
        rbind(
            read.table(i, sep = "\t") %>% 
                row_to_names(1) %>%
                #arrange(PValue) %>%
                slice(1:10) %>%
                mutate(Model = paste0(str_remove(i, "./LD/") %>% str_remove("/gwas.*")))
        )
}

top_10_warped %>%
    write_tsv("top_10_snps_of_each_model_warped.tsv")

# All lambda values

lambda_path_list <- list.files( pattern = "^lambda", recursive = TRUE, full.names = TRUE)
all_lambdas <- data.frame(matrix(ncol = 2, nrow = 0))
names <- c("model", "lambda")
colnames(all_lambdas) <- names

for (i in lambda_path_list){
    all_lambdas <- all_lambdas %>%
        rbind(
            read.table(i, sep = "\t") %>% 
                rename(model = V1, lambda = V2) %>%
                mutate(model = paste0(str_remove(i, "./") %>% str_remove("/lamda.*"))) %>%
                separate(model,  into = c("MAF", "combination"), sep = "/")
        )
}

all_lambdas %>%
    write_tsv("all_lambdas.tsv")


