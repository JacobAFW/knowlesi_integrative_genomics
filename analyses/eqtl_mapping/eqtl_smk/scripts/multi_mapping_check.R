# eQTL SNP Overlap and Multi-mapping Check

library(tidyverse)

# Set paths
results_dir <- "results/eqtl_results"
eqtl_files <- list.files(results_dir, pattern = "\\.tsv$", full.names = TRUE)

# Load all eQTL result files into a single data frame
all_eqtl <- eqtl_files %>%
  set_names(~ str_remove(basename(.x), "\\.tsv$")) %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE), .id = "transcript_id")

# List of SNPs associated with TRINITY_DN8677_c0_g1
snp_hits <- c(
  "14_v2:78033:T:A:",
  "10_v2:1055355:T:*:",
  "09_v2:168307:C:A:",
  "12_v2:1306228:A:G:",
  "13_v2:2053530:A:G:",
  "12_v2:2719576:T:C:"
)

# Filter for SNPs of interest
shared_snp_hits <- all_eqtl %>%
  filter(SNP %in% snp_hits) %>%
  arrange(SNP, PValue)

# Save output for later review
write_tsv(shared_snp_hits, "shared_snp_across_transcripts.tsv")

# Optional: summarize how many transcripts each SNP is linked to
shared_snp_hits %>%
  count(SNP, name = "num_transcripts") %>%
  arrange(desc(num_transcripts)) %>%
  write_tsv("snp_transcript_counts.tsv")
