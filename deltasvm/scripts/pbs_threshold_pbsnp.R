Args <- commandArgs(TRUE)
pbs_input <- Args[1]
deltasvm_input <- Args[2]
output <- Args[3]
pbsnp_ratio <- as.numeric(Args[4]) / 100

library(tidyverse)

pbs <- read_tsv(pbs_input, col_names = c("id", "obs_auc", "obs_pval", "pbs", "pbs_pval", "pbs_fdr"))  %>%
  separate(id, c("tf", "snp"), extra = "drop", fill = "right", sep = ":") 
deltasvm <- read_tsv(deltasvm_input)  %>%
  separate(exp, c("tf", "domain", "cycle", "protein", "well"), sep = '_', remove = F)

scores  <- inner_join(deltasvm, pbs, by = c("tf", "snp")) %>% 
  group_by(snp) %>% 
  filter(pbs_pval==min(pbs_pval)) %>%
  ungroup()

pbsnp_scores <- abs(scores %>%
  filter(pbs_pval < 0.01) %>%
  pull(deltasvm))

threshold <- scores %>%
  select(exp, kmer) %>%
  distinct() %>%
  separate(exp, c("tf", "domain", "cycle", "protein", "well"), sep = '_', remove = F) %>%
  mutate(threshold = quantile(pbsnp_scores, pbsnp_ratio))

write_tsv(threshold, output)
