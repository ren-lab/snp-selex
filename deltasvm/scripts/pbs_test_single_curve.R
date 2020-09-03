Args <- commandArgs(TRUE)
pbs_input <- Args[1]
deltasvm_input <- Args[2]
output <- Args[3]

library(PRROC)
library(tidyverse)

pbs <- read_tsv(pbs_input, col_names = c("id", "obs_auc", "obs_pval", "pbs", "pbs_pval", "pbs_fdr"))  %>%
  separate(id, c("tf", "snp"), extra = "drop", fill = "right", sep = ":") 
deltasvm <- read_tsv(deltasvm_input)  %>%
  separate(exp, c("tf", "domain", "cycle", "protein", "well"), sep = '_', remove = F)

scores  <- inner_join(deltasvm, pbs, by = c("tf", "snp")) %>% 
  group_by(snp) %>% 
  filter(pbs_pval==min(pbs_pval)) %>%
  ungroup()

corr <- scores %>%
  group_by(exp, kmer) %>%
  summarise(pcc=cor(deltasvm, pbs), pcc_pval=cor.test(deltasvm, pbs)$p.value, scc=cor(deltasvm, pbs, method = "spearman"), scc_pval=cor.test(deltasvm, pbs, method = "spearman")$p.value)  %>%
  separate(exp, c("tf", "domain", "cycle", "protein", "well"), sep = '_', remove = F)
write_tsv(corr, paste0(output, ".corr.tsv"))

pdf(paste0(output, ".corr.pdf"))
scores %>%
  ggplot(aes(pbs, deltasvm)) +
  geom_hex(aes(fill=log10(..count..)), bins = 100) +
  geom_smooth(method ='lm', colour="red", se = FALSE) +
  facet_wrap(~tf) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=20)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("PBS") + 
  xlab("deltaSVM") +
  guides(fill=guide_colorbar(title="Density")) 
dev.off()

scores <-  scores %>%
  mutate(pbsnp = ifelse(pbs_pval < 0.01, 1, 0), deltasvm = ifelse(pbs>0, -deltasvm, deltasvm))  %>%
  filter(pbsnp == 1 | pbs_pval > 0.5)

auc <- scores %>%
  group_by(exp, kmer) %>%
  summarize(
    aupr = pr.curve(scores.class0 = deltasvm, weights.class0 = pbsnp)$auc.integral,
    auroc = roc.curve(scores.class0 = deltasvm, weights.class0 = pbsnp)$auc,
    count = sum(pbsnp)
  ) %>%
  separate(exp, c("tf", "domain", "cycle", "protein", "well"), sep = '_', remove = F) 

write_tsv(auc, paste0(output, ".auc.tsv"))

wroc <- roc.curve(scores.class0 = scores$deltasvm, weights.class0 = scores$pbsnp, curve = T)
wpr <- pr.curve(scores.class0 = scores$deltasvm, weights.class0 = scores$pbsnp, curve = T)

pdf(paste0(output, ".curve.pdf"))
plot(wroc)
plot(wpr)
dev.off()

wroc_c <- wroc$curve
colnames(wroc_c) <- c("x", "y", "c")
as_tibble(wroc_c) %>%
  mutate(tf = scores$tf[1]) %>%
  write_tsv(paste0(output, ".roc.tsv"))

wpr_c <- wpr$curve
colnames(wpr_c) <- c("x", "y", "c")
as_tibble(wpr_c) %>%
  mutate(tf = scores$tf[1]) %>%
  write_tsv(paste0(output, ".pr.tsv"))

set.seed(0)
scores_bal <- NULL

tmp_pos <- scores %>% 
  filter(pbsnp == 1)
tmp_neg <- scores %>% 
  filter(pbsnp == 0) 
if (nrow(tmp_neg) < nrow(tmp_pos)) next
tmp_neg <- tmp_neg %>%
  sample_n(nrow(tmp_pos))
scores_bal <- bind_rows(scores_bal, tmp_pos, tmp_neg)

auc_bal <- scores_bal %>%
  group_by(exp, kmer) %>%
  summarize(
    aupr = pr.curve(scores.class0 = deltasvm, weights.class0 = pbsnp)$auc.integral,
    auroc = roc.curve(scores.class0 = deltasvm, weights.class0 = pbsnp)$auc,
    count = sum(pbsnp)
  ) %>%
  separate(exp, c("tf", "domain", "cycle", "protein", "well"), sep = '_', remove = F) 

write_tsv(auc_bal, paste0(output, ".auc_bal.tsv"))

pdf(paste0(output, ".curve_bal.pdf"))
t <- scores_bal
plot(roc.curve(scores.class0 = scores_bal$deltasvm, weights.class0 = scores_bal$pbsnp, curve = T))
plot(pr.curve(scores.class0 = scores_bal$deltasvm, weights.class0 = scores_bal$pbsnp, curve = T))
dev.off()
