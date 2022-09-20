library(tidyverse)
library(pheatmap)
library(Hmisc)
library(corrplot)

ab_titer_joined <- readRDS("ab_data/reshaped/ab_titer_joined.rds")

spread_dat <- ab_titer_joined %>%
  filter(timepoint %in% c("D2_4weeks")) %>%
  select(`Code ID`, S2P_num, RBD_num, Readout, timepoint, date_samp) %>%
  gather(key = antigen, value = titer, -c(`Code ID`, Readout, timepoint, date_samp)) %>%
  mutate(antigen = gsub("_num", "", antigen)) %>%
  mutate(antigen = gsub("S2P", "S-2P", antigen)) %>%
  mutate(antigen_isotype = paste(antigen, Readout, sep = " ")) %>%
  group_by(`Code ID`, antigen_isotype) %>%
  summarise(titer = mean(titer, na.rm = TRUE)) %>%
  mutate(titer = log(titer)) %>%
  spread(key = antigen_isotype, value = titer) %>%
  ungroup()

ab_mat <- spread_dat %>%
  select(-`Code ID`) %>%
  as.matrix()

#corr <- cor(ab_mat, use = "pairwise.complete.obs", method = "spearman")
corr <- rcorr(ab_mat, type = "spearman")

dir.create("plots/paper_figures/")
pdf("plots/paper_figures/ab_titer_cor.pdf", height = 4, width = 4)
col <- colorRampPalette(rev(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")))
corrplot(corr$r, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = corr$P, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         cl.cex = .4
)

corrplot(corr$r, method="color", col=col(200),  
         type="full", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = corr$P, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=TRUE,
         cl.cex = .4
)
dev.off()

