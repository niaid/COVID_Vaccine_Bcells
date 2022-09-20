library(tidyverse)
library(ggbeeswarm)


anova_dat <- readRDS("analysis_out/anova_results_clusters.rds")

cell_freq <- readRDS("analysis_out/FlowSOM/cluster_freq_and_antigen_specificity_summary.rds")

cluster_anno <- anova_dat %>%
  select(fs_meta_cluster, clust_with_anno) %>%
  distinct()
levels(cluster_anno$clust_with_anno) <- paste0("C", levels(cluster_anno$clust_with_anno))


cell_freq <- left_join(cell_freq, cluster_anno)

cell_freq <- cell_freq %>%
  mutate(Timepoint = timepoint)

source("workflow/paper_figures/util.convert_timepoints.R")
levels(cell_freq$Timepoint) <- new_tp(levels(cell_freq$timepoint))

#plot cells that change over time

keep_clusters1 <- anova_dat %>%
  filter(feature == "cluster_fraction_total" & p.adjusted < .05) %>%
  pull(clust_with_anno)

keep_clusters1 <- as.character(c(9, 13, 6, 4, 19, 27))
keep_clusters1 <- as.character(c(6, 4, 27, 9, 13, 19))
keep_clusters1_with_anno <- cluster_anno$clust_with_anno[match(keep_clusters1, cluster_anno$fs_meta_cluster)]

plot_dat1 <- cell_freq %>%
  filter(fs_meta_cluster %in% keep_clusters1) %>%
  mutate(clust_with_anno = factor(clust_with_anno, levels = keep_clusters1_with_anno))



p1 <- plot_dat1 %>%
  ggplot(aes(x = Timepoint, y = cluster_fraction_total)) +
  geom_boxplot(outlier.shape = NA, color = "red") +
  geom_beeswarm(alpha = .5, size = .5) +
  scale_y_continuous(trans='sqrt') +
  facet_wrap(~clust_with_anno, nrow = 1, scales = "free") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("plots/paper_figures/cluster_timecourse_fraction_total.pdf", height = 2, width = 9)
print(p1)
dev.off()

keep_clusters2 <- as.character(c(2, 4, 11, 9, 5, 3, 6, 13))
keep_clusters2_with_anno <- cluster_anno$clust_with_anno[match(keep_clusters2, cluster_anno$fs_meta_cluster)]

plot_dat2 <- cell_freq %>%
  filter(fs_meta_cluster %in% keep_clusters2) %>%
  mutate(clust_with_anno = factor(clust_with_anno, levels = keep_clusters2_with_anno))

p2 <- plot_dat2 %>%
  filter(fs_meta_cluster %in% keep_clusters2) %>%
  ggplot(aes(x = Timepoint, y = double_fraction)) +
  geom_boxplot(outlier.shape = NA, color = "red") +
  geom_beeswarm(alpha = .5, size = .5) +
  scale_y_continuous(trans='sqrt') +
  facet_wrap(~clust_with_anno, nrow = 2, scales = "free_y") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("plots/paper_figures/cluster_timecourse_fraction_rbds1plus.pdf", height = 3.2, width = 9)
print(p2)
dev.off()

