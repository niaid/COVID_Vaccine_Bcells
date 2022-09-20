library(tidyverse)

plot_dat <- readRDS("analysis_out/anova_results_clusters.rds")

plot_dat <- plot_dat %>% filter(!is.na(p.adjusted))

levels(plot_dat$clust_with_anno) <- paste0("C", levels(plot_dat$clust_with_anno))

plot_dat <- plot_dat %>%
  filter(clust_with_anno != "C8 NB")

pdf("plots/paper_figures/cluster_time_associations_anova_main_fig.pdf", width = 9, height = 2.5)
p <- plot_dat %>% 
  filter(term == "timepoint" & feature %in% c("cluster_fraction_total", "double_fraction")) %>%
  mutate(feature = factor(feature, levels = c("double_fraction", "cluster_fraction_total"))) %>%  
  ggplot(aes(x = clust_with_anno, y = feature)) +
  geom_point(aes(size = -log10(p.adjusted), alpha = p.adjusted < .05)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_viridis_c() +
  ylab("") +
  xlab("Cluster")
print(p)

dev.off()



