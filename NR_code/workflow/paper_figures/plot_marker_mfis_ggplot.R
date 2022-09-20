library(tidyverse)
library(FlowSOM)
library(pheatmap)
library(cowplot)
library(readxl)

fSOM <- readRDS("analysis_out/FlowSOM/fsom.rds")

param_dat <- readRDS("processed_data/param_dat.rds")

freq_dat <- readRDS("analysis_out/FlowSOM/cluster_freq_and_antigen_specificity_summary.rds")

fs_meta_cluster <- fSOM$metaclustering[fSOM$FlowSOM$map$mapping[, 1]]
sort(table(fs_meta_cluster))  

clust_anno <- read_excel("Clusters v050621 Annotated 051321.xlsx", n_max = 31)
clust_anno <- clust_anno %>%
  rename(fs_meta_cluster = Cluster) %>%
  mutate(fs_meta_cluster = factor(fs_meta_cluster, levels = levels(freq_dat$fs_meta_cluster)))


meta_mfi <- sapply(unique(fs_meta_cluster), function(clust){
    apply(fSOM$FlowSOM$data[fs_meta_cluster == clust, ], 2, median)
  }) %>% t()
rownames(meta_mfi) <- unique(fs_meta_cluster)
cnames <- param_dat$desc[match(colnames(meta_mfi), param_dat$name)] 
colnames(meta_mfi) <- gsub("Fix ", "", cnames)

meta_mfi <- meta_mfi[, !is.na(colnames(meta_mfi))]

colnames(meta_mfi) <- gsub(" Fix", "", colnames(meta_mfi))

meta_mfi <- meta_mfi[, setdiff(colnames(meta_mfi), c("S1", "RBD"))]

plot_dat <- meta_mfi %>%
  as.data.frame() %>%
  rownames_to_column(var = "cluster") %>%
  gather(key = "feature", value = "mfi", -cluster)

plot_dat <- plot_dat %>% 
  left_join(rename(clust_anno, cluster = fs_meta_cluster))

hc_row <- meta_mfi %>% 
  dist() %>% hclust()
plot_dat$cluster_ordered <- factor(plot_dat$cluster, levels = hc_row$labels[hc_row$order])

hc_col <- meta_mfi %>% 
  t() %>%
  dist() %>% hclust()

plot_dat$feature_ordered <- factor(plot_dat$feature, levels = hc_col$labels[hc_col$order])

plot_dat <- plot_dat %>%
  mutate(feature_type = ifelse(startsWith(feature, "Ig"), "Ig", "CD"))

freq_summ <- freq_dat %>%
  left_join(clust_anno) %>%
  group_by(fs_meta_cluster, Population) %>%
  summarise(n_cells = sum(n_cells_in_cluster),
            S1 = sum(s1_count) / sum(n_cells_in_cluster),
            RBD = sum(rbd_count) / sum(n_cells_in_cluster),
            `RBD+S1+` = sum(double_count) / sum(n_cells_in_cluster)) %>%
  ungroup()
  

freq_summ$cluster_ordered <- factor(freq_summ$fs_meta_cluster, levels = hc_row$labels[hc_row$order])

pop_factor_levels <- levels(factor(clust_anno$Population %>% .[.!= "NB"]))
color_palette_pop <- ggthemes::colorblind_pal()(length(pop_factor_levels))
names(color_palette_pop) <- pop_factor_levels

levels(plot_dat$cluster_ordered) <- paste0("C", levels(plot_dat$cluster_ordered))
levels(freq_summ$cluster_ordered) <- paste0("C", levels(freq_summ$cluster_ordered))

freq_summ <- freq_summ %>%
  filter(Population != "NB")

plot_dat <- plot_dat %>%
  filter(Population != "NB")

freq_dat_long <- freq_summ %>%
  select(cluster_ordered, S1, RBD, `RBD+S1+`, Population) %>%
  gather(key = specificity, value = fraction, -c(cluster_ordered, Population))

p0 <- 
  ggplot(plot_dat, aes(x = 1, y = cluster_ordered)) +
  geom_tile(aes(fill = Population)) +
  xlab("") + ylab("Cluster") +
  facet_grid(3 ~ 1, scales = "free", space = "free") +
  scale_fill_manual(values = color_palette_pop) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background.y = element_blank(), strip.text.y = element_blank(),
        legend.position = "left",
        strip.text = element_blank(),
        panel.spacing.y=unit(0, "lines"))

p1 <- ggplot(plot_dat, aes(x = feature_ordered, y = cluster_ordered)) +
  geom_tile(aes(fill = mfi)) +
  #facet_grid(Population ~ feature_type, scales = "free", space = "free") +
  facet_grid(3 ~ feature_type, scales = "free", space = "free") +
  scale_fill_viridis_c() +
  xlab("Feature") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background.y = element_blank(), strip.text.y = element_blank(),
        legend.position = "top",
        panel.spacing.y=unit(0, "lines"))

p2_v2 <- ggplot(freq_dat_long %>% filter(specificity == "RBD+S1+"), aes(x = factor(specificity), y = cluster_ordered)) +
  geom_tile(aes(fill = fraction)) + theme_minimal() + 
  facet_grid(2 ~ 1, scales = "free", space = "free") +
  scale_fill_viridis_c(option = "B", trans = "sqrt") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(size=8), legend.position =  "top",
        legend.text=element_text(size=4),
        panel.spacing.y=unit(0, "lines")) +
  xlab("")

p3 <- ggplot(freq_summ, aes(x = n_cells, y = cluster_ordered)) +
  geom_col() + theme_minimal() + 
  facet_grid(2 ~ 1, scales = "free", space = "free") +
  #ggtitle("Ag Specificity") +
  scale_x_continuous(trans='log10') +
  xlab("Number of cells") +
  theme(#axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(size=8),
        panel.spacing.y=unit(0, "lines"))

pdf("plots/paper_figures/FlowSOM_cluster_MFI_heatmap.pdf", height = 5, width = 7)
p <- plot_grid(p0, p1, p2_v2, p3, align = "h", nrow = 1, axis = "tb", rel_widths = c(.38, .6, .08, .37))
print(p)
dev.off()

