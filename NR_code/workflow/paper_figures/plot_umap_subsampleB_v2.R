library(tidyverse)
library(readxl)
library(ggrepel)
library(ggthemes)

source("workflow/paper_figures/util.high_contrast_palette.R")
source("workflow/paper_figures/util.convert_timepoints.R")

dir.create("plots/paper_figures/umap")

umap_in <- readRDS("analysis_out/umap/umap_subsample_listB.rds")

umap_df <- umap_in$umap_df


clust_anno <- read_excel("Clusters v050621 Annotated 051321.xlsx", n_max = 31)

clust_anno <- clust_anno %>%
  rename(fs_meta_cluster = Cluster) %>%
  mutate(clust_with_anno = paste(fs_meta_cluster, Population)) %>%
  mutate(fs_meta_cluster = factor(fs_meta_cluster, levels = levels(umap_df$fs_meta_cluster)))


#https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
pop_factor_levels <- levels(factor(clust_anno$Population %>% .[.!= "NB"]))
color_palette_pop <- ggthemes::colorblind_pal()(length(pop_factor_levels))
names(color_palette_pop) <- pop_factor_levels

umap_df <- umap_df %>%
  mutate(timepoint = umap_in$cell_meta_sub$timepoint) %>%
  mutate(Timepoint = timepoint)


levels(umap_df$Timepoint) <- new_tp(levels(umap_df$Timepoint))

umap_df <- umap_df %>% left_join(clust_anno)

umap_df$clust_with_anno <- factor(umap_df$clust_with_anno)

levels(umap_df$clust_with_anno) <- paste0("C", levels(umap_df$clust_with_anno))
levels(umap_df$fs_meta_cluster) <- paste0("C", levels(umap_df$fs_meta_cluster))


rbd_s1_dat <- readRDS("processed_data/spike_rbd_dat.rds")
rbd_s1_dat <- rbd_s1_dat[umap_in$keep_ix,]
table(rbd_s1_dat$gated_rbdplus, umap_df$rbd_plus)


umap_df <- umap_df %>%
  filter(clust_with_anno != "C8 NB")

clust_summ_dat <- umap_df %>%
  group_by(clust_with_anno, fs_meta_cluster) %>%
  summarise(x = median(x), y = median(y),
            pcnt_double = sum(rbd_plus & s1_plus) / n())

clust_summ_timepoint_dat <- umap_df %>%
  group_by(clust_with_anno, fs_meta_cluster, Timepoint) %>%
  summarise(x = median(x), y = median(y),
            pcnt_double = sum(rbd_plus & s1_plus) / n(),
            n_double = sum(rbd_plus & s1_plus)) %>%
  ungroup()

population_summ_dat <- umap_df %>%
  group_by(Population) %>%
  summarise(x = mean(x), y = median(y))

population_summ_timepoint_dat <- umap_df %>%
  group_by(Population, Timepoint) %>%
  summarise(x = median(x), y = median(y),
            pcnt_double = sum(rbd_plus & s1_plus) / n(),
            n_double = sum(rbd_plus & s1_plus)) %>%
  ungroup()

memory_pop_kmeans <- umap_df %>%
  filter(Population == "MBC") %>%
  select(x,y) %>%
  as.matrix() %>%
  kmeans(centers = 3)

memory_pop_kmeans_centers_pos <- 
  memory_pop_kmeans$centers %>%
  `colnames<-`(c("x", "y")) %>%
  as.data.frame()
  
color_palette <- high_contrast_palette()[1:length(levels(umap_df$clust_with_anno))]
names(color_palette) <- levels(umap_df$clust_with_anno)

png("plots/paper_figures/umap/umap_w_clusters_super_small.png", width = 1900, height = 1500, res = 300)
p <- ggplot(umap_df, aes(x = x, y = y)) +
  geom_point(aes(color = clust_with_anno), size = .001) +
  geom_label_repel(data = clust_summ_dat, aes(label = as.character(fs_meta_cluster)), size = 2.5) +
  theme_void() + 
  scale_color_manual(values = color_palette) +
  guides(colour = guide_legend(override.aes = list(size=6)))
  
print(p)
dev.off()

#file.remove("plots/paper_figures/umap/umap_w_annotations.png")
png("plots/paper_figures/umap/umap_w_annotations.png", width = 1800, height = 1500, res = 300)
p <- ggplot(umap_df, aes(x = x, y = y)) +
  geom_point(aes(color = Population), size = .001) +
  geom_label_repel(data = population_summ_dat %>% filter(Population != "MBC"), aes(label = as.character(Population))) +
  geom_label_repel(data = memory_pop_kmeans_centers_pos, aes(label = "MBC")) +
  scale_color_manual(values = color_palette_pop) +
  theme_void() + guides(colour = guide_legend(override.aes = list(size=10)))
print(p)
dev.off()


keep_tp <- c("v1D0", "v1D7", "v1D10", "v1D14", "v2D0", "v2D7", "v2D10", "v2D28")
plotlist <- lapply(keep_tp, function(tp){
  tmp_doublepos_dat <- umap_df %>% filter(rbd_plus & s1_plus & 
                                            Timepoint == tp)
  
  highlight_clusters <- clust_summ_timepoint_dat %>% 
    filter(pcnt_double > .003, Timepoint == tp) %>%
    pull(clust_with_anno)
  
  tmp_cluster_pos_dat <- umap_df %>%
    filter(clust_with_anno %in% highlight_clusters)
  
  p <- ggplot(umap_df, aes(x = x, y = y)) +
    stat_summary_2d(aes(z = x), fun = function(x) 1, bins = 300) +
    scale_fill_gradient(low = "grey", high = "grey") +
    geom_point(data = tmp_doublepos_dat, fill="darkolivegreen1", color = "blue4", size = 1, shape = 23, stroke = 1) +
    geom_label_repel(data = clust_summ_dat %>% filter(clust_with_anno %in% highlight_clusters), aes(label = as.character(fs_meta_cluster))) +
    theme_classic()+
    theme(legend.position = "none", axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 18, face = "bold")) +
    ggtitle(tp) + 
    guides(colour = guide_legend(override.aes = list(size=6)))

})
timeplot1 <- cowplot::plot_grid(plotlist = plotlist, nrow = 2)

png("plots/paper_figures/umap/umap_rbds1_plus_timepoints.png", width = 4000, height = 2000, res = 300)
print(timeplot1)
dev.off()
