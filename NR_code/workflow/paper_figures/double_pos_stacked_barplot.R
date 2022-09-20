library(tidyverse)
library(readxl)

cell_freq <- readRDS("analysis_out/FlowSOM/cluster_freq_and_antigen_specificity_summary.rds")

clust_anno <- read_excel("Clusters v050621 Annotated 051321.xlsx", n_max = 31)

clust_anno <- clust_anno %>%
  rename(fs_meta_cluster = Cluster) %>%
  mutate(clust_with_anno = paste(fs_meta_cluster, Population)) %>%
  mutate(fs_meta_cluster = factor(fs_meta_cluster, levels = levels(cell_freq$fs_meta_cluster)))

levels(clust_anno$clust_with_anno) <- paste0("C", levels(clust_anno$clust_with_anno))


cell_freq <- left_join(cell_freq, clust_anno)

cell_freq <- cell_freq %>%
  mutate(Timepoint = timepoint)

source("workflow/paper_figures/util.convert_timepoints.R")
levels(cell_freq$Timepoint) <- new_tp(levels(cell_freq$timepoint))

timepoint_summary_dat <- cell_freq %>%
  group_by(Timepoint, clust_with_anno, fs_meta_cluster) %>%
  summarise(n_double = sum(double_count),
            sum_total = sum(n_cd19_total),
            double_div_cd19 = sum(double_count) / sum(n_cd19_total)) %>% ungroup()


timepoint_summary_dat <- timepoint_summary_dat %>%
  mutate(clust_with_anno = factor(clust_with_anno))

levels(timepoint_summary_dat$clust_with_anno) <- paste0("C", levels(timepoint_summary_dat$clust_with_anno))

timepoint_summary_dat <- timepoint_summary_dat %>%
  filter(clust_with_anno != "8 NB")

source("workflow/paper_figures/util.high_contrast_palette.R")
color_palette <- high_contrast_palette()[1:length(levels(timepoint_summary_dat$clust_with_anno))]
names(color_palette) <- levels(timepoint_summary_dat$clust_with_anno)

timepoint_summary_dat <- timepoint_summary_dat %>%
  filter(Timepoint != "v2D5")

pdf("plots/paper_figures/double_positive_barplot.pdf", height =4, width = 8)
p <- ggplot(timepoint_summary_dat, 
       aes(x = Timepoint, y = double_div_cd19)) +
  geom_col(aes(fill = clust_with_anno)) +
  geom_text(position = position_stack(vjust = .5), size = 2.85, 
            color = "white",
            fontface = "bold",
            #aes(label = ifelse(double_div_cd19 > .0001, as.character(clust_with_anno), ""), group = clust_with_anno)) +
            aes(label = ifelse(double_div_cd19 > .00005, as.character(clust_with_anno), ""), group = clust_with_anno)) +
  ylab("RBD+S1+ \n (fraction total CD19)") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  theme(legend.title = element_blank())
print(p)
dev.off()

