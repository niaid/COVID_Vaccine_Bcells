library(tidyverse)
library(readxl)

summ <- readRDS("analysis_out/FlowSOM/cluster_freq_and_antigen_specificity_summary.rds")

sample_meta <- readRDS("processed_data/subject_level_meta.rds")

clust_anno <- read_excel("Clusters v050621 Annotated 051321.xlsx", n_max = 31)

summ <- left_join(summ, sample_meta)

feats <- c("s1_fraction", "rbd_fraction", "double_fraction", 
           "s1_fraction_total_cd19", "rbd_fraction_total_cd19", 
           "double_fraction_total_cd19", "cluster_fraction_total")
keep_cols <- c("fs_meta_cluster", "Code ID", "timepoint", "date_samp", "study_day", feats)

summ_long <- summ %>% select(keep_cols) %>%
  gather(key = feature, value = value, 
         -c(fs_meta_cluster, `Code ID`, timepoint, date_samp, study_day)) %>% 
  rename(CodeID = `Code ID`) %>%
  filter(!is.na(study_day)) %>%
  mutate(CodeID = factor(CodeID))

do_anova <- function(dat){
  full_mod <- lm(value ~ timepoint + CodeID, data = dat)
  res <- anova(full_mod) %>% broom::tidy()
  res
}


mod_results <- summ_long %>%
  group_by(fs_meta_cluster, feature) %>%
  do(do_anova(.)) %>% 
  ungroup()

mod_results <- mod_results %>%
  mutate(p.adjusted = p.adjust(p.value, method = "fdr"))

clust_anno <- clust_anno %>%
  rename(fs_meta_cluster = Cluster) %>%
  mutate(fs_meta_cluster = factor(fs_meta_cluster, levels = levels(mod_results$fs_meta_cluster)))

plot_dat <- mod_results %>%
  mutate(clust_with_anno = fs_meta_cluster)

pop_prefix1 <- clust_anno$Population[match(levels(plot_dat$clust_with_anno), as.character(clust_anno$fs_meta_cluster))]
levels(plot_dat$clust_with_anno) <- paste(levels(plot_dat$clust_with_anno), pop_prefix1)

saveRDS(plot_dat, file = "analysis_out/anova_results_clusters.rds")

