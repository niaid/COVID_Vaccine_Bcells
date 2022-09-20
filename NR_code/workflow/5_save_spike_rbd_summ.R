library(tidyverse)


param_dat <- readRDS("processed_data/param_dat.rds")

cell_meta <- readRDS("processed_data/cell_meta.rds")

spike_rbd_dat <- readRDS("processed_data/spike_rbd_dat.rds")

sample_meta <- cell_meta %>%
        select(`Code ID`, timepoint, date_samp) %>%
        distinct()

# Summarize cluster counts and percentages for modeling ---------------

# If using flowSOM clusters generated from 3_run_flowsom.R
# fSOM <- readRDS("analysis_out/FlowSOM/fsom.rds")
# fs_meta_cluster <- fSOM$metaclustering[fSOM$FlowSOM$map$mapping[, 1]]
# If using the pre-clustered data
fs_meta_cluster <- readRDS("flowsom_clusters.rds")


cell_meta <- cell_meta %>%
  mutate(S1 = spike_rbd_dat$gated_s1plus, RBD = spike_rbd_dat$gated_rbdplus) %>%
  mutate(fs_meta_cluster = fs_meta_cluster)

summ <- cell_meta %>%
  group_by(fs_meta_cluster, `Code ID`, timepoint, date_samp, .drop = FALSE) %>% 
  summarize(s1_count = sum(S1), rbd_count = sum(RBD), n_cells_in_cluster = n(),
            double_count = sum(S1 & RBD)) %>%
  mutate(s1_fraction = s1_count /n_cells_in_cluster, rbd_fraction = rbd_count / n_cells_in_cluster,
         double_fraction = double_count/n_cells_in_cluster) %>% 
  ungroup() %>%
  group_by(`Code ID`, timepoint, date_samp) %>%
  mutate(n_cd19_total = sum(n_cells_in_cluster)) %>%
  mutate(cluster_fraction_total = n_cells_in_cluster / n_cd19_total) %>%
  mutate(s1_fraction_total_cd19 = s1_count / n_cd19_total, 
         rbd_fraction_total_cd19 = rbd_count / n_cd19_total,
         double_fraction_total_cd19 = double_count / n_cd19_total) %>%
  ungroup()

saveRDS(summ, "analysis_out/FlowSOM/cluster_freq_and_antigen_specificity_summary.rds")
