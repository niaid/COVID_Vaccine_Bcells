library(tidyverse)


cell_freq <- readRDS("analysis_out/FlowSOM/cluster_freq_and_antigen_specificity_summary.rds")

anova_res <- readRDS("analysis_out/anova_results_clusters.rds")

cell_freq <- cell_freq %>%
  filter(!is.na(date_samp) & fs_meta_cluster != "8") %>%
  select(`Code ID`, timepoint, fs_meta_cluster, n_cells_in_cluster, cluster_fraction_total,n_cd19_total, double_count, double_fraction) %>%
  rename(Cluster = fs_meta_cluster,
    cluster_fraction_total_CD10 = cluster_fraction_total, n_total_CD19 = n_cd19_total, `n_RBD+S1+` = double_count, 
         `fraction_RBD+S1+` = double_fraction)
levels(cell_freq$Cluster) <- paste0("C", levels(cell_freq$Cluster))


dir.create("supp_tables/")
write_csv(cell_freq, "supp_tables/cluster_frequencies_and_spike_specificity.csv")

anova_res <- anova_res %>%
  rename(Cluster = fs_meta_cluster) %>%
  select(-clust_with_anno)

levels(anova_res$Cluster) <- paste0("C", levels(anova_res$Cluster))

anova_res2 <- anova_res %>%
  filter(Cluster != "8") %>%
  filter(feature %in% c("cluster_fraction_total", "double_fraction")) %>%
  mutate(p.adjusted2 = p.adjust(p.value, method = "fdr")) %>%
  group_by(feature) %>%
  mutate(p.adjusted3 = p.adjust(p.value, method = "fdr"))

anova_res2 %>%
  select(-c(p.adjusted, p.adjusted2)) %>%
  mutate(feature = gsub("total", "total_CD19+", feature)) %>%
  mutate(feature = gsub("double", "RBD+S1+", feature)) %>%
  rename(p.adjusted = p.adjusted3) %>%
  write_csv("supp_tables/anova_results.csv")
