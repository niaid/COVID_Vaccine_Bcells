library(tidyverse)
library(uwot)

set.seed(1)

fSOM <- readRDS("analysis_out/FlowSOM/fsom.rds")

param_dat <- readRDS("processed_data/param_dat.rds")

fsom_list <- list(fSOM = fSOM)

cell_meta <- readRDS("processed_data/cell_meta.rds")

sample_meta <- cell_meta %>%
        select(`Code ID`, timepoint, date_samp) %>%
        distinct()


# Summarize cluster counts and percentages for modeling ---------------

#fs_meta_cluster <- fSOM$metaclustering[fSOM$FlowSOM$map$mapping[, 1]]
fs_meta_cluster <- readRDS("flowsom_clusters.rds")

rbd_s1_dat <- readRDS("processed_data/spike_rbd_dat.rds")

cell_meta <- cell_meta %>%
  bind_cols(rbd_s1_dat)

mat_fluor_name <- fSOM$FlowSOM$data
cnames <- param_dat$desc[match(colnames(mat_fluor_name), param_dat$name)] 
colnames(mat_fluor_name) <- gsub(" Fix", "", cnames)

all_ix <- seq_len(nrow(cell_meta)) %>% split(paste(cell_meta$sample_id, cell_meta$timepoint, cell_meta$date_samp))

n_keep <- min(sapply(all_ix, length))

keep_ix <- lapply(all_ix, sample, size = n_keep)

keep_ix <- unlist(keep_ix, use.names =  F)
keep_ix <- union(keep_ix, which(rbd_s1_dat$gated_rbdplus | rbd_s1_dat$gated_s1plus))

keep_feat <- c("CD20 Fix", "CD138 Fix", "CD38 Fix", "CD10 Fix", "CD11c",
               "CD19 Fix", "CD27 Fix", "CD21 Fix", "IgD Fix", "IgM Fix",
               "IgG Fix", "IgA Fix")
keep_feat <- gsub(" Fix", "", keep_feat)

umap_input_mat <- mat_fluor_name[keep_ix, keep_feat]
#transform with asinh
umap_input_mat <- asinh(umap_input_mat)
embedding <- umap(umap_input_mat)

fsom_metaclust_sub <- fs_meta_cluster[keep_ix]
cell_meta_sub <- cell_meta[keep_ix, ]


umap_df = data.frame(x = embedding[,1], y = embedding[,2], 
                     rbd_plus = cell_meta_sub$gated_rbdplus,
                     s1_plus = cell_meta_sub$gated_s1plus,
                     fs_meta_cluster = fsom_metaclust_sub
                   )

out_list <- list(umap_df = umap_df, cell_meta_sub = cell_meta_sub, keep_ix = keep_ix)

dir.create("analysis_out/umap")
saveRDS(out_list, "analysis_out/umap/umap_subsample_listB.rds")




