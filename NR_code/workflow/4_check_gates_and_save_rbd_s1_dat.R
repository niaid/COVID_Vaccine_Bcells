library(tidyverse)

gate_dat <- readRDS("processed_data/gate_indices.rds")


gated_s1plus <- rowSums(gate_dat[, grepl("S1\\+", colnames(gate_dat)) |
                                   grepl("S1 Fix\\+", colnames(gate_dat))]) > 0
gated_rbdplus <- rowSums(gate_dat[, grepl("RBD\\+", colnames(gate_dat)) | 
                           grepl("RBD Fix\\+", colnames(gate_dat))]) > 0

spike_rbd_dat <- data.frame(gated_rbdplus = gated_rbdplus, gated_s1plus = gated_s1plus)
saveRDS(spike_rbd_dat, "processed_data/spike_rbd_dat.rds")



