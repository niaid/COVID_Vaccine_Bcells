
library(data.table)

e_mat <- fread("./all_subjects_cd19_positive_and_keys/all_subjects_cd19_positive.csv")

saveRDS(e_mat, "processed_data/e_mat.rds")
