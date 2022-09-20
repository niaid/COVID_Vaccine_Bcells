library(flowCore)
library(readxl)
library(tidyverse)


in_path <- "processed_data/ff_list.rds"

ff_list <- readRDS(in_path)

timepoints_dat <- read_excel("metadata/COVID vaccine schedule for Nick & William revised 041021.xlsx")

param_dat_list <- lapply(ff_list, function(x) parameters(x)@data)
feat_name_list <- lapply(param_dat_list, `[[`, "name")

#There are some files that have a sample ID column that I have to remove to 
#  make it a valid flowSet and work with FlowSom
#tmp <- parameters(ff_list[[1]])@data$name
#ff_tmp <- ff_list[[1]][, tmp[1:5]]

#save expression values in one big matrix for flowSOM ----------------

feature_intersect <- Reduce(intersect, feat_name_list)

ff_list <- lapply(ff_list, function(x){
  x[, feature_intersect]
})

e_list <- lapply(ff_list, exprs) 

e_mat <- do.call(rbind, e_list)

saveRDS(e_mat, "processed_data/e_mat.rds", compress = FALSE)  

# metadata attached to each cell--------------------

subj_meta <- data.frame(sample_id = names(ff_list)) %>%
  mutate(subj = sapply(strsplit(sample_id, split = "__"), `[[`, 1)) %>%
  mutate(date_samp_char_orig = sapply(strsplit(sample_id, split = "__"), `[[`, 2))

subj_meta <- subj_meta %>%
  mutate(date_samp_char = replace(date_samp_char_orig, date_samp_char_orig == "31-Dec-2021", "31-Dec-2020"))

subj_meta <- subj_meta %>%
  separate(subj, into = c("initials", "Code ID"), sep = " ")

subj_meta <- subj_meta %>%
  mutate(date_samp = as.Date(subj_meta$date_samp_char, "%d-%B-%Y"))

keep_cols <- 
  c("ID by initials", 
  "Code ID", 
  "Gender at Birth",
  "DOB", "Age at day of vaccine",
  #"Dose 1 Vaccine Date",
  "Pre...7",
  "1 week...9",
  "10 days...11",
  "2 weeks...13",
  #"Dose 2 Vaccine Date",
  "Pre...16",
  "Day 5",
  "1 week...20",
  "2 weeks...24",
  "10 days...22",
  "Visit day post dose 2...26",
  "4 weeks")

timepoint_order <- c(
  "D1_pre",
  #"D1_0days",
  "D1_1week",
  "D1_10days",
  "D1_2week",
  "D2_pre",
  #"D2_0days", 
  "D2_5days",
  "D2_1week",
  "D2_10days",
  "D2_2weeks",
  "D2_4weeks"
)

timepoints_dat_long <- timepoints_dat %>%
  select(keep_cols) %>%
  rename(D1_1week =  `1 week...9`,
         #D1_0days = `Dose 1 Vaccine Date`,
         D1_pre = Pre...7,
         D1_10days = `10 days...11`,
         D1_2week = `2 weeks...13`,
         D2_pre = Pre...16,
         #D2_0days = `Dose 2 Vaccine Date`,
         D2_5days = `Day 5`,
         D2_1week = `1 week...20`,
         D2_10days = `10 days...22`,
         D2_2weeks = `2 weeks...24`,
         D2_2weeks2 = `Visit day post dose 2...26`,
         D2_4weeks = `4 weeks`) %>%
  gather(key = timepoint, value = date_samp, -c("ID by initials", 
                                      "Code ID", 
                                      "Gender at Birth",
                                      "DOB", "Age at day of vaccine")) %>%
  mutate(timepoint = gsub("D2_2weeks2", "D2_2weeks", timepoint)) %>%
  mutate(timepoint = factor(timepoint, levels = timepoint_order))
  
stopifnot(sum(is.na(timepoints_dat_long$timepoint)) == 0)

important_dates_dat <- timepoints_dat %>%
  select(`Code ID`,
         `Dose 1 Vaccine Date`, `Dose 2 Vaccine Date`,
         )

timepoints_dat_long <- timepoints_dat_long %>%
  left_join(important_dates_dat)

#have to account for person with baseline at -17
timepoints_dat_long <- 
  timepoints_dat_long %>%
  mutate(days_diff_D1_raw = 
           as.numeric(difftime(date_samp, `Dose 1 Vaccine Date`, units="days"))) %>%
  mutate(days_diff_D2_raw = 
           as.numeric(difftime(date_samp, `Dose 2 Vaccine Date`,units="days"))) %>%
  mutate(days_diff_D1 = 
           replace(days_diff_D1_raw, days_diff_D1_raw < -10 & timepoint == "D1_pre", 0)) %>%
  mutate(days_diff_D2 = 
           replace(days_diff_D2_raw, days_diff_D1_raw < -10 & timepoint == "D1_pre", 
                   days_diff_D2_raw[which(days_diff_D1_raw < -10 & timepoint == "D1_pre")] -
                   days_diff_D1_raw[which(days_diff_D1_raw < -10 & timepoint == "D1_pre")] ))

timepoints_dat_long %>%
  filter(days_diff_D1_raw < 0) %>% as.data.frame()

timepoints_dat_long <- timepoints_dat_long %>% filter(!is.na(date_samp))

subj_meta <- subj_meta %>%
  full_join(timepoints_dat_long)

subj_meta <- subj_meta %>% filter(!is.na(date_samp_char_orig))

subj_meta <- subj_meta[match(names(e_list), subj_meta$sample_id), ]

subj_meta <- subj_meta %>%
  mutate(study_day = days_diff_D1)

subj_meta$study_day[startsWith(as.character(subj_meta$timepoint), "D2")] <- 
  subj_meta$days_diff_D2[startsWith(as.character(subj_meta$timepoint), "D2")] + 28

saveRDS(subj_meta, "processed_data/subject_level_meta.rds")

cell_meta <- lapply(names(e_list), function(samp_id){
  n_cells <- nrow(e_list[[samp_id]])
  samp_index <- which(names(e_list) == samp_id)
  subj_meta[rep.int(samp_index, n_cells), ]
  
}) 
cell_meta <- bind_rows(cell_meta)

stopifnot(nrow(e_mat) == nrow(cell_meta))

saveRDS(cell_meta, "processed_data/cell_meta.rds")

param_dat <- param_dat_list[[1]]
saveRDS(param_dat, "processed_data/param_dat.rds")
