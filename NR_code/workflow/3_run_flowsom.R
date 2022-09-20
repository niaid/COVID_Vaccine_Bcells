library(FlowSOM)
library(flowCore)

e_mat <- readRDS("processed_data/e_mat.rds")
param_dat <- readRDS("processed_data/param_dat.rds")
#cell_meta <- readRDS("processed_data/cell_meta.rds")

keep_feat <- c("CD20 Fix", "CD138 Fix", "CD38 Fix", "CD10 Fix", "CD11c",
               "CD19 Fix", "CD27 Fix", "CD21 Fix", "IgD Fix", "IgM Fix",
               "IgG Fix", "IgA Fix")

# have serum antibodies binding the surface receptors
keep_feat_fluor_name <- param_dat$name[match(keep_feat, param_dat$desc)]
names(keep_feat_fluor_name) <- keep_feat

#e_mat <- asinh(e_mat/150)

nmarkers <- length(keep_feat_fluor_name)


ff <- flowFrame(exprs = e_mat)
set.seed(3345)
fSOM <- FlowSOM(ff, colsToUse = keep_feat_fluor_name, nClus = 30, transform = FALSE,
                scale = FALSE)

dir.create("analysis_out/FlowSOM", recursive = TRUE)

saveRDS(fSOM, "analysis_out/FlowSOM/fsom.rds", compress = FALSE)

