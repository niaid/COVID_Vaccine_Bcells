fsom <- readRDS("~/Documents/covid_bcell_vax/20210505/analysis_out/FlowSOM/fsom.rds")

library(FlowSOM)
clusts <- GetMetaclusters(fsom)
saveRDS(clusts, "~/Documents/covid_bcell_vax/20210505/workflow_for_upload/flowsom_clusters.rds")
