#Tried to copy scaling and what not from
#https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.23030
#Comparison of clustering methods for high‐dimensional single‐cell flow and mass cytometry data Lukas M. Weber  Mark D. Robinson

#from the paper
#The clustering algorithms were run on all remaining single, live cells; no additional pre‐gating was performed, since our aim is to evaluate performance in maximally automated settings. In addition, we did not perform any standardization of individual protein marker dimensions. This was unnecessary since the arcsinh already transforms all dimensions to comparable scales; and importantly, standardization of dimensions that do not contain a true signal could amplify the effect of noise and outliers, adversely affecting clustering performance.

library(flowCore)
library(dplyr)
library(FlowSOM)
library(CytoML)
library(flowWorkspace)

wsp_files <- "~/Documents/covid_bcell_vax/20210505/102521 Belgian not pregnant females_FCS & FlowJo files/21-Oct-2021 Belgian nonP females CoV2.wsp"

ws <- open_flowjo_xml(wsp_files)

gs <- flowjo_to_gatingset(ws, transform = FALSE, name = 1, skip_faulty_gate = TRUE)

#save expression values -----------------------------------------------------------------

#get CD19 cytoframe
ff_list <- lapply(gs, function(x){
  gh_pop_get_data(x, "CD19+")
})

#makes a flowFrame which is the same as a cytoframe, but stored in memory instead of with pointer
ff_list_inmem <- lapply(ff_list, cytoframe_to_flowFrame)
names(ff_list_inmem)

#save flowFrame list to be used in flowSOM
out_path <- "processed_data/ff_list_out_path.rds"
dir.create(dirname(out_path))
saveRDS(ff_list_inmem, out_path)

#examples of how to access data from flowFrames ---------------------

#how to see the markers and names
param_dat_list <- lapply(ff_list_inmem, function(x) parameters(x)@data)
str(param_dat_list[[1]])

#how to get the matrix of expression values
expression_matrix <- exprs(ff_list_inmem[[1]])
str(expression_matrix)

#get gates ---------------------------------------------------------
nodelists <- lapply(gs, function(x){
  gs_get_pop_paths(x, path = "auto")
})

shared_gates <- Reduce(intersect, nodelists)

nodelists_sub <- lapply(nodelists, function(x){
  x[x %in% shared_gates]
})

nodelist_mat <- do.call(rbind, nodelists_sub)


gate_indices <- lapply(seq_along(nodelists), function(i){
  nodel <- nodelists[[i]]
  print(i)
  
  cd19_plus <- flowWorkspace::gh_pop_get_indices(gs[[i]][[1]], "CD19+")
  
  dat <- lapply(nodel, function(node){
    flowWorkspace::gh_pop_get_indices(gs[[i]][[1]], node)[cd19_plus]
  }) %>% as.data.frame()
  colnames(dat) <- nodel
  dat
})
names(gate_indices) <- names(nodelists)

gate_indices_combined <- bind_rows(gate_indices, .id = "sample_id")

gates_out_path <- "/processed_data/gates_list.rds"
saveRDS(gate_indices, gates_out_path)

gates_combined_out_path <- "processed_data/gates_dataframe.rds"
saveRDS(gate_indices_combined, gates_combined_out_path)
