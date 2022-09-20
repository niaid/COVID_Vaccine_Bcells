
#These have already been run and you will use the output but are there for your reference
#workflow/NOT_RUN_1_read_in_data2_wsp_files.R
#workflow/NOT_RUN_2_save_expression_mat_and_cellmeta.R

## download the data and convert to rds ----------------------

Rscript workflow/0_prep_zenodo_data.R

# These will generate the outputs necessary for the figures --

Rscript workflow/3_run_flowsom.R

Rscript workflow/4_check_gates_and_save_rbd_s1_dat.R

Rscript workflow/5_save_spike_rbd_summ.R

Rscript workflow/6_umap_subsampleB.R

Rscript workflow/7_anova_models.R

# making the figures ----------------------------------------

Rscript workflow/paper_figures/antibody_plots/ab_cor_plots.R
Rscript workflow/paper_figures/antibody_plots/ab_d2d28_density.R
Rscript workflow/paper_figures/cluster_timecourse_plots.R
Rscript workflow/paper_figures/double_pos_stacked_barplot.R
Rscript workflow/paper_figures/plot_anova.R
Rscript workflow/paper_figures/plot_marker_mfis_ggplot.R
Rscript workflow/paper_figures/plot_umap_subsampleB_v2.R
Rscript workflow/paper_figures/save_supp_table.R
