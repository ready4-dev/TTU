## -----------------------------------------------------------------------------
library(TTU)


## -----------------------------------------------------------------------------
path_to_write_to_1L_chr <- "Output" 
dir.create(path_to_write_to_1L_chr)


## -----------------------------------------------------------------------------
path_to_inp_data_1L_chr <- "Data" 
dir.create(path_to_inp_data_1L_chr)
write.csv(replication_popl_tb %>% dplyr::select(-BADS,-GAD7, -OASIS,-SCARED),"Data/data.csv", row.names=FALSE, col.names = T)
write.csv(repln_ds_dict_r3,"Data/dictionary.csv", row.names=FALSE)


## ----message = F--------------------------------------------------------------
data_tb <- readr::read_csv("Data/data.csv") 
dictionary_r3 <- readr::read_csv("Data/dictionary.csv") %>% ready4use::ready4_dictionary()


## ----message=FALSE, results='hide', warning=FALSE-----------------------------
scored_data_tb <- data_tb %>% add_adol6d_scores() %>%
  ready4use::add_labels_from_dictionary(dictionary_r3) 


## -----------------------------------------------------------------------------
ds_smry_ls <- list(candidate_predrs_chr = c("K6","PHQ9"),
                   candidate_covar_nms_chr = c("d_sex_birth_s", "SOFAS", "c_p_diag_s", "c_clinical_staging_s", "d_age",  "d_sexual_ori_s", "d_country_bir_s", "d_relation_s", "d_studying_working"),
                   dep_var_nm_1L_chr = "aqol6d_total_w",
                   id_var_nm_1L_chr = "fkClientID",
                   round_var_nm_1L_chr = "round",
                   round_bl_val_1L_chr = "Baseline")


## -----------------------------------------------------------------------------
ds_smry_ls$predictors_lup <- TTU_predictors_lup(make_pt_TTU_predictors_lup(short_name_chr = ds_smry_ls$candidate_predrs_chr,
                                              long_name_chr = c("K6 total score", "PHQ9 total score"),
                                              min_val_dbl = 0,
                                              max_val_dbl = c(24,27),
                                              class_chr = "integer",
                                              increment_dbl = 1,
                                              class_fn_chr = "integer",
                                              mdl_scaling_dbl = 0.01,
                                              covariate_lgl = F))


## -----------------------------------------------------------------------------
data("mdl_types_lup",package = "TTU") # LUP needs to be a class of TTU
mdl_types_lup


## -----------------------------------------------------------------------------
mdl_smry_ls <- list(mdl_types_lup = mdl_types_lup,
                    mdl_types_chr = mdl_types_lup$short_name_chr,
                    choose_from_pfx_chr = c("GLM","OLS","BET"))


## -----------------------------------------------------------------------------
seed_1L_int = 12345
set.seed(seed_1L_int)


## -----------------------------------------------------------------------------
bl_tb <- transform_ds_for_tstng(scored_data_tb, 
                                dep_var_nm_1L_chr = ds_smry_ls$dep_var_nm_1L_chr,
                                candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr,
                                dep_var_max_val_1L_dbl = 0.999,
                                round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr, 
                                round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr)
bl_tb 


## -----------------------------------------------------------------------------
ds_smry_ls$candidate_predrs_chr <- reorder_cndt_predrs_chr(ds_smry_ls$candidate_predrs_chr,
        data_tb = bl_tb, dep_var_nm_1L_chr = ds_smry_ls$dep_var_nm_1L_chr)


## -----------------------------------------------------------------------------
mdl_smry_ls$predr_var_nm_1L_chr <- ds_smry_ls$candidate_predrs_chr[1]
mdl_smry_ls$predr_var_desc_1L_chr <- ds_smry_ls$predictors_lup %>% ready4fun::get_from_lup_obj(match_value_xx = mdl_smry_ls$predr_var_nm_1L_chr, match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = "long_name_chr",
        evaluate_lgl = F)
mdl_smry_ls$predr_var_desc_1L_chr


## -----------------------------------------------------------------------------
mdl_smry_ls$predr_vals_dbl <- make_predr_vals(mdl_smry_ls$predr_var_nm_1L_chr, candidate_predrs_lup = ds_smry_ls$predictors_lup)
mdl_smry_ls$predr_vals_dbl


## -----------------------------------------------------------------------------
mdl_smry_ls$n_folds_1L_int = 10L
mdl_smry_ls$smry_of_sngl_predr_mdls_tb <- write_sngl_predr_multi_mdls_outps(data_tb = bl_tb,
        n_folds_1L_int = mdl_smry_ls$n_folds_1L_int, mdl_types_chr = mdl_smry_ls$mdl_types_chr,
        dep_var_nm_1L_chr = ds_smry_ls$dep_var_nm_1L_chr, predr_var_nm_1L_chr = mdl_smry_ls$predr_var_nm_1L_chr,
        predr_var_desc_1L_chr = mdl_smry_ls$predr_var_nm_1L_chr, predr_vals_dbl = mdl_smry_ls$predr_vals_dbl,
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, new_dir_nm_1L_chr =  "A_Candidate_Mdls_Cmprsns",
        mdl_types_lup = mdl_smry_ls$mdl_types_lup) # Takes more than vignette time limit to render. Need to add and link to report summarising output.
mdl_smry_ls$smry_of_sngl_predr_mdls_tb # Make this a class with a print method


## -----------------------------------------------------------------------------
mdl_smry_ls$prefd_mdl_types_chr <- make_prefd_mdls_vec(mdl_smry_ls$smry_of_sngl_predr_mdls_tb,
                                           choose_from_pfx_chr = mdl_smry_ls$choose_from_pfx_chr)
mdl_smry_ls$prefd_mdl_types_chr


## -----------------------------------------------------------------------------
mdl_smry_ls$prefd_mdl_types_chr <- c("GLM_GSN_LOG","OLS_CLL")


## -----------------------------------------------------------------------------
mdl_smry_ls$max_nbr_of_boruta_mdl_runs_int = 300L
mdl_smry_ls$predr_cmprsns_tb <- write_predr_cmprsn_outps(data_tb = bl_tb,
                                             path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                             new_dir_nm_1L_chr = "B_Candidate_Predrs_Cmprsns",
                                             dep_var_nm_1L_chr = ds_smry_ls$dep_var_nm_1L_chr, 
                                             candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr,
                                             max_nbr_of_boruta_mdl_runs_int = mdl_smry_ls$max_nbr_of_boruta_mdl_runs_int)
mdl_smry_ls$predr_cmprsns_tb # Add report and class with print method



## -----------------------------------------------------------------------------
    mdl_smry_ls$smry_of_mdl_sngl_predrs_tb <- write_mdl_type_multi_outps(data_tb = bl_tb,
        n_folds_1L_int = mdl_smry_ls$n_folds_1L_int, predrs_var_nms_chr = mdl_smry_ls$predr_cmprsns_tb$predr_chr,
        mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1], dep_var_nm_1L_chr = ds_smry_ls$dep_var_nm_1L_chr,
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, new_dir_nm_1L_chr = "C_Predrs_Sngl_Mdl_Cmprsns",
        fl_nm_pfx_1L_chr = "C_PREDR", mdl_types_lup = mdl_smry_ls$mdl_types_lup) 
mdl_smry_ls$smry_of_mdl_sngl_predrs_tb # Add report and class with print method


## -----------------------------------------------------------------------------
bl_tb <- scored_data_tb %>% transform_ds_for_tstng(candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr,
        covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr, remove_all_mssng_1L_lgl = T,
        round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr)
mdl_smry_ls$mdls_with_covars_smry_tb <- write_mdl_type_covars_mdls(bl_tb,
        dep_var_nm_1L_chr = ds_smry_ls$dep_var_nm_1L_chr, predrs_var_nms_chr = ds_smry_ls$candidate_predrs_chr,
        covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr, mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1],
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, new_dir_nm_1L_chr = "D_Predr_Covars_Cmprsn",
        fl_nm_pfx_1L_chr = "D_CT",
        mdl_types_lup = mdl_smry_ls$mdl_types_lup)
mdl_smry_ls$mdls_with_covars_smry_tb # Add report and class with print method. Note negative values


## -----------------------------------------------------------------------------
mdl_smry_ls$signt_covars_chr <- get_signft_covars(mdls_with_covars_smry_tb = mdl_smry_ls$mdls_with_covars_smry_tb,
                                      covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr) # THIS FN NEEDS FIXING FOR CAT VARS
mdl_smry_ls$signt_covars_chr


## -----------------------------------------------------------------------------
mdl_smry_ls$prefd_covars_chr <- NA_character_
    if (is.na(mdl_smry_ls$prefd_covars_chr))
        mdl_smry_ls$prefd_covars_chr <- mdl_smry_ls$signt_covars_chr


## -----------------------------------------------------------------------------
    empty_tb <- write_mdl_type_multi_outps(data_tb = bl_tb, n_folds_1L_int = NULL,
        start_1L_chr = NA_character_, predrs_var_nms_chr = mdl_smry_ls$predr_cmprsns_tb$predr_chr,
        covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr, mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[2],
        dep_var_nm_1L_chr = ds_smry_ls$dep_var_nm_1L_chr, path_to_write_to_1L_chr = path_to_write_to_1L_chr,
        new_dir_nm_1L_chr = "E_Predrs_W_Covars_Sngl_Mdl_Cmprsns", mdl_types_lup = mdl_smry_ls$mdl_types_lup, fl_nm_pfx_1L_chr = "E_CK_CV")
    mdl_smry_ls$predr_vars_nms_ls <- make_predr_vars_nms_ls(main_predrs_chr = mdl_smry_ls$predr_cmprsns_tb$predr_chr,
        covars_ls = list(mdl_smry_ls$prefd_covars_chr))


## -----------------------------------------------------------------------------
mdl_smry_ls$mdl_nms_ls <- make_mdl_nms_ls(mdl_smry_ls$predr_vars_nms_ls, mdl_types_chr = mdl_smry_ls$prefd_mdl_types_chr)


## -----------------------------------------------------------------------------
outp_smry_ls <- list(scored_data_tb = scored_data_tb, 
                     smry_of_sngl_predr_mdls_tb = mdl_smry_ls$smry_of_sngl_predr_mdls_tb,
                     prefd_mdl_types_chr = mdl_smry_ls$prefd_mdl_types_chr, 
                     predr_cmprsns_tb = mdl_smry_ls$predr_cmprsns_tb,
                     smry_of_mdl_sngl_predrs_tb = mdl_smry_ls$smry_of_mdl_sngl_predrs_tb,
                     mdls_with_covars_smry_tb = mdl_smry_ls$mdls_with_covars_smry_tb,
                     signt_covars_chr = mdl_smry_ls$signt_covars_chr, 
                     prefd_covars_chr = mdl_smry_ls$prefd_covars_chr,
                     dep_var_nm_1L_chr = ds_smry_ls$dep_var_nm_1L_chr, 
                     predr_vars_nms_ls = mdl_smry_ls$predr_vars_nms_ls,
                     mdl_nms_ls = mdl_smry_ls$mdl_nms_ls, 
                     id_var_nm_1L_chr = ds_smry_ls$id_var_nm_1L_chr,
                     round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr,
                     round_bl_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr,
                     path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
                     seed_1L_int = seed_1L_int,
                     n_folds_1L_int = mdl_smry_ls$n_folds_1L_int, 
                     max_nbr_of_boruta_mdl_runs_int = mdl_smry_ls$max_nbr_of_boruta_mdl_runs_int,
                     mdl_types_lup = mdl_smry_ls$mdl_types_lup, 
                     file_paths_chr = list.files(path_to_write_to_1L_chr, recursive = T))


## -----------------------------------------------------------------------------
outp_smry_ls <- write_ts_mdls_from_alg_outp(outp_smry_ls,
                                            fn_ls = list(fit_gsn_log_lnk,fit_clg_log_tfmn),
                                            backend_1L_chr = "cmdstanr",
                                            iters_1L_int = 500L)


## ----eval=FALSE---------------------------------------------------------------
## outp_smry_ls$dv_ls <- list(dv_nm_1L_chr = "firstbounce",
##                                ds_url_1L_chr = "https://doi.org/10.7910/DVN/JC6PTV",
##                                parent_dv_dir_1L_chr = "../../../../Data/Dataverse")


## ----eval=T-------------------------------------------------------------------
outp_smry_ls <- write_shareable_mdls(outp_smry_ls,
                                     new_dir_nm_1L_chr = "G_Shareable",
                                     shareable_title_detail_1L_chr = "The model that can be used to predict adolescent AQoL6D. Note this model is a placeholder as it has been estimated from synthetic data.")


## ----eval=T-------------------------------------------------------------------
outp_smry_ls$fk_data_tb <- TTU::make_fake_ts_data(outp_smry_ls)


## -----------------------------------------------------------------------------
saveRDS(outp_smry_ls,paste0(outp_smry_ls$path_to_write_to_1L_chr,"/I_ALL_OUTPUT_.RDS"))

