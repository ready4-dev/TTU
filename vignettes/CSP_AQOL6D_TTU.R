## ----echo = TRUE, message=FALSE-----------------------------------------------
library(magrittr)
library(TTU)

## -----------------------------------------------------------------------------
header_yaml_args_ls <- make_header_yaml_args_ls(authors_tb = ready4show::authors_tb,
                                                institutes_tb = ready4show::institutes_tb,
                                                title_1L_chr = "A hypothetical study using fake data for instructional purposes only",
                                                keywords_chr = c("this","is","a","replication","using","fake","data","do", "not","cite"))

## -----------------------------------------------------------------------------
output_format_ls <- make_output_format_ls(manuscript_outp_1L_chr = "Word",
                                               manuscript_digits_1L_int = 2L,
                                               supplementary_outp_1L_chr = "PDF",
                                               supplementary_digits_1L_int = 2L)

## ----message=FALSE------------------------------------------------------------
ds_tb <- youthvars::replication_popl_tb %>% 
            youthvars::transform_raw_ds_for_analysis() 

## -----------------------------------------------------------------------------
use_fake_data_1L_lgl <-  TRUE

## ----echo = TRUE--------------------------------------------------------------
ds_descvs_ls <- make_ds_descvs_ls(candidate_predrs_chr = c("K6","PHQ9"),#
                     cohort_descv_var_nms_chr = c("d_age", "d_relation_s","d_studying_working","c_p_diag_s","c_clinical_staging_s"),
                     dictionary_tb = youthvars::make_final_rpln_ds_dict(), 
                     id_var_nm_1L_chr = "fkClientID", 
                     msrmnt_date_var_nm_1L_chr = "d_interview_date",
                     round_var_nm_1L_chr = "round", 
                     round_vals_chr = c("Baseline", "Follow-up"),
                     maui_item_pfx_1L_chr = "aqol6d_q", 
                     utl_wtd_var_nm_1L_chr = "aqol6d_total_w", 
                     utl_unwtd_var_nm_1L_chr = "aqol6d_total_c")

## -----------------------------------------------------------------------------
predictors_lup <- TTU_predictors_lup(make_pt_TTU_predictors_lup(short_name_chr = ds_descvs_ls$candidate_predrs_chr,
                                              long_name_chr = c("K6 total score", "PHQ9 total score"),
                                              min_val_dbl = 0,
                                              max_val_dbl = c(24,27),
                                              class_chr = "integer",
                                              increment_dbl = 1,
                                              class_fn_chr = "as.integer",
                                              mdl_scaling_dbl = 0.01,
                                              covariate_lgl = F))

## -----------------------------------------------------------------------------
maui_params_ls <- make_maui_params_ls(maui_itm_short_nms_chr = c("Household tasks", "Getting around","Morbility","Self care","Enjoy close rels","Family rels", "Community involvement","Despair","Worry", "Sad", "Agitated","Energy level","Control", "Coping","Frequency of pain", "Degree of pain","Pain interference","Vision", "Hearing","Communication"),
                                      maui_scoring_fn = youthvars::add_adol6d_scores,
                    utl_min_val_1L_dbl = 0.03)

## -----------------------------------------------------------------------------
mdl_smry_ls <- make_mdl_smry_ls(mdl_types_lup = get_cndts_for_mxd_mdls(),
                                     folds_1L_int = 10L,
                                     max_nbr_of_boruta_mdl_runs_int = 300L)

## ----echo = TRUE--------------------------------------------------------------
analysis_core_params_ls <- make_analysis_core_params_ls(ds_descvs_ls = ds_descvs_ls,
                                                             mdl_smry_ls = mdl_smry_ls,
                                                             output_format_ls = output_format_ls,
                                                             predictors_lup = predictors_lup,
                                                             use_fake_data_1L_lgl = use_fake_data_1L_lgl)

## ---- warning=F, message=F, eval=F--------------------------------------------
#  path_params_ls <- make_path_params_ls(use_fake_data_1L_lgl = use_fake_data_1L_lgl,
#                                        write_new_dir_1L_lgl = T)

## ----eval = F-----------------------------------------------------------------
#  valid_params_ls_ls <- make_valid_params_ls_ls(analysis_core_params_ls,
#                                                     candidate_covar_nms_chr = c("d_age", "SOFAS", "c_p_diag_s", "c_clinical_staging_s", "d_relation_s", "d_studying_working"),
#                                                     ds_tb = ds_tb,
#                                                     maui_params_ls = maui_params_ls,
#                                                     path_params_ls = path_params_ls,
#                                                     prefd_mdl_types_chr = c("GLM_GSN_LOG",
#                                                                             "OLS_CLL"))

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  write_report(params_ls = valid_params_ls_ls$params_ls,
#                    paths_ls = path_params_ls$paths_ls,
#                    rprt_nm_1L_chr = "Main_Analysis_Rprt",
#                    abstract_args_ls = NULL,
#                    header_yaml_args_ls = header_yaml_args_ls)

