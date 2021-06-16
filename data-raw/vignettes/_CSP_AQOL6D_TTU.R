## ----echo = TRUE--------------------------------------------------------------
library(magrittr)

## ----message=FALSE------------------------------------------------------------
ds_tb <- youthvars::replication_popl_tb %>% 
            youthvars::transform_raw_ds_for_analysis() 

## ----inputds, eval = knitr::is_html_output(), results='asis'------------------
ds_tb %>%
    head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Input dataset",
                          mkdn_tbl_ref_1L_chr = "tab:inputds")

## -----------------------------------------------------------------------------
use_fake_data_1L_lgl <-  TRUE

## -----------------------------------------------------------------------------
dictionary_tb <-  youthvars::make_final_rpln_ds_dict() #youthvars::make_tfd_repln_ds_dict_r3()

## ----dictionary, eval = knitr::is_html_output(), results='asis'---------------
dictionary_tb %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Data dictionary",
                          mkdn_tbl_ref_1L_chr = "tab:dictionary")

## ----echo = TRUE--------------------------------------------------------------
ds_descvs_ls <- TTU::make_ds_descvs_ls(candidate_predrs_chr = c("K6","PHQ9"),
                     cohort_descv_var_nms_chr = c("d_age", "d_relation_s","d_studying_working","c_p_diag_s","c_clinical_staging_s"),
                     dictionary_tb = dictionary_tb, 
                     id_var_nm_1L_chr = "fkClientID", 
                     msrmnt_date_var_nm_1L_chr = "d_interview_date",
                     round_var_nm_1L_chr = "round", 
                     round_vals_chr = c("Baseline", "Follow-up"),
                     maui_item_pfx_1L_chr = "aqol6d_q", 
                     utl_wtd_var_nm_1L_chr = "aqol6d_total_w", 
                     utl_unwtd_var_nm_1L_chr = "aqol6d_total_c")

## -----------------------------------------------------------------------------
predictors_lup <- TTU::TTU_predictors_lup(TTU::make_pt_TTU_predictors_lup(short_name_chr = ds_descvs_ls$candidate_predrs_chr,
                                              long_name_chr = c("K6 total score", "PHQ9 total score"),
                                              min_val_dbl = 0,
                                              max_val_dbl = c(24,27),
                                              class_chr = "integer",
                                              increment_dbl = 1,
                                              class_fn_chr = "as.integer",
                                              mdl_scaling_dbl = 0.01,
                                              covariate_lgl = F))

## -----------------------------------------------------------------------------
maui_params_ls <- TTU::make_maui_params_ls(maui_itm_short_nms_chr = c("Household tasks", "Getting around","Morbility","Self care","Enjoy close rels","Family rels", "Community involvement","Despair","Worry", "Sad", "Agitated","Energy level","Control", "Coping","Frequency of pain", "Degree of pain","Pain interference","Vision", "Hearing","Communication"),
                                      maui_scoring_fn = youthvars::add_adol6d_scores)

## -----------------------------------------------------------------------------
mdl_types_lup <- TTU::get_cndts_for_mxd_mdls()

## ----mdltypeslup, eval = knitr::is_html_output(), results='asis'--------------
mdl_types_lup %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Candidate model types lookup table",
                          mkdn_tbl_ref_1L_chr = "tab:mdltypeslup")

## -----------------------------------------------------------------------------
mdl_smry_ls <- TTU::make_mdl_smry_ls(mdl_types_lup = mdl_types_lup,
                                     folds_1L_int = 10L,
                                     max_nbr_of_boruta_mdl_runs_int = 300L)

