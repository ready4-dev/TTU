## ----echo = TRUE, message=FALSE-----------------------------------------------
library(magrittr)
library(TTU)

## -----------------------------------------------------------------------------
authors_tb <- ready4show::authors_tb

## ----authorsds, eval = knitr::is_html_output(), echo=F, results='asis'--------
authors_tb %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Authors table",
                          mkdn_tbl_ref_1L_chr = "tab:authorsds")

## -----------------------------------------------------------------------------
institutes_tb <- ready4show::institutes_tb

## ----institutesds, eval = knitr::is_html_output(), echo = F,results='asis'----
institutes_tb %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Author institutions table",
                          mkdn_tbl_ref_1L_chr = "tab:institutesds")

## -----------------------------------------------------------------------------
header_yaml_args_ls <- make_header_yaml_args_ls(authors_tb = authors_tb,
                                                     institutes_tb = institutes_tb,
                                                     title_1L_chr = "A hypothetical study using fake data for instructional purposes only",
                                                     keywords_chr = c("this","is","a","replication","using","fake","data","do", "not","cite"))

## -----------------------------------------------------------------------------
output_format_ls <- make_output_format_ls(manuscript_outp_1L_chr = "Word",
                                               manuscript_digits_1L_int = 2L,
                                               supplementary_outp_1L_chr = "PDF",
                                               supplementary_digits_1L_int = 2L)

## ----message=FALSE------------------------------------------------------------
ds_tb <- readRDS(url("https://dataverse.harvard.edu/api/access/datafile/4750597")) 

## ----inputds, echo=F, eval = knitr::is_html_output(), results='asis'----------
ds_tb %>%
    head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Input dataset",
                          mkdn_tbl_ref_1L_chr = "tab:inputds")

## -----------------------------------------------------------------------------
use_fake_data_1L_lgl <-  TRUE

## -----------------------------------------------------------------------------
dictionary_tb <- make_eq5d_ds_dict(ds_tb) 

## ----dictionary, echo=F, eval = knitr::is_html_output(), results='asis'-------
dictionary_tb %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Data dictionary",
                          mkdn_tbl_ref_1L_chr = "tab:dictionary")

## ----echo = TRUE--------------------------------------------------------------
ds_descvs_ls <- make_ds_descvs_ls(candidate_predrs_chr = c("k10_int","psych_well_int"),#
                     cohort_descv_var_nms_chr = c("d_age", "Gender", "d_relation_s","d_sexual_ori_s", "Region", "CALD", "d_studying_working"),
                     dictionary_tb = dictionary_tb, 
                     id_var_nm_1L_chr = "uid", 
                     msrmnt_date_var_nm_1L_chr = "data_collection_dtm",
                     round_var_nm_1L_chr = "Timepoint", 
                     round_vals_chr = c("BL", "FUP"),
                     maui_item_pfx_1L_chr = "eq5dq_", 
                     utl_wtd_var_nm_1L_chr = "EQ5D_total_dbl", 
                     utl_unwtd_var_nm_1L_chr = "EQ5d_cumulative_dbl")

## -----------------------------------------------------------------------------
predictors_lup <- TTU_predictors_lup(make_pt_TTU_predictors_lup(short_name_chr = ds_descvs_ls$candidate_predrs_chr,
                                              long_name_chr = c("Kessler Psychological Distress - 10 Item Total Score",
                                                   "Overall Psychological Wellbeing (Winefield et al. 2012)"),
                                              min_val_dbl = c(10,18),
                                              max_val_dbl = c(50,90),
                                              class_chr = "integer",
                                              increment_dbl = 1,
                                              class_fn_chr = "as.integer",
                                              mdl_scaling_dbl = 0.01,
                                              covariate_lgl = F))

## -----------------------------------------------------------------------------
maui_params_ls <- make_maui_params_ls(maui_itm_short_nms_chr = c("Mobility", "Self-care", "Activities","Pain","Anxiety"),
                                      maui_scoring_fn = function(data_tb,
                                             maui_item_pfx_1L_chr,
                                             id_var_nm_1L_chr, 
                                             utl_wtd_var_nm_1L_chr,
                                             utl_unwtd_var_nm_1L_chr){
                    require(eq5d)
                    data_tb %>% 
                      dplyr::rename_with(~stringr::str_replace(.x,maui_item_pfx_1L_chr,""),
                                         dplyr::starts_with(maui_item_pfx_1L_chr)) %>%
                      dplyr::mutate(`:=`(!!rlang::sym(utl_wtd_var_nm_1L_chr), 
                                         eq5d::eq5d(., country="UK", version = "5L", type = "CW"))) %>%
                      dplyr::rename_with(~paste0(maui_item_pfx_1L_chr,.x),c("MO","SC","UA","PD","AD")) %>%
                      dplyr::mutate(`:=`(!!rlang::sym(utl_unwtd_var_nm_1L_chr),                                          rowSums(dplyr::across(dplyr::starts_with(maui_item_pfx_1L_chr))))) %>%
                      dplyr::filter(!is.na(!!rlang::sym(utl_unwtd_var_nm_1L_chr)))
                    },
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
                                                             use_fake_data_1L_lgl = use_fake_data_1L_lgl,
                                                             control_ls = list(adapt_delta = 0.99))

## ---- warning=F, message=F, eval=F--------------------------------------------
#  path_params_ls <- make_path_params_ls(use_fake_data_1L_lgl = use_fake_data_1L_lgl,
#                                        write_new_dir_1L_lgl = T)

## ----eval = F-----------------------------------------------------------------
#  valid_params_ls_ls <- make_valid_params_ls_ls(analysis_core_params_ls,
#                                                     candidate_covar_nms_chr = c("d_sex_birth_s", "d_age",  "d_sexual_ori_s", "d_relation_s", "d_studying_working"),
#                                                     ds_tb = ds_tb,
#                                                     maui_params_ls = maui_params_ls,
#                                                     path_params_ls = path_params_ls,
#                                                     prefd_mdl_types_chr = c("OLS_NTF", "GLM_GSN_LOG", "BET_CLL"))

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  write_report(params_ls = valid_params_ls_ls$params_ls,
#                    paths_ls = path_params_ls$paths_ls,
#                    rprt_nm_1L_chr = "AAA_PMRY_ANLYS_MTH",
#                    abstract_args_ls = NULL,
#                    header_yaml_args_ls = header_yaml_args_ls)

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  rprt_lups_ls <- write_mdl_smry_rprt(header_yaml_args_ls,
#                           path_params_ls = path_params_ls,
#                           use_fake_data_1L_lgl = use_fake_data_1L_lgl,
#                           output_format_ls = output_format_ls,
#                           use_shareable_mdls_1L_lgl = T)

## ----echo = TRUE--------------------------------------------------------------
dv_ds_nm_and_url_chr <- c("fakes", # REPLACE WITH DATAVERSE IN WHICH YOUR DATASET IS LOCATED
                          "https://doi.org/10.7910/DVN/612HDC") # REPLACE WITH YOUR DATASET'S DOI
                                 

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  if(!is.null(dv_ds_nm_and_url_chr)){
#    write_study_outp_ds(dv_ds_nm_and_url_chr = dv_ds_nm_and_url_chr,
#                        rprt_lups_ls = rprt_lups_ls,
#                        output_format_ls = output_format_ls,
#                        path_params_ls = path_params_ls,
#                        use_fake_data_1L_lgl = use_fake_data_1L_lgl)
#  }

## ----echo = FALSE, eval=FALSE-------------------------------------------------
#  # Render Scientific Summary
#  # if(!use_fake_data_1L_lgl){
#    # var_nm_change_lup <- tibble::tibble(old_nms_chr = c("k10","psychwell"),
#    #                                     new_nms_chr = c("K10", "Winefield wellbeing"))
#  study_descs_ls <- make_study_descs_ls(health_utl_nm_1L_chr = "EQ-5D",
#                                             time_btwn_bl_and_fup_1L_chr = "three months",
#                                             predr_ctgs_ls = list(`psychological distress` = c("k10"),
#                                                                  `psychological wellbeing` =  c("psychwell")))
#  spine_of_results_ls <- make_results_ls_spine(output_data_dir_1L_chr = path_params_ls$paths_ls$output_data_dir_1L_chr,
#                                                    #var_nm_change_lup = var_nm_change_lup,
#                                                    study_descs_ls = study_descs_ls,
#                                                    nbr_of_digits_1L_int = output_format_ls$manuscript_digits_1L_int)
#  cs_ts_ratios_tb = tibble::tibble(predr_nm_chr = c("k10","psychwell"),
#                                   ratios_chr = c(paste0("about ",
#                                                         round(spine_of_results_ls$mdl_coef_ratios_ls$`psychological distress`,2)),
#                                                  paste0(round(spine_of_results_ls$mdl_coef_ratios_ls$`psychological wellbeing`,2))
#                                   ))
#  # ctgl_vars_regrouping_ls = list('Primary Diagnosis' = list(anxdpr = list(name_1L_chr = "anxiety/depression",
#  #                                                                           ctgs_chr = c("Anxiety", "Depression", "Depression and Anxiety"))),
#  #                                  'Clinical Stage' = list(early = list(name_1L_chr = "early (prior to first episode of a serious mental disorder) clinical stages",
#  #                                                                       ctgs_chr = c("0-1a","1b"))))
#  sig_covars_some_predrs_mdls_tb = tibble::tibble(covars_ls = list(c("No other confounding factor")
#                                                                   # ,c("sex at birth"),c("primary diagnosis", "clinical staging","age")
#                                                                   ),
#                                                                   predrs_mdls_ls = list(c("either predictor")#, c("K6"),c("anxiety and depression measurements other than PHQ-9")
#                                                                                         ),
#                                                                   sig_thresh_dbl = c(NA_real_#,0.01,NA_real_
#                                                                                      ),
#                                                                   assn_strgth_chr = c(NA_character_#,NA_character_,"weakly"
#                                                                                       ))
#  results_ls <- make_results_ls(spine_of_results_ls,
#                                cs_ts_ratios_tb = cs_ts_ratios_tb,
#                                ctgl_vars_regrouping_ls = NULL,
#                                sig_covars_some_predrs_mdls_tb = sig_covars_some_predrs_mdls_tb,
#                                sig_thresh_covars_1L_chr = "0.005")
#  saveRDS(results_ls,paste0(path_params_ls$paths_ls$output_data_dir_1L_chr,"/results_ls.RDS"))
#    rmarkdown::render("AQoL_VIH/AQoL_VIH.Rmd",
#                      output_format = NULL,
#                      params = list(figures_in_body_lgl = F,
#                                    output_type_1L_chr = "Word",
#                                    results_ls = results_ls,
#                                    tables_in_body_lgl = F),
#                      output_file = paste0("AQoL_VIH", ".docx"),
#                      output_dir = path_params_ls$paths_ls$reports_dir_1L_chr)
#  # }

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  write_to_delete_ds_copies(path_params_ls$paths_ls)

