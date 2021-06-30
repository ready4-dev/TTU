## ----message=FALSE, warning=FALSE---------------------------------------------
library(TTU)

## -----------------------------------------------------------------------------
seed_1L_int = 12345
set.seed(seed_1L_int)

## ----warning = F--------------------------------------------------------------
path_to_write_to_1L_chr <- "Output" 
dir.create(path_to_write_to_1L_chr)

## ----message=FALSE------------------------------------------------------------
data_tb <- readRDS(url("https://dataverse.harvard.edu/api/access/datafile/4750597"))

## ----inputds, echo=F, eval = knitr::is_html_output(), results='asis'----------
data_tb %>%
    head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Input dataset",
                          mkdn_tbl_ref_1L_chr = "tab:inputds")

## -----------------------------------------------------------------------------
predictors_lup <- TTU_predictors_lup(make_pt_TTU_predictors_lup(short_name_chr = c("k10_int","psych_well_int"),
                                              long_name_chr = c("Kessler Psychological Distress - 10 Item Total Score",
                                                   "Overall Wellbeing Measure (Winefield et al. 2012)"),
                                              min_val_dbl = c(10,18),
                                              max_val_dbl = c(50,90),
                                              class_chr = "integer",
                                              increment_dbl = 1,
                                              class_fn_chr = "as.integer",
                                              mdl_scaling_dbl = 0.01,
                                              covariate_lgl = F))

## -----------------------------------------------------------------------------
dictionary_tb <- make_eq5d_ds_dict(data_tb,
                                   predictors_lup = predictors_lup)

## ----message=FALSE, results='hide', warning=FALSE-----------------------------
library(eq5d)
data_tb <- data_tb %>% 
        dplyr::rename_with(~stringr::str_replace(.x,"eq5dq_",""),dplyr::starts_with("eq5dq_")) %>%
  dplyr::mutate(`:=`(EQ5D_total_dbl, 
                     eq5d::eq5d(., country="UK", version = "5L", type = "CW"))) %>%
    dplyr::rename_with(~paste0("eq5dq_",.x),c("MO","SC","UA","PD","AD"))

## -----------------------------------------------------------------------------
data_tb <- data_tb %>%
  ready4use::add_labels_from_dictionary(dictionary_tb = dictionary_tb,
                                        remove_old_lbls_1L_lgl = T) 

## ----scoredds, eval = knitr::is_html_output(), echo=F,results='asis'----------
data_tb %>%
    head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Scored and labelled dataset",
                          use_lbls_as_col_nms_1L_lgl = T,
                          mkdn_tbl_ref_1L_chr = "tab:scoredds")

## -----------------------------------------------------------------------------
ds_smry_ls <- make_ds_smry_ls(candidate_predrs_chr = predictors_lup$short_name_chr,
                              candidate_covar_nms_chr = c("d_sex_birth_s", "d_age",  "d_sexual_ori_s", "d_relation_s", "d_studying_working"),
                              depnt_var_nm_1L_chr = "EQ5D_total_dbl",
                              dictionary_tb = dictionary_tb,
                              id_var_nm_1L_chr = "uid",
                              round_var_nm_1L_chr = "Timepoint",
                              round_bl_val_1L_chr = "BL",
                              predictors_lup = predictors_lup)

## -----------------------------------------------------------------------------
mdl_types_lup <- TTU::get_cndts_for_mxd_mdls()
mdl_types_lup %>%
    ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Model types lookup table")

## -----------------------------------------------------------------------------
mdl_smry_ls <- mdl_types_lup %>%
  make_mdl_smry_ls()

## ----warning=F, message=F-----------------------------------------------------
cmprsns_ls <- write_mdl_cmprsn(scored_data_tb = data_tb,
                                ds_smry_ls = ds_smry_ls,
                                mdl_smry_ls = mdl_smry_ls,
                                output_data_dir_1L_chr = path_to_write_to_1L_chr,
                                seed_1L_int = seed_1L_int)

## ----mdl_cmprsn, eval = knitr::is_html_output(), results='asis'---------------
cmprsns_ls$mdl_smry_ls$smry_of_sngl_predr_mdls_tb %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Comparison of model types for highest correlated predictor using baseline data",
                          mkdn_tbl_ref_1L_chr = "tab:mdl_cmprsn")

## -----------------------------------------------------------------------------
cmprsns_ls$mdl_smry_ls$prefd_mdl_types_chr

## -----------------------------------------------------------------------------
cmprsns_ls$mdl_smry_ls$prefd_mdl_types_chr <- c("OLS_NTF", "GLM_GSN_LOG", "BET_CLL")

## ----message=FALSE, results='hide', warning=FALSE-----------------------------
cmprsns_ls <- write_predr_and_covars_cmprsn(scored_data_tb = data_tb,
                                            bl_tb = cmprsns_ls$bl_tb,
                                            ds_smry_ls = cmprsns_ls$ds_smry_ls,
                                            mdl_smry_ls = cmprsns_ls$mdl_smry_ls,
                                            output_data_dir_1L_chr = path_to_write_to_1L_chr,
                                            seed_1L_int = seed_1L_int)

## -----------------------------------------------------------------------------
cmprsns_ls$mdl_smry_ls$predr_cmprsn_tb %>%
    ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Comparison of all candidate predictors using preferred model")

## -----------------------------------------------------------------------------
cmprsns_ls$mdl_smry_ls$smry_of_mdl_sngl_predrs_tb %>%
    ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Preferred single predictor model performance by candidate predictor")

## -----------------------------------------------------------------------------
cmprsns_ls$mdl_smry_ls$mdls_with_covars_smry_tb %>%
    ready4show::print_table(output_type_1L_chr = "HTML",
                            caption_1L_chr = "Preferred performance of models with covariate predictors by candidate predictor", 
                          use_lbls_as_col_nms_1L_lgl = T)

## -----------------------------------------------------------------------------
cmprsns_ls$mdl_smry_ls$signt_covars_chr

## -----------------------------------------------------------------------------
# cmprsns_ls$mdl_smry_ls$prefd_covars_chr <- c("COVARIATE WE WANT TO USE", "ANOTHER COVARIATE")

## ----message=FALSE, results='hide', warning=FALSE, fig.show ='hide'-----------
outp_smry_ls <- write_mdls_with_covars_cmprsn(scored_data_tb = data_tb,
                                           bl_tb = cmprsns_ls$bl_tb,
                                           ds_smry_ls = cmprsns_ls$ds_smry_ls,
                                           mdl_smry_ls = cmprsns_ls$mdl_smry_ls,
                                           output_data_dir_1L_chr = path_to_write_to_1L_chr,
                                           seed_1L_int = seed_1L_int)

## ----eval = F-----------------------------------------------------------------
#  outp_smry_ls <- write_ts_mdls_from_alg_outp(outp_smry_ls,
#                                              predictors_lup = predictors_lup ,
#                                              utl_min_val_1L_dbl = 0.03,
#                                              backend_1L_chr = "cmdstanr",
#                                              new_dir_nm_1L_chr = "F_TS_Mdls",
#                                              iters_1L_int = 4000L,
#                                              prior_ls = NULL,
#                                              control_ls = list(adapt_delta = 0.99))

## ----eval=F-------------------------------------------------------------------
#  write_to_delete_mdl_fls(outp_smry_ls)
#  outp_smry_ls$scored_data_tb <- NULL

## ----eval=F-------------------------------------------------------------------
#  saveRDS(outp_smry_ls,paste0(outp_smry_ls$path_to_write_to_1L_chr,"/I_ALL_OUTPUT_.RDS"))

