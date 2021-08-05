## ----echo = TRUE--------------------------------------------------------------
library(magrittr)

## ----message=FALSE------------------------------------------------------------
ds_tb <- readRDS(url("https://dataverse.harvard.edu/api/access/datafile/4750597")) 

## ----inputds, eval = knitr::is_html_output(), results='asis'------------------
ds_tb %>%
    head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Input dataset",
                          mkdn_tbl_ref_1L_chr = "tab:inputds")

## -----------------------------------------------------------------------------
use_fake_data_1L_lgl <-  TRUE

## -----------------------------------------------------------------------------
authors_tb <- ready4show::authors_tb
institutes_tb <- ready4show::institutes_tb

## -----------------------------------------------------------------------------
header_yaml_args_ls = list(authors_tb = authors_tb,
                           institutes_tb = institutes_tb,
                           fl_nm_1L_chr = "header_common.yaml",
                           title_1L_chr = "A hypothetical study using fake data for instructional purposes only",
                           keywords_chr = c("this","is","a","replication","using","fake","data","do", "not","cite"))
nbr_of_digits_1L_int = 2L
output_types_ls = list(manuscript_1L_chr = "Word",
                       supplementary_1L_chr = "PDF")

## -----------------------------------------------------------------------------
path_params_ls <- list(path_from_top_level_1L_chr = normalizePath("../") %>% strsplit("\\\\") %>% purrr::pluck(1) %>% tail(1),
                  path_to_data_from_top_level_chr = "fake_data.rds", 
                  path_to_current_1L_chr = normalizePath(".") %>% strsplit("\\\\") %>% purrr::pluck(1) %>% tail(1))

## ----echo = TRUE--------------------------------------------------------------
ds_descvs_ls <- list(candidate_predrs_chr = c("K10_int","Psych_well_int"),
                     cohort_descv_var_nms_chr = c("d_age", "Gender", "d_relation_s","d_sexual_ori_s", "Region", "CALD", "d_studying_working"),
                     dictionary_tb = TTU::make_eq5d_ds_dict(ds_tb), 
                     id_var_nm_1L_chr = "uid", 
                     msrmnt_date_var_nm_1L_chr = "data_collection_dtm",
                     round_var_nm_1L_chr = "Timepoint", 
                     round_vals_chr = c("BL", "FUP"),
                     maui_item_pfx_1L_chr = "eq5dq_", 
                     utl_wtd_var_nm_1L_chr = "EQ5D_total_dbl", 
                     utl_unwtd_var_nm_1L_chr = "EQ5d_cumulative_dbl")

## -----------------------------------------------------------------------------
predictors_lup <- TTU::TTU_predictors_lup(TTU::make_pt_TTU_predictors_lup(short_name_chr = ds_descvs_ls$candidate_predrs_chr,
                                              long_name_chr = c("Kessler Psychological Distress - 10 Item Total Score",
                                                   "Overall Wellbeing Measure (Winefield et al. 2012)"),
                                              min_val_dbl = c(10,18),
                                              max_val_dbl = c(50,90),
                                              class_chr = "integer",
                                              increment_dbl = 1,
                                              class_fn_chr = "as.integer",
                                              mdl_scaling_dbl = 0.01,
                                              covariate_lgl = F))

## ----echo = TRUE--------------------------------------------------------------
analysis_params_ls <- list(ds_descvs_ls = ds_descvs_ls,
                           iters_1L_int = 4000L,
                           mdl_smry_ls = list(mdl_types_lup = TTU::mdl_types_lup,
                                              mdl_types_chr = TTU::mdl_types_lup$short_name_chr,
                                              choose_from_pfx_chr = c("GLM","OLS","BET"),
                                              folds_1L_int = 10L,
                                              max_nbr_of_boruta_mdl_runs_int = 300L),
                           nbr_of_digits_1L_int = nbr_of_digits_1L_int,
                           output_type_1L_chr = output_types_ls$supplementary_1L_chr,
                           predictors_lup = predictors_lup,
                           seed_1L_int = 12345,
                           use_fake_data_1L_lgl =  use_fake_data_1L_lgl,
                           prior_ls = NULL, 
                           control_ls = list(adapt_delta = 0.99))

## ----echo = TRUE--------------------------------------------------------------
rprt_lup <- TTU::rprt_lup 
file.create("aaaaaaaaaa.txt") 
paths_ls <- TTU::write_main_oupt_dir(path_params_ls, 
                                     use_fake_data_1L_lgl = use_fake_data_1L_lgl)
paths_ls$path_to_current_1L_chr <- path_params_ls$path_to_current_1L_chr

## ----echo = TRUE--------------------------------------------------------------
dv_ls <- list(primary_dv_chr = c("fakes","https://doi.org/10.7910/DVN/612HDC"))

