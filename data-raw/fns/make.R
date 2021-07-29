make_abstract_args_ls <- function(results_ls,
                                  fl_nm_1L_chr = "abstract.txt"){
  mdl_cmprsns_ls <- get_mdl_cmprsns(results_ls, as_list_1L_lgl = T)
  abstract_args_ls <- list(abstract_ls = list(Background = get_background_text(results_ls),
                                              Objectives = paste0("We aimed to identify the best transfer to utility (TTU) regression models to predict ",
                                                                  get_hlth_utl_nm(results_ls, short_nm_1L_lgl = F),
                                                                  " (",
                                                                  get_hlth_utl_nm(results_ls),
                                                                  ") utility and evaluate the predictive ability of ",
                                                                  get_nbr_of_predrs(results_ls),
                                                                  " candidate measure",
                                                                  ifelse(get_nbr_of_predrs(results_ls, as_words_1L_lgl = F)>1,"s",""),
                                                                  " of ",
                                                                  get_predr_ctgs(results_ls),
                                                                  "."),
                                              Methods = paste0(results_ls$study_descs_ls$sample_desc_1L_chr,
                                                               " Follow-up measurements were ",
                                                               results_ls$study_descs_ls$time_btwn_bl_and_fup_1L_chr,
                                                               " after baseline. ",
                                                               paste0(length(mdl_cmprsns_ls$OLS) %>%
                                                                        xfun::numbers_to_words() %>%
                                                                        Hmisc::capitalize()),
                                                               " Ordinary Least Squares (OLS) and ",
                                                               length(mdl_cmprsns_ls$GLM) %>%
                                                                 xfun::numbers_to_words(),
                                                               " generalised linear models (GLMs) were explored to identify the best TTU algorithm. ",
                                                               " Predictive ability of ",
                                                               get_nbr_of_predrs(results_ls),
                                                               " candidate measure",
                                                               ifelse(get_nbr_of_predrs(results_ls, as_words_1L_lgl = F)>1,"s",""),
                                                               " of ",
                                                               get_predr_ctgs(results_ls),
                                                               " were assessed using ten fold cross validation",
                                                               ifelse(get_nbr_of_predrs(results_ls, as_words_1L_lgl = F)>1," and forest models",""),
                                                               ". Linear / generalised linear mixed effect models were then used to construct longitudinal predictive models for ",
                                                               get_hlth_utl_nm(results_ls),
                                                               " change."),
                                              Results = paste0(make_ten_fold_text(results_ls, for_abstract_1L_lgl = T),
                                                               make_random_forest_text(results_ls, for_abstract_1L_lgl = T),
                                                               ". ",
                                                               make_selected_mdl_text(results_ls, for_abstract_1L_lgl = T),
                                                               paste0(" The mean ratio between the within-person and between-person associated coefficients was ",
                                                                      make_within_between_ratios_text(results_ls),
                                                                      ".")),
                                              Conclusions = get_conclusion_text(results_ls),
                                              Data = make_data_availability_text(results_ls)),
                           fl_nm_1L_chr = fl_nm_1L_chr)

 return(abstract_args_ls)
}
make_all_mdl_types_smry_tbl <- function(outp_smry_ls,
                                        mdls_tb){
  mdls_ls <- make_mdls_ls(outp_smry_ls,
                          mdls_tb = mdls_tb)
  all_mdl_types_smry_tbl_tb <- 1:length(mdls_ls) %>% # Modified from 1:2 - need to test.
    purrr::map_dfc(
      ~ {
        make_mdl_type_smry_tbl(mdls_tb = mdls_tb,
                               mdl_nms_chr = mdls_ls[[.x]],
                               mdl_type_1L_chr = outp_smry_ls$prefd_mdl_types_chr[.x],
                               add_mdl_nm_sfx_1L_lgl = T)
      }
    ) %>%
    dplyr::select(-paste0("Parameter_",outp_smry_ls$prefd_mdl_types_chr[-1])) %>%
    dplyr::rename_with(~"Parameter",.cols = paste0("Parameter_",outp_smry_ls$prefd_mdl_types_chr[1]))
  return(all_mdl_types_smry_tbl_tb)
}
make_analysis_core_params_ls <- function(ds_descvs_ls,
                                         mdl_smry_ls = make_mdl_smry_ls(),
                                         output_format_ls = make_output_format_ls(),
                                         predictors_lup,
                                         control_ls = NULL,
                                         iters_1L_int = 4000L,
                                         prefd_covars_chr = NULL,
                                         prefd_mdl_types_chr = NULL,
                                         prior_ls = NULL,
                                         seed_1L_int = 12345,
                                         candidate_covar_nms_chr = NULL,#NA_character_,
                                         use_fake_data_1L_lgl = NULL # F
                                         ){
  if(missing(candidate_covar_nms_chr)){
    candidate_covar_nms_chr <- ds_descvs_ls$candidate_covar_nms_chr
  }else{
    warning("candidate_covar_nms_chr is soft deprecated - it is recommended to specify the candiate covariates as part of the list object passed to the ds_descvs_ls argument")
  }
  if(missing(use_fake_data_1L_lgl)){
    use_fake_data_1L_lgl <- ds_descvs_ls$is_fake_1L_lgl
  }else{
    warning("use_fake_data_1L_lgl is soft deprecated - it is recommended to specify whether the dataset is fake in the is_fake_1L_lgl element of the list obkect passed to the ds_descvs_ls argument")
  }
  analysis_core_params_ls <- list(candidate_covar_nms_chr = candidate_covar_nms_chr,
                                  ds_descvs_ls = ds_descvs_ls,
                                  iters_1L_int = iters_1L_int,
                                  mdl_smry_ls = mdl_smry_ls,
                                  nbr_of_digits_1L_int = output_format_ls$supplementary_digits_1L_int,
                                  output_type_1L_chr = output_format_ls$supplementary_outp_1L_chr,
                                  predictors_lup = predictors_lup,
                                  prefd_covars_chr = prefd_covars_chr,
                                  prefd_mdl_types_chr = prefd_mdl_types_chr,
                                  seed_1L_int = seed_1L_int,
                                  use_fake_data_1L_lgl = use_fake_data_1L_lgl,
                                  prior_ls = prior_ls,
                                  control_ls = control_ls)
  return(analysis_core_params_ls)
}
make_analysis_ds_smry_ls <- function(ds_descvs_ls,
                                     candidate_covar_nms_chr,
                                     predictors_lup){
  analysis_ds_smry_ls <- list(candidate_predrs_chr = ds_descvs_ls$candidate_predrs_chr,
                              candidate_covar_nms_chr = candidate_covar_nms_chr,
                              depnt_var_nm_1L_chr = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                              id_var_nm_1L_chr = ds_descvs_ls$id_var_nm_1L_chr,
                              predictors_lup = predictors_lup,
                              round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                              round_bl_val_1L_chr = ds_descvs_ls$round_vals_chr[1],
                              dictionary_tb = ds_descvs_ls$dictionary_tb)

  return(analysis_ds_smry_ls)
}
make_bl_fup_add_to_row_ls <-  function(df,
                                       n_at_bl_1L_int,
                                       n_at_fup_1L_int){
  add_to_row_ls <- list(pos = list(-1, nrow(df)),
                        command = c(paste("\\toprule \n",
                                          paste0("\\multicolumn{2}{c}{} & \\multicolumn{2}{c}{\\textbf{Baseline (N=",n_at_bl_1L_int,")}} & \\multicolumn{2}{c}{\\textbf{Follow-up (N=",n_at_fup_1L_int,")}} \\\\\n")),
                                    paste("\\bottomrule \n"
                                    )
                        )
  )
  return(add_to_row_ls)
}
make_brms_mdl_print_ls <- function (mdl_ls, label_stub_1L_chr, caption_1L_chr, output_type_1L_chr = "PDF",
                                    digits_1L_dbl = 2, big_mark_1L_chr = " ")
{
    smry_mdl_ls <- summary(mdl_ls, digits = 4)
    mdl_smry_chr <- smry_mdl_ls %>% utils::capture.output()
    idx_dbl <- c("Formula: ", "Samples: ", "Group-Level Effects: ",
                 "Population-Level Effects: ", "Family Specific Parameters: ",
                 "Samples were drawn using ") %>% purrr::map_dbl(~mdl_smry_chr %>%
                                                                     startsWith(.x) %>% which())
    data_tb <- make_brms_mdl_smry_tbl(smry_mdl_ls, grp_1L_chr = mdl_smry_chr[idx_dbl[3]],
                                      popl_1L_chr = mdl_smry_chr[idx_dbl[4]], fam_1L_chr = mdl_smry_chr[idx_dbl[5]])
    bold_lgl <- data_tb$Parameter %in% c(mdl_smry_chr[idx_dbl[3]],
                                         mdl_smry_chr[idx_dbl[4]], mdl_smry_chr[idx_dbl[5]])
    if (output_type_1L_chr == "PDF") {
        data_tb <- data_tb %>% dplyr::mutate(Parameter = purrr::map_chr(Parameter,
                                                                        ~.x %>% Hmisc::latexTranslate()))
    }
    data_tb <- data_tb %>% dplyr::mutate(Parameter = Parameter %>%
                                             purrr::map2_chr(dplyr::all_of(bold_lgl), ~ifelse(.y &
                                                                                                  output_type_1L_chr == "PDF", paste0("\\textbf{",
                                                                                                                                      .x, "}"), .x)))
    if (output_type_1L_chr != "PDF") {
        data_tb <- data_tb %>% dplyr::mutate(dplyr::across(c(Bulk_ESS,
                                                             Tail_ESS), ~format(., big.mark = big_mark_1L_chr)))
    }
    if (output_type_1L_chr == "HTML") {
        data_tb <- data_tb %>% dplyr::mutate(dplyr::across(where(is.numeric),
                                                           ~format(round(., digits = digits_1L_dbl), digits = digits_1L_dbl,
                                                                   nsmall = digits_1L_dbl)))
    }
    data_tb <- data_tb %>% dplyr::mutate(dplyr::across(where(is.character),
                                                       ~dplyr::case_when(is.na(.) ~ "", . == "NA" ~ "", endsWith(.,
                                                                                                                 " NA") ~ "", TRUE ~ .)))
    end_matter_1L_chr <- trimws(mdl_smry_chr[idx_dbl[6]:length(mdl_smry_chr)]) %>%
        paste0(collapse = " ")
    brms_mdl_print_ls <- list(part_1 = mdl_smry_chr[idx_dbl[1]],
                              part_2 = "\n\n", part_3 = c(trimws(mdl_smry_chr[1:(idx_dbl[2] -
                                                                                     1)][-idx_dbl[1]]), paste0(trimws(mdl_smry_chr[idx_dbl[2]]),
                                                                                                               " ", trimws(mdl_smry_chr[idx_dbl[2] + 1]), collapse = " ")) %>%
                                  paste0(collapse = ifelse(output_type_1L_chr == "PDF",
                                                           "\n\n", "\n")), part_4 = "\n\n", part_5 = list(data_tb = data_tb,
                                                                                                          output_type_1L_chr = output_type_1L_chr, caption_1L_chr = caption_1L_chr,
                                                                                                          mkdn_tbl_ref_1L_chr = paste0("tab:", label_stub_1L_chr), merge_row_idx_int = as.integer(which(bold_lgl)),
                                                                                                          digits_dbl = c(ifelse(output_type_1L_chr == "PDF",
                                                                                                                                0, NA_real_) %>% purrr::discard(is.na), names(data_tb) %>%
                                                                                                                             purrr::map_dbl(~ifelse(.x %in% c("Bulk_ESS",
                                                                                                                                                              "Tail_ESS"), 0, digits_1L_dbl))), big_mark_1L_chr = big_mark_1L_chr,
                                                                                                          hline_after_ls = c(-1, 0), sanitize_fn = force, footnotes_chr = NA_character_),
                              part_6 = end_matter_1L_chr)
    if (output_type_1L_chr != "PDF") {
        brms_mdl_print_ls$part_5$footnotes_chr <- c(paste0(brms_mdl_print_ls$part_1,
                                                           ifelse(output_type_1L_chr == "Word", "", "\n")),
                                                    brms_mdl_print_ls$part_3, brms_mdl_print_ls$part_6)
        brms_mdl_print_ls$part_6 <- NULL
    }
    else {
        footnotes_chr <- c(mdl_smry_chr[idx_dbl[1]], trimws(mdl_smry_chr[1:(idx_dbl[2] -
                                                                                1)][-idx_dbl[1]]), trimws(mdl_smry_chr[idx_dbl[2]]),
                           trimws(mdl_smry_chr[idx_dbl[2] + 1]), trimws(mdl_smry_chr[idx_dbl[6]:length(mdl_smry_chr)])) %>%
            Hmisc::latexTranslate()
        footnotes_chr[1] <- footnotes_chr[1] %>% stringr::str_replace("~",
                                                                      "\\\\textasciitilde")
        brms_mdl_print_ls$part_5$add_to_row_ls <- list(pos = purrr::map(c(0,
                                                                     rep(nrow(data_tb), length(footnotes_chr) + 1)), ~.x),
                                                  command = c(names(data_tb) %>% Hmisc::latexTranslate() %>%
                                                                  paste0(collapse = " & ") %>% paste0("\\\\\n"),
                                                              c("\\toprule\n"), footnotes_chr %>% purrr::map_chr(~paste0("\\multicolumn{",
                                                                                                                         ncol(data_tb), "}{l}{", paste0("{\\footnotesize ",
                                                                                                                                                        .x, "}\n", collapse = ","), "}\\\\\n"))))
    }
    return(brms_mdl_print_ls)
}
make_brms_mdl_smry_tbl <- function (smry_mdl_ls, grp_1L_chr, popl_1L_chr, fam_1L_chr)
{
    brms_mdl_smry_tb <- purrr::map(1:length(smry_mdl_ls$random),
                                   ~make_mdl_smry_elmt_tbl(ctg_chr = c(ifelse(.x == 1, grp_1L_chr,
                                                                              character(0)),
                                                                       paste0(names(smry_mdl_ls$ngrps)[.x],
                                                                              " (Number of levels: ",
                                                                              smry_mdl_ls$ngrps[.x][[1]],
                                                                              ")")),
                                                           mat = smry_mdl_ls$random[.x][[1]])) %>%
      dplyr::bind_rows(make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$fixed,
                                              ctg_chr = popl_1L_chr),
                       make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$spec_pars,
                                              ctg_chr = fam_1L_chr))
    return(brms_mdl_smry_tb)
}
make_cmpst_sctr_and_dnsty_plt <- function(outp_smry_ls,
                                          output_data_dir_1L_chr,
                                          predr_var_nms_chr,
                                          labels_chr = c("A","B","C","D"),
                                          label_x_1L_dbl = 0.1,
                                          label_y_1L_dbl = 0.9,
                                          label_size_1L_dbl = 22){
  filtered_paths_chr <- outp_smry_ls$file_paths_chr %>% purrr::discard(~endsWith(.x,"_sim_sctr.png")|endsWith(.x,"_sim_dnst.png")|endsWith(.x,"_cnstrd_sctr_plt.png")|endsWith(.x,"_cnstrd_dnst.png"))
  plot_ls <- paste0(output_data_dir_1L_chr,"/",filtered_paths_chr[filtered_paths_chr %>% purrr::map_lgl(~stringr::str_detect(.x,paste0(predr_var_nms_chr,"_1")) & (stringr::str_detect(.x,"_dnst.png") | stringr::str_detect(.x,"_sctr_plt.png")))])  %>% purrr::map(~cowplot::ggdraw() + cowplot::draw_image(.x))
  composite_plt <- cowplot::plot_grid(plot_ls[[1]],plot_ls[[2]],plot_ls[[3]],plot_ls[[4]],nrow = 2, labels = labels_chr, label_x = label_x_1L_dbl,label_y = label_y_1L_dbl, label_size = label_size_1L_dbl)
}
make_cndt_predr_text <- function(results_ls,
                                 type_1L_chr = "description"){
  nbr_of_predrs_1L_int <- get_nbr_of_predrs(results_ls,
                                            as_words_1L_lgl = F)
  if(type_1L_chr == "description"){
    text_1L_chr <- paste0(get_nbr_of_predrs(results_ls) %>% Hmisc::capitalize(),
                          " measure",
                          ifelse(nbr_of_predrs_1L_int>1,
                                 "s",
                                 ""),
                          " of ",
                          get_nbr_of_predrs_by_ctg(results_ls),
                          ifelse(nbr_of_predrs_1L_int>1," were"," was"),
                          " used as ",
                          ifelse(nbr_of_predrs_1L_int>1,
                                 "",
                                 "a "),
                          "candidate predictor",
                          ifelse(nbr_of_predrs_1L_int>1,
                                 "s",
                                 ""),
                          " to construct TTU models.")
  }
  if(type_1L_chr == "comparison"){
    text_1L_chr <- paste0(ifelse(nbr_of_predrs_1L_int > 1,
                                 paste0("We compared the usefulness of the candidate predictors by using a random forest model including ",
                                        ifelse(nbr_of_predrs_1L_int > 2,
                                               paste0("all ",
                                                      get_nbr_of_predrs(results_ls)),
                                               "both"),
                                        " candidate predictors and by evaluating the independent predictive ability of different candidate predictors using 10-fold cross-validation."),
                                 ""))
  }
  return(text_1L_chr)
}
make_cohort_ls <- function(descv_tbls_ls,
                           ctgl_vars_regrouping_ls = NULL,
                           nbr_of_digits_1L_int = 2L){
  numeric_vars_chr <- descv_tbls_ls$cohort_desc_tb %>% dplyr::filter(label ==
                                                                       "Median (Q1, Q3)") %>% dplyr::pull(variable)
  ctgl_vars_chr <- unique(descv_tbls_ls$cohort_desc_tb$variable[!descv_tbls_ls$cohort_desc_tb$variable %in% numeric_vars_chr])
  nbr_by_round_dbl <- paste0(ds_descvs_ls$round_vals_chr,"_val_1_dbl") %>% purrr::map_dbl(~descv_tbls_ls$cohort_desc_tb %>%
                                                                                         dplyr::filter(variable == ctgl_vars_chr[1]) %>%
                                                                                         dplyr::pull(.x) %>%
                                                                                         as.numeric() %>% purrr::map_dbl(~.x[[1]]) %>% sum())
  numeric_vars_smry_ls <- numeric_vars_chr %>%
    purrr::map(~{
      var_smry_tb <- descv_tbls_ls$cohort_desc_tb %>%
        dplyr::filter(variable == .x)
      list(bl_min_1L_dbl = var_smry_tb %>%
             dplyr::filter(label == "Min - Max") %>%
             dplyr::pull(!!rlang::sym(paste0(ds_descvs_ls$round_vals_chr[1],"_val_1_dbl"))) %>%
             as.numeric(),
           bl_max_1L_dbl = var_smry_tb %>%
             dplyr::filter(label == "Min - Max") %>%
             dplyr::pull(!!rlang::sym(paste0(ds_descvs_ls$round_vals_chr[1],"_val_2_ls"))) %>%
             as.numeric(),
           bl_mean_1L_dbl = round(var_smry_tb %>% dplyr::filter(label == "Mean (SD)") %>%
                                    dplyr::pull(!!rlang::sym(paste0(ds_descvs_ls$round_vals_chr[1],"_val_1_dbl"))) %>%
                                    as.numeric(),
                                  nbr_of_digits_1L_int),
           bl_sd_1L_dbl = round(var_smry_tb %>% dplyr::filter(label == "Mean (SD)") %>%
                                  dplyr::pull(!!rlang::sym(paste0(ds_descvs_ls$round_vals_chr[1],"_val_2_ls"))) %>%
                                  stringr::str_remove("\\(") %>%
                                  stringr::str_remove("\\)") %>%
                                  as.numeric(),
                                nbr_of_digits_1L_int))

    }) %>% stats::setNames(numeric_vars_chr)
  cohort_ls <- list(n_all_1l_dbl = descv_tbls_ls$ds_descvs_ls$nbr_participants_1L_int,
                    n_inc_1L_dbl = nbr_by_round_dbl[1],
                    n_fup_1L_dbl = nbr_by_round_dbl[2],
                    numeric_vars_smry_ls = numeric_vars_smry_ls)
  if(!is.null(ctgl_vars_regrouping_ls)){
    append_ls <- ctgl_vars_regrouping_ls %>%
      purrr::map2(names(ctgl_vars_regrouping_ls),
                  ~{

                    var_nm_1L_chr <- .y
                    .x %>%
                      purrr::map(~{
                        list(name_1L_chr = .x$name_1L_chr,
                             n_in_group_1L_dbl =  descv_tbls_ls$cohort_desc_tb %>%
                               dplyr::filter(variable == var_nm_1L_chr) %>%
                               dplyr::filter(label %in% .x$ctgs_chr) %>%
                               dplyr::pull(Baseline_val_1_dbl)  %>%
                               as.numeric() %>% purrr::map_dbl(~.x) %>% sum(),
                             n_msng_1L_dbl = descv_tbls_ls$cohort_desc_tb %>%
                               dplyr::filter(variable == var_nm_1L_chr) %>%
                               dplyr::filter(label == "Missing") %>%
                               dplyr::pull(Baseline_val_1_dbl)  %>%
                               as.numeric())
                      }
                      )

                  })
    cohort_ls <- append(cohort_ls,append_ls)
  }
  return(cohort_ls)
}
make_coi_text <- function(results_ls){
  text_1L_chr <- ifelse(is.null(results_ls$study_descs_ls$coi_1L_chr),
                        "",
                        results_ls$study_descs_ls$coi_1L_chr)
  return(text_1L_chr)
}
make_correlation_text <- function(results_ls){
  correlation_text_1L_chr <- ifelse(length(results_ls$hlth_utl_and_predrs_ls$cor_seq_dscdng_chr)<2,
                                    "",
                                    paste0(results_ls$hlth_utl_and_predrs_ls$cor_seq_dscdng_chr[1],
                                           " was found to have the highest correlation with utility score both at baseline and follow-up followed by ",
                                           results_ls$hlth_utl_and_predrs_ls$cor_seq_dscdng_chr[2],
                                           ifelse(length(results_ls$hlth_utl_and_predrs_ls$cor_seq_dscdng_chr)<3,"",paste0(" and ", results_ls$hlth_utl_and_predrs_ls$cor_seq_dscdng_chr[3])),
                                           ifelse(length(results_ls$hlth_utl_and_predrs_ls$cor_seq_dscdng_chr)<4,"",paste0("; baseline and follow-up ",results_ls$hlth_utl_and_predrs_ls$cor_seq_dscdng_chr[length(results_ls$hlth_utl_and_predrs_ls$cor_seq_dscdng_chr)]," was found to have the lowest correlation coefficients with utility score")),
                                           "."))
  return(correlation_text_1L_chr)
}
make_covar_ttu_tbl_refs <- function(params_ls){
  results_ls <- params_ls$results_ls
  n_mdls_1L_int <- length(results_ls$ttu_lngl_ls$best_mdls_tb$model_type)
  n_covars_1L_int <- length(results_ls$ttu_lngl_ls$incld_covars_chr)
  text_1L_chr <- paste0(ifelse(n_covars_1L_int < 1,
                               "",
                               paste0(" (see ",
                                      ifelse(params_ls$output_type_1L_chr == "Word","","Table"),
                                      "s ",
                                      paste0("\\@ref(tab:coefscovarstype",1:n_mdls_1L_int,")", collapse = ", ") %>%
                                        stringi::stri_replace_last(fixed = ",", " and"),
                                      ")"
                               ))
  )
  return(text_1L_chr)
}
make_covariates_text <- function(results_ls){
  if(!is.null(results_ls$candidate_covars_ls)){
    if(length(results_ls$candidate_covars_ls)<1){
      text_1L_chr <- ""
    }else{
      n_predrs_1L_int <- get_nbr_of_predrs(results_ls, as_words_1L_lgl = F)
      text_1L_chr <- paste0("The confounding effect of other participant characteristics when using the candidate predictors in predicting utility score were also evaluated. Using the baseline data, ",
                            ifelse(is.na(results_ls$ttu_cs_ls$sig_covars_all_predrs_mdls_chr[1]),"no confounding factor",results_ls$ttu_cs_ls$sig_covars_all_predrs_mdls_chr %>%
                                     paste0(collapse = ", ") %>%
                                     stringi::stri_replace_last(fixed = ",", " and")),
                            " ",
                            ifelse(is.na(results_ls$ttu_cs_ls$sig_covars_all_predrs_mdls_chr[1]),"was",ifelse(length(results_ls$ttu_cs_ls$sig_covars_all_predrs_mdls_chr)==1,"was","were")),
                            " found to independently predict utility scores in models for ",
                            ifelse(n_predrs_1L_int == 1,
                                   "the ",
                                   ifelse(n_predrs_1L_int == 2,
                                          "both ",
                                          paste0("all ", get_nbr_of_predrs(results_ls)))),
                            "candidate predictor",
                            ifelse(n_predrs_1L_int == 1," ","s "),
                            "*(p<0.01)*.")
      mdls_with_signft_covars_ls <- results_ls$mdls_with_signft_covars_ls
      duplicates_int <- which(duplicated(mdls_with_signft_covars_ls))
      if(!identical(integer(0), duplicates_int)){
        unduplicated_ls <- mdls_with_signft_covars_ls[!duplicated(mdls_with_signft_covars_ls)]
        duplicated_ls <- mdls_with_signft_covars_ls[duplicates_int]
        add_to_chr <- duplicates_int %>%
          purrr::map_chr(~{
            match_chr <- mdls_with_signft_covars_ls[[.x]]
            mdls_with_signft_covars_ls[1:(.x-1)] %>%
              purrr::map_lgl(~identical(.x,match_chr)) %>%
              names()
          })
        signft_covars_chr <- names(unduplicated_ls) %>%
          purrr::map(~{
            vars_chr <- c(.x,names(duplicated_ls)[which(.x == add_to_chr)])
            paste0(vars_chr %>%
                     purrr::map_chr(~ transform_names(.x,
                                                      rename_lup = results_ls$var_nm_change_lup)) %>%
                     paste0(collapse = ", ") %>%
                     stringi::stri_replace_last(fixed = ",", " and"),
                   ifelse(length(vars_chr)>1," were significant covariates *(p<0.01)*"," was a significant covariate *(p<0.01)*"),
                   " in the ")
          }) %>%
          purrr::flatten_chr()
        mdls_ls <- unduplicated_ls
      }else{
        signft_covars_chr <- names(mdls_with_signft_covars_ls) %>%
          purrr::map_chr(~paste0(transform_names(.x,
                                                 rename_lup = results_ls$var_nm_change_lup),
                                 " was a significant covariate *(p<0.01)* in the "))
        mdls_ls <- mdls_with_signft_covars_ls
      }
      nbr_predrs_1L_int <- get_nbr_of_predrs(results_ls, as_words_1L_lgl = F)
      sig_for_some_int <- which(mdls_ls %>% purrr::map_lgl(~length(.x) != nbr_predrs_1L_int)) %>% unname()
      if(!identical(sig_for_some_int, integer(0))){
        mdl_ls <- mdl_ls[sig_for_some_int]
        signft_covars_chr <- signft_covars_chr[sig_for_some_int]
        text_1L_chr <- paste0(text_1L_chr,
                              " ",
                              mdls_ls %>% purrr::map_chr(~ .x %>%
                                                           purrr::map_chr(~ transform_names(.x,
                                                                                            rename_lup = results_ls$var_nm_change_lup)) %>%
                                                           paste0(collapse = ", ") %>%
                                                           stringi::stri_replace_last(fixed = ",", " and")) %>%
                                purrr::map2_chr(signft_covars_chr,
                                                ~ paste0(.y,
                                                         .x,
                                                         " model",
                                                         ifelse(stringr::str_detect(.x," and "),"s. ",". "))) %>%
                                paste0(collapse = ""))
      }
    }
  }else{
    text_1L_chr <- ""
  }
  return(text_1L_chr)
}
make_covar_ttu_tbl_title <- function(results_ls,
                                     ref_1L_int = 1){
  title_1L_chr <- paste0('Estimated coefficients from longitudinal TTU models based on individual candidate predictors with ',results_ls$ttu_lngl_ls$incld_covars_chr %>% paste0(collapse = ", ") %>% stringi::stri_replace_last(fixed = ",", " and"),' using ', results_ls$ttu_lngl_ls$best_mdls_tb[[ref_1L_int,"model_type"]], ' (', results_ls$ttu_lngl_ls$best_mdls_tb[[ref_1L_int,"link_and_tfmn_chr"]],')')
  return(title_1L_chr)
}
make_cs_ts_ratios_tb <- function(predr_ctgs_ls,
                                 mdl_coef_ratios_ls,
                                 nbr_of_digits_1L_int = 2L,
                                 fn_ls = NULL){
  if(is.null(fn_ls))
    fn_ls <- purrr::map(1:length(predr_ctgs_ls),
                        ~ make_mdl_coef_range_text)
  cs_ts_ratios_tb <- 1:length(predr_ctgs_ls) %>%
    purrr::map_dfr(
           ~ {
             if(length(predr_ctgs_ls %>% purrr::pluck(.x))>1){
               predr_nm_1L_chr <- paste0(names(predr_ctgs_ls)[.x] %>% tolower(),
                                         " measurements")
             }else{
               predr_nm_1L_chr <- predr_ctgs_ls %>% purrr::pluck(.x)
             }
             tibble::tibble(predr_nm_chr = predr_nm_1L_chr,
                            ratios_chr = ifelse(identical(make_mdl_coef_range_text,fn_ls %>% purrr::pluck(.x)),
                                                make_mdl_coef_range_text(mdl_coef_ratios_ls %>%
                                                                           purrr::pluck(names(predr_ctgs_ls)[.x]),
                                                                         nbr_of_digits_1L_int = nbr_of_digits_1L_int),
                                                paste0(round(rlang::exec(fn_ls %>%
                                                                    purrr::pluck(.x),
                                                                  mdl_coef_ratios_ls %>%
                                                                    purrr::pluck(names(predr_ctgs_ls)[.x])),
                                                      nbr_of_digits_1L_int),
                                                ifelse(identical(min,fn_ls %>% purrr::pluck(.x)),
                                                       " or over",
                                                       ifelse(identical(max,fn_ls %>% purrr::pluck(.x)),
                                                              " or under",
                                                              ""))))
                            )
             })
  return(cs_ts_ratios_tb)
}
make_data_availability_text <- function(results_ls){
  text_1L_chr <- ifelse(is.null(results_ls$dv_ds_nm_and_url_chr),
                        "None available",
                        paste0("Detailed results in the form of catalogues of the TTU models produced by this study and other supporting information are available in the online repository: ",results_ls$dv_ds_nm_and_url_chr[2]))
  return(text_1L_chr)
}
make_dnsty_and_sctr_plt_title <- function(results_ls){
  title_1L_chr <- paste0("Comparison of observed and predicted ",
                         results_ls$study_descs_ls$health_utl_nm_1L_chr,
                         " utility score from longitudinal TTU model using ",
                         results_ls$predr_var_nms_chr %>%
                           paste0(collapse = ", ") %>%
                           stringi::stri_replace_last(fixed = ",",
                                                      " and"),
                         " (A) Density plots of observed and predicted utility scores (",
                         results_ls$ttu_lngl_ls$best_mdls_tb %>%
                           purrr::pmap_chr(~paste0(..1,
                                                   " (",
                                                   ..2,
                                                   ")")) %>%
                           purrr::pluck(1),
                         ") (B) Scatter plots of observed and predicted utility scores by timepoint (",
                         results_ls$ttu_lngl_ls$best_mdls_tb %>%
                           purrr::pmap_chr(~paste0(..1,
                                                   " (",
                                                   ..2,
                                                   ")")) %>%
                           purrr::pluck(1),
                         ") (C) Density plots of observed and predicted utility scores (",
                         ifelse(nrow(results_ls$ttu_lngl_ls$best_mdls_tb)>1,
                                paste0(results_ls$ttu_lngl_ls$best_mdls_tb %>%
                                         purrr::pmap_chr(~paste0(..1,
                                                                 "(",
                                                                 ..2,
                                                                 ")")) %>%
                                         purrr::pluck(2),
                                       ") (D) Scatter plots of observed and predicted utility scores by timepoint (",
                                       results_ls$ttu_lngl_ls$best_mdls_tb %>%
                                         purrr::pmap_chr(~paste0(..1,
                                                                 " (",
                                                                 ..2,
                                                                 ")")) %>%
                                         purrr::pluck(2),")"),
                                ""))
  return(title_1L_chr)
}
make_ds_descvs_ls <- function(candidate_predrs_chr,
                              cohort_descv_var_nms_chr,
                              dictionary_tb,
                              id_var_nm_1L_chr,
                              msrmnt_date_var_nm_1L_chr,
                              round_var_nm_1L_chr,
                              round_vals_chr,
                              maui_item_pfx_1L_chr,
                              utl_wtd_var_nm_1L_chr = "wtd_utl_dbl",
                              utl_unwtd_var_nm_1L_chr = "unwtd_utl_dbl",
                              candidate_covar_nms_chr = NULL,
                              is_fake_1L_lgl = NULL){
  ds_descvs_ls <- list(candidate_covar_nms_chr = candidate_covar_nms_chr,
                       candidate_predrs_chr = candidate_predrs_chr,
                       cohort_descv_var_nms_chr = cohort_descv_var_nms_chr,
                       dictionary_tb = dictionary_tb,
                       id_var_nm_1L_chr = id_var_nm_1L_chr,
                       is_fake_1L_lgl = is_fake_1L_lgl,
                       msrmnt_date_var_nm_1L_chr = msrmnt_date_var_nm_1L_chr,
                       round_var_nm_1L_chr = round_var_nm_1L_chr,
                       round_vals_chr = round_vals_chr,
                       maui_item_pfx_1L_chr = maui_item_pfx_1L_chr,
                       utl_wtd_var_nm_1L_chr = utl_wtd_var_nm_1L_chr,
                       utl_unwtd_var_nm_1L_chr = utl_unwtd_var_nm_1L_chr)
  return(ds_descvs_ls)
}
make_ds_smry_ls <- function(candidate_predrs_chr,
                            candidate_covar_nms_chr,
                            depnt_var_nm_1L_chr,
                            dictionary_tb,
                            id_var_nm_1L_chr,
                            round_var_nm_1L_chr,
                            round_bl_val_1L_chr,
                            predictors_lup){
  ds_smry_ls <- list(candidate_predrs_chr = candidate_predrs_chr,
                     candidate_covar_nms_chr = candidate_covar_nms_chr,
                     depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                     dictionary_tb = dictionary_tb,
                     id_var_nm_1L_chr = id_var_nm_1L_chr,
                     round_var_nm_1L_chr = round_var_nm_1L_chr,
                     round_bl_val_1L_chr = round_bl_val_1L_chr,
                     predictors_lup = predictors_lup)
  return(ds_smry_ls)
}
make_eq5d_ds_dict <- function(data_tb = make_fake_eq5d_ds(),
                              predictors_lup = make_psych_predrs_lup()){
  dictionary_tb <- youthvars::make_tfd_repln_ds_dict_r3() %>%
    dplyr::filter(var_nm_chr %in% names(data_tb)) %>%
    youthvars::make_final_rpln_ds_dict(additions_tb = ready4use::make_pt_ready4_dictionary(var_nm_chr = c("uid",
                                                                                                          "Timepoint", "data_collection_dtm",
                                                                                                          paste0("eq5dq_",c("MO", "SC", "UA", "PD", "AD")),
                                                                                                          "EQ5D_total_dbl",
                                                                                                          "EQ5d_cumulative_dbl",
                                                                                                          predictors_lup$short_name_chr),
                                                                                           var_ctg_chr = c("Identifier",
                                                                                                           rep("Temporal",2),
                                                                                                           rep("Multi-Attribute Utility Instrument Question",5),
                                                                                                           rep("Multi-Attribute Utility Instrument Score",2),
                                                                                                           rep("Clinical",2)),
                                                                                           var_desc_chr = c("Unique identifier",
                                                                                                            "Data collection round",
                                                                                                            "Date of data collection",
                                                                                                            "EQ5D - Mobility domain score",
                                                                                                            "EQ5D - Self-Care domain score",
                                                                                                            "EQ5D - Usual Activities domain score",
                                                                                                            "EQ5D - Pain / Discomfort domain score",
                                                                                                            "EQ5D - Anxiety / Depression domain score",
                                                                                                            "EQ5D - total weighted score",
                                                                                                            "EQ5D - total unweighted score",
                                                                                                            predictors_lup$long_name_chr),
                                                                                           var_type_chr = c("integer",
                                                                                                            "character","date",
                                                                                                            rep("integer",5),
                                                                                                            "double",
                                                                                                            "integer",
                                                                                                            predictors_lup$class_chr
                                                                                           ))) %>%
    dplyr::arrange(var_ctg_chr)
  return(dictionary_tb)
}
make_ethics_text <- function(results_ls){
  text_1L_chr <- ifelse(is.null(results_ls$study_descs_ls$ethics_1L_chr),
                        "",
                        results_ls$study_descs_ls$ethics_1L_chr)
  return(text_1L_chr)
}

make_fake_eq5d_ds <- function(country_1L_chr = "UK",
                              version_1L_chr = "5L",
                              type_1L_chr = "CW",
                              force_attach_1L_lgl = T,
                              prop_with_fup_data_1L_dbl = 0.65,
                              seed_1L_int = 1234,
                              sample_from_1L_int = 10000){
  set.seed(seed_1L_int)
  requireNamespace("eq5d")
  if(force_attach_1L_lgl)
  attachNamespace("eq5d")
  data_tb <- purrr::map(c("MO","SC","UA","PD","AD"),
                        ~list(1:5) %>% stats::setNames(.x)) %>%
    purrr::flatten_dfr() %>%
    tidyr::expand(MO,SC,UA,PD,AD) %>%
    dplyr::mutate(total_eq5d = eq5d::eq5d(., country = country_1L_chr, version = version_1L_chr, type = type_1L_chr))
  k10_lup_tb <- tibble::tibble(k10_int = (9.5 + rexp(sample_from_1L_int,rate=0.18)) %>% purrr::map_int(~as.integer(min(max(.x,10),50))),
                               pred_eq5d_dbl = purrr::map2_dbl(k10_int,rnorm(sample_from_1L_int,0,0.075), ~ predict_utl_from_k10(.x,
                                                                                                                       eq5d_error_1L_dbl = .y)[2]),
                               match_idx_int = purrr::map_int(pred_eq5d_dbl, ~which.min(abs(data_tb$total_eq5d-.x))))
  data_tb <- dplyr::left_join(k10_lup_tb,
                              data_tb %>% dplyr::mutate(match_idx_int = 1:dplyr::n()))
  data_tb <- data_tb %>%
    dplyr::mutate(Timepoint = total_eq5d %>%
                    purrr::map2_chr(runif(sample_from_1L_int),
                                    ~{
                      p_1L_dbl <- ifelse(.x>median(data_tb$total_eq5d),0.55,0.45)
                      ifelse(.y>p_1L_dbl,"FUP","BL")
                    })) %>%
    dplyr::group_by(Timepoint) %>%
    dplyr::mutate(uid = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(uid, Timepoint)
  bl_uids_chr <- data_tb %>% dplyr::filter(Timepoint == "BL") %>% dplyr::pull(uid)
  fup_uids_chr <- sample(bl_uids_chr,round((prop_with_fup_data_1L_dbl * length(bl_uids_chr)),0))
  data_tb <- data_tb %>%
    dplyr::filter(uid %in% bl_uids_chr) %>%
    dplyr::filter(Timepoint == "BL" | uid %in% fup_uids_chr)
  data_tb <- data_tb %>%
    dplyr::mutate(psych_well_int = faux::rnorm_pre(k10_int, mu = 69.9, sd = 9.9, r = -0.56) %>%
                    round(0) %>% purrr::map_int(~min(.x,90) %>% max(18) %>%
                                                  as.integer()))
  data_tb <- data_tb %>% dplyr::select(uid,Timepoint,MO,SC,UA,PD,AD,k10_int,psych_well_int)
  data("replication_popl_tb", package = "youthvars", envir = environment())
  demog_data_tb <- replication_popl_tb %>%
    youthvars::transform_raw_ds_for_analysis() %>%
    dplyr::filter(round == "Baseline") %>%
    dplyr::mutate(uid = 1:dplyr::n()) %>%
    dplyr::select(uid,
                  d_interview_date,
                  d_age,
                  Gender,
                  d_sex_birth_s,
                  d_sexual_ori_s,
                  d_relation_s,
                  d_ATSI,
                  CALD,
                  Region,
                  d_studying_working) %>%
    dplyr::rename(data_collection_dtm =
                    d_interview_date)
  data_tb <- dplyr::left_join(data_tb %>%
                                dplyr::filter(uid %in% demog_data_tb$uid),
                              demog_data_tb) %>%
    dplyr::select(uid,
                  Timepoint,
                  data_collection_dtm,
                  d_age,
                  Gender,
                  d_sex_birth_s,
                  d_sexual_ori_s,
                  d_relation_s,
                  d_ATSI,
                  CALD,
                  Region,
                  d_studying_working,
                  dplyr::everything()) %>%
    dplyr::rename_with(~stringr::str_c("eq5dq_", .), .cols = c("MO", "SC", "UA", "PD", "AD")) %>%
    dplyr::mutate(data_collection_dtm = dplyr::case_when(Timepoint == "FUP" ~ data_collection_dtm + lubridate::days(120),
                                                         T ~ data_collection_dtm))
  data_tb <- data_tb %>%
    dplyr::mutate(d_age = dplyr::case_when(Timepoint == "FUP" ~ d_age %>% purrr::map2_int(runif(nrow(data_tb)), ~ ifelse((.y + 120/365.25)>=1,.x+1L,.x)),
                                           T ~ d_age))
  return(data_tb)
}
make_fake_ts_data <- function (outp_smry_ls,
                               dep_vars_are_NA_1L_lgl = T)
{
    data_tb <- outp_smry_ls$scored_data_tb %>% transform_tb_to_mdl_inp(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
        predr_vars_nms_chr = outp_smry_ls$predr_vars_nms_ls %>% purrr::flatten_chr() %>% unique(),
        id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
        round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr)
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = names(data_tb)[names(data_tb) !=
        outp_smry_ls$id_var_nm_1L_chr], seed = outp_smry_ls$seed_1L_int)
    fk_data_tb <- fk_data_ls$syn %>%
      dplyr::mutate(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr) := as.character(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr))) %>%
      dplyr::group_by(!!rlang::sym(outp_smry_ls$id_var_nm_1L_chr)) %>%
      dplyr::mutate(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr) := !!rlang::sym(outp_smry_ls$round_var_nm_1L_chr) %>%
                      transform_timepoint_vals(timepoint_levels_chr = outp_smry_ls$scored_data_tb %>%
                                                 dplyr::pull(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr)) %>%
                                                 unique(),
                                               bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr)
                      )%>%
      dplyr::ungroup()
    if(dep_vars_are_NA_1L_lgl){
      dep_vars_chr <- names(fk_data_tb)[names(fk_data_tb) %>%
                                              purrr::map_lgl(~startsWith(.x, outp_smry_ls$depnt_var_nm_1L_chr))]
      fk_data_tb <- fk_data_tb %>% dplyr::mutate(dplyr::across(dplyr::all_of(dep_vars_chr),
                                                                   ~NA_real_))
    }
    return(fk_data_tb)
}
make_folds_ls <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", folds_1L_int = 10L)
{
    folds_ls <- caret::createFolds(data_tb %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
                                   k = folds_1L_int, list = TRUE, returnTrain = FALSE)
    return(folds_ls)
}
make_funding_text <- function(results_ls){
  text_1L_chr <- ifelse(is.null(results_ls$study_descs_ls$funding_1L_chr),
                        "",
                        results_ls$study_descs_ls$funding_1L_chr)
  return(text_1L_chr)
}
make_header_yaml_args_ls <- function(authors_tb,
                                     institutes_tb,
                                     title_1L_chr,
                                     keywords_chr,
                                     fl_nm_1L_chr = "header_common.yaml",
                                     use_fake_data_1L_lgl = F){
  if(!use_fake_data_1L_lgl){
    header_yaml_args_ls <- list(authors_tb = authors_tb,
                                institutes_tb = institutes_tb,
                                fl_nm_1L_chr = "header_common.yaml",
                                title_1L_chr = title_1L_chr,
                                keywords_chr = keywords_chr)
  }else{
    data("authors_tb", package = "ready4show", envir = environment())
    data("institutes_tb", package = "ready4show", envir = environment())
    header_yaml_args_ls <- make_header_yaml_args_ls(authors_tb = authors_tb,
                                                    institutes_tb = institutes_tb,
                                                    title_1L_chr = "A hypothetical study using fake data for instructional purposes only",
                                                    keywords_chr = c("this","is","a","replication","using","fake","data","do", "not","cite"))
  }
  return(header_yaml_args_ls)
}
make_hlth_utl_and_predrs_ls <- function(outp_smry_ls,
                                        descv_tbls_ls,
                                        nbr_of_digits_1L_int = 2L,
                                        old_nms_chr = NULL,
                                        new_nms_chr = NULL){
  ranked_predrs_ls <- make_ranked_predrs_ls(descv_tbls_ls,
                                            old_nms_chr = old_nms_chr,
                                            new_nms_chr = new_nms_chr)
  var_nm_1L_chr <- descv_tbls_ls$ds_descvs_ls$dictionary_tb %>%
    ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "var_nm_chr",
                                match_value_xx = outp_smry_ls$depnt_var_nm_1L_chr,
                                target_var_nm_1L_chr = "var_desc_chr",
                                evaluate_lgl = F) %>% as.vector()


  hlth_utl_and_predrs_ls = list(bl_hu_mean_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                              match_value_xx = var_nm_1L_chr,
                                                              target_var_nm_1L_chr = paste0(ds_descvs_ls$round_vals_chr[1],
                                                                                            "_val_1_dbl"),
                                                              evaluate_lgl = F) %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                bl_hu_sd_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                              match_value_xx = var_nm_1L_chr,
                                                              target_var_nm_1L_chr = paste0(ds_descvs_ls$round_vals_chr[1],
                                                                                            "_val_2_ls"),
                                                              evaluate_lgl = F) %>%
                                  stringr::str_remove("\\(") %>%
                                  stringr::str_remove("\\)") %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                fup_hu_mean_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                              match_value_xx = var_nm_1L_chr,
                                                              target_var_nm_1L_chr = paste0(ds_descvs_ls$round_vals_chr[2],
                                                                                            "_val_1_dbl"),
                                                              evaluate_lgl = F) %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                fup_hu_sd_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                              match_value_xx = var_nm_1L_chr,
                                                              target_var_nm_1L_chr = paste0(ds_descvs_ls$round_vals_chr[2],
                                                                                            "_val_2_ls"),
                                                              evaluate_lgl = F) %>%
                                  stringr::str_remove("\\(") %>%
                                  stringr::str_remove("\\)") %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                predrs_nartv_seq_chr = ranked_predrs_ls$unranked_predrs_chr,
                                cor_seq_dscdng_chr =  ranked_predrs_ls$ranked_predrs_chr)
  return(hlth_utl_and_predrs_ls)

}
make_indpnt_predrs_lngl_tbls_ref <- function(params_ls){
  results_ls <- params_ls$results_ls
  n_mdls_1L_int <- length(results_ls$ttu_lngl_ls$best_mdls_tb$model_type)
  text_1L_chr <- paste0(ifelse(params_ls$output_type_1L_chr=="Word","","Table"),
                        ifelse(n_mdls_1L_int < 3,
                               " \\@ref(tab:cfscl)",
                               paste0("s ",
                                      paste0("\\@ref(tab:cfscl",1:n_mdls_1L_int,")", collapse = ", ") %>%
                                        stringi::stri_replace_last(fixed = ",", " and")
                               ))
  )
  return(text_1L_chr)
}
make_indpnt_predrs_lngl_tbl_title <- function (results_ls, ref_1L_int = 1)
{
  title_1L_chr <- paste0("Estimated coefficients for single predictor longitudinal TTU models using ",
                         results_ls$ttu_lngl_ls$best_mdls_tb[[ref_1L_int,
                                                              "model_type"]], " (", results_ls$ttu_lngl_ls$best_mdls_tb[[ref_1L_int,
                                                                                                                         "link_and_tfmn_chr"]], ")")
  return(title_1L_chr)
}
make_input_params <- function(ds_tb,
                              ds_descvs_ls,
                              header_yaml_args_ls,
                              maui_params_ls,
                              predictors_lup,
                              control_ls = NULL,
                              dv_ds_nm_and_url_chr = NULL,
                              iters_1L_int = 4000L,
                              mdl_smry_ls = make_mdl_smry_ls(),
                              output_format_ls = make_output_format_ls(),
                              path_params_ls = NULL,
                              prefd_covars_chr = NULL,
                              prefd_mdl_types_chr = NULL,
                              prior_ls = NULL,
                              seed_1L_int = 12345,
                              scndry_anlys_params_ls = NULL,
                              write_new_dir_1L_lgl = T){
  path_params_ls <- make_path_params_ls(use_fake_data_1L_lgl = ds_descvs_ls$is_fake_1L_lgl,
                                        dv_ds_nm_and_url_chr = dv_ds_nm_and_url_chr,
                                        write_new_dir_1L_lgl = write_new_dir_1L_lgl)
  params_ls_ls <- make_analysis_core_params_ls(ds_descvs_ls = ds_descvs_ls,
                                               output_format_ls = output_format_ls,
                                               predictors_lup = predictors_lup,
                                               prefd_covars_chr = prefd_covars_chr,
                                               prefd_mdl_types_chr = prefd_mdl_types_chr,
                                               mdl_smry_ls = mdl_smry_ls,
                                               control_ls = control_ls,
                                               iters_1L_int = iters_1L_int,
                                               prior_ls = prior_ls,
                                               seed_1L_int = seed_1L_int) %>%
    make_valid_params_ls_ls(ds_tb = ds_tb,
                            maui_params_ls = maui_params_ls,
                            path_params_ls = path_params_ls)
  params_ls_ls$header_yaml_args_ls <- header_yaml_args_ls
  params_ls_ls$output_format_ls <- output_format_ls
  params_ls_ls$scndry_anlys_params_ls <- scndry_anlys_params_ls
  return(params_ls_ls)
}
make_knit_pars_ls <- function (rltv_path_to_data_dir_1L_chr, mdl_types_chr, predr_vars_nms_ls,
                               output_type_1L_chr = "HTML", mdl_types_lup = NULL, plt_types_lup = NULL,
                               plt_types_chr = NA_character_, section_type_1L_chr = "#")
{
    if (is.null(mdl_types_lup))
        utils::data(mdl_types_lup, envir = environment())
    if (is.null(plt_types_lup))
        utils::data(plt_types_lup, envir = environment())
  if(is.na(plt_types_chr)){
    plt_types_chr <- plt_types_lup$short_name_chr
  }
   paths_to_all_data_fls_chr <- list.files(rltv_path_to_data_dir_1L_chr,full.names = T)
    lab_idx_dbl <- 1:(length(mdl_types_chr) * length(predr_vars_nms_ls))
    knit_pars_ls <- purrr::pmap(list(predr_vars_nms_ls,
                                     split(lab_idx_dbl,
                                           ceiling(seq_along(lab_idx_dbl)/length(mdl_types_chr))),
                                     make_mdl_nms_ls(predr_vars_nms_ls,
                                                     mdl_types_chr = mdl_types_chr)),
                                ~{
                                  mdl_nms_chr <- ..3
                                  mdl_data_paths_ls <- mdl_nms_chr %>%
                                    purrr::map(~paths_to_all_data_fls_chr[stringr::str_detect(paths_to_all_data_fls_chr,.x)]) %>%
                                    stats::setNames(mdl_nms_chr)
                                  paths_to_mdls_chr <- mdl_data_paths_ls %>%
                                    purrr::map_chr(~ifelse(identical(.x[endsWith(.x,".RDS")],character(0)),
                                                           NA_character_,
                                                           .x[endsWith(.x,".RDS")]
                                                           )) %>% unname()
                                  paths_to_mdl_plts_ls <- mdl_data_paths_ls %>%
                                    purrr::map(~{
                                      paths_to_all_plots_chr <- .x[endsWith(.x,".png")]
                                      plt_types_chr %>%
                                        purrr::map(~ {
                                          sfx_1L_chr <- paste0(.x,".png")
                                                     paths_to_all_plots_chr[paths_to_all_plots_chr %>%
                                                                                 purrr::map_lgl(~endsWith(.x,sfx_1L_chr))][paths_to_all_plots_chr[paths_to_all_plots_chr %>%
                                                                                                                                                        purrr::map_lgl(~endsWith(.x,sfx_1L_chr))] %>%
                                                                                                                              nchar() %>%
                                                                                                                              which.min()]
                                                   }) %>%
                                        purrr::flatten_chr() %>% unique()
                                      })
                                  mdl_ttls_chr <- paste0(..1[1], ifelse(is.na(..1[2]),
                                                                        "",
                                                                        paste(" with ", ..1[2])),
                                                         " ",
                                                         mdl_types_chr %>% purrr::map_chr(~ready4fun::get_from_lup_obj(mdl_types_lup,
                                                                                     match_var_nm_1L_chr = "short_name_chr",
                                                                                     match_value_xx = .x,
                                                                                     target_var_nm_1L_chr = "long_name_chr",
                                                                                     evaluate_lgl = F)))
                                  section_ttls_chr <- paste0(section_type_1L_chr, " ", mdl_ttls_chr)
                                  plt_nms_ls <- paths_to_mdl_plts_ls %>% purrr::map2(mdl_ttls_chr,~{
                                    paths_to_mdl_plts_chr <- .x
                                    mdl_ttl_1L_chr <- .y
                                    transform_1L_lgl <- paths_to_mdl_plts_chr %>% endsWith("_coefs.png") %>% any()
                                    if(paths_to_mdl_plts_chr %>% endsWith("_hetg.png") %>% any())
                                      transform_1L_lgl <- F # ALL OF THIS TFMN LOGIC CAN BE BINNED WHEN NEW CATEGORY OF PLOT TYPE IS ADDED.
                                    plt_types_chr %>% purrr::map(~{
                                      if(endsWith(paths_to_mdl_plts_chr,paste0("_",.x,".png")) %>% any()){
                                        paste0(mdl_ttl_1L_chr,
                                               " ",
                                               ifelse(transform_1L_lgl & .x == "coefs",
                                                      "population and group level effects",
                                                      ready4fun::get_from_lup_obj(plt_types_lup,
                                                                           match_var_nm_1L_chr = "short_name_chr",
                                                                           match_value_xx = .x,
                                                                           target_var_nm_1L_chr = "long_name_chr",
                                                                           evaluate_lgl = F)))

                                      }else{
                                        character(0)
                                      }}) %>% purrr::flatten_chr()
                                  })
                                  list(plt_nms_ls = plt_nms_ls,
                                       paths_to_mdls_chr = paths_to_mdls_chr,
                                       tbl_captions_chr = mdl_ttls_chr,
                                       label_stubs_chr = paste0("lab", ..2),
                                       output_type_1L_chr = output_type_1L_chr,
                                       section_ttls_chr = section_ttls_chr,
                                       paths_to_mdl_plts_ls = paths_to_mdl_plts_ls)
                                  }) %>%
      stats::setNames(predr_vars_nms_ls %>% purrr::map_chr(~paste(.x,collapse="_")))
    return(knit_pars_ls)
}
make_lngl_ttu_with_covars_text <- function(results_ls){
  text_1L_chr <- ifelse((is.na(results_ls$ttu_lngl_ls$incld_covars_chr) | length(results_ls$ttu_lngl_ls$incld_covars_chr) == 0),
                        "",
                        paste0("We also evaluated models with ",
                               results_ls$ttu_lngl_ls$incld_covars_chr %>%
                                 paste0(collapse = ", ") %>%
                                 stringi::stri_replace_last(fixed = ",", " and") %>%
                                 paste0(" at baseline"),
                               " and ",
                               results_ls$ttu_lngl_ls$incld_covars_chr %>%
                                 paste0(collapse = ", ") %>%
                                 stringi::stri_replace_last(fixed = ",", " and") %>%
                                 paste0(" change from baseline"),
                               " added to ",
                               names(results_ls$study_descs_ls$predr_ctgs_ls) %>%
                                 tolower() %>%
                                 paste0(collapse = ", ") %>%
                                 stringi::stri_replace_last(fixed = ",", " and"),
                               " predictors"))
  return(text_1L_chr )
}
make_lngl_ttu_r2_text <- function(results_ls,
                                  part_int = 1){
  if(1 %in% part_int){
    text_1L_chr <- ifelse(length(results_ls$candidate_predrs_chr) == 1,
                          "",
                          paste0("In ",
                                 ifelse(length(results_ls$ttu_lngl_ls$best_mdls_tb$name_chr %>% unique())>1,
                                        "",
                                        ifelse(nrow(results_ls$ttu_lngl_ls$best_mdls_tb)==2,
                                               "both ",
                                               "all ")),
                                 results_ls$ttu_lngl_ls$best_mdls_tb$model_type %>%
                                   paste0(collapse = ", ") %>%
                                   stringi::stri_replace_last(fixed = ",", " and"),
                                 " models, the prediction models using ",
                                 ifelse(length(results_ls$ttu_lngl_ls$best_mdls_tb$name_chr %>% unique())==1,
                                        results_ls$ttu_lngl_ls$best_mdls_tb$name_chr %>% unique(),
                                        results_ls$ttu_lngl_ls$best_mdls_tb$name_chr %>%
                                          paste0(collapse = ", ") %>%
                                          stringi::stri_replace_last(fixed = ",", " and") %>%
                                          paste0(" respectively")),
                                 " had the highest R^2^ (",
                                 results_ls$ttu_lngl_ls$best_mdls_tb$r2_dbl %>%
                                   stringr::str_squish() %>%
                                   paste0(collapse = ", ") %>%
                                   stringi::stri_replace_last(fixed = ",", " and"),
                                 ")"))
  }else{
    text_1L_chr <- ""
  }
  if(2 %in% part_int){
    mdl_types_chr <- results_ls$ttu_lngl_ls$best_mdls_tb$model_type
    mdl_r2_var_nms_chr <- intersect(results_ls$tables_ls$ind_preds_coefs_tbl %>% names(),
                                    c("R2_OLS_CLL","R2_GLM_GSN_LOG"))
    mdl_r2_mins_dbl <- mdl_r2_var_nms_chr %>%
      purrr::map_dbl(~results_ls$tables_ls$ind_preds_coefs_tbl %>%
                       dplyr::filter(!is.na(!!rlang::sym(.x))) %>%
                       dplyr::pull(!!rlang::sym(.x)) %>%
                       as.numeric() %>%
                       min())
    mdl_r2_maxs_dbl <- mdl_r2_var_nms_chr %>%
      purrr::map_dbl(~results_ls$tables_ls$ind_preds_coefs_tbl %>%
                       dplyr::filter(!is.na(!!rlang::sym(.x))) %>%
                       dplyr::pull(!!rlang::sym(.x)) %>%
                       as.numeric() %>%
                       max())
    text_1L_chr <- paste0(text_1L_chr,
                          ifelse(1 %in% part_int,". ",""),
                          list(mdl_types_chr,
                               mdl_r2_mins_dbl,
                               mdl_r2_maxs_dbl,
                               1:length(mdl_types_chr)) %>%
                            purrr::pmap_chr(~ paste0(ifelse(..4 == 1,"R^2^ was ",""),
                                                     ifelse(length(results_ls$candidate_predrs_chr) == 1,
                                                            paste0(,..2, " for the "),
                                                            paste0(ifelse(..2==..3,
                                                                          paste0("",..2),
                                                                          paste0("between ",..2," and ",..3," for all ")))),
                                                     ..1,
                                                     ifelse(length(results_ls$candidate_predrs_chr) == 1,
                                                            paste0(""),
                                                            paste0("s"))
                            )) %>%
                            paste0(collapse = ", ") %>%
                            stringi::stri_replace_last(fixed = ",", " and"),
                          ".")
  }
  return(text_1L_chr)
}

make_maui_params_ls <- function(maui_itm_short_nms_chr,
                                maui_domains_pfcs_1L_chr = NULL,
                                maui_scoring_fn = NULL,
                                short_and_long_nm = NULL,
                                utl_min_val_1L_dbl = -1){
  if(is.null(maui_scoring_fn)){
    maui_scoring_fn <- function(data_tb,
                                maui_item_pfx_1L_chr,
                                id_var_nm_1L_chr,
                                utl_wtd_var_nm_1L_chr,
                                utl_unwtd_var_nm_1L_chr){
      data_tb %>%
        dplyr::mutate(`:=`(!!rlang::sym(utl_unwtd_var_nm_1L_chr),                                          rowSums(dplyr::across(dplyr::starts_with(maui_item_pfx_1L_chr))))) %>%
        dplyr::filter(!is.na(!!rlang::sym(utl_unwtd_var_nm_1L_chr)))
    }
  }
  maui_params_ls <- list(maui_domains_pfcs_1L_chr = maui_domains_pfcs_1L_chr,
                         maui_itm_short_nms_chr = maui_itm_short_nms_chr,
                         maui_scoring_fn = maui_scoring_fn,
                         short_and_long_nm = short_and_long_nm,
                         utl_min_val_1L_dbl = utl_min_val_1L_dbl)
  return(maui_params_ls)
}
make_mdl <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF",
    predr_var_nm_1L_chr, covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF",
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    if (is.null(start_1L_chr)) {
        start_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
            target_var_nm_1L_chr = "start_chr", evaluate_lgl = F)
    }
    if (!is.na(control_1L_chr)) {
        idx_1L_int <- 1 + stringi::stri_locate_last_fixed(mdl_type_1L_chr,
            "_")[1, 1] %>% as.vector()
        link_1L_chr <- get_link_from_tfmn(stringr::str_sub(mdl_type_1L_chr, start = idx_1L_int))
    }
    mdl_1L_chr <- paste0(ready4fun::get_from_lup_obj(mdl_types_lup,
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
        target_var_nm_1L_chr = "fn_chr", evaluate_lgl = F), "(",
        transform_depnt_var_nm(depnt_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr),
        " ~ ", predr_var_nm_1L_chr, ifelse(is.na(covar_var_nms_chr[1]),
            "", paste0(" + ", paste0(covar_var_nms_chr, collapse = " + "))),
        ", data = data_tb", ifelse(!is.na(ready4fun::get_from_lup_obj(mdl_types_lup,
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
            target_var_nm_1L_chr = "family_chr", evaluate_lgl = F)),
            paste0(", family = ", ready4fun::get_from_lup_obj(mdl_types_lup,
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                target_var_nm_1L_chr = "family_chr", evaluate_lgl = F)),
            ""), ifelse((!is.na(start_1L_chr)|(is.na(start_1L_chr) & !is.na(control_1L_chr))), ", ", ""), ifelse(!is.na(control_1L_chr),
            paste0("link=\"", link_1L_chr, "\",control=", control_1L_chr,
                "("), ""), ifelse(!is.na(start_1L_chr), paste0("start=c(",
            start_1L_chr, ")"), ""), ifelse(!is.na(control_1L_chr),
            ")", ""), ")")
    model_mdl <- eval(parse(text = mdl_1L_chr))
    return(model_mdl)
}
make_mdl_coef_ratio_ls <- function(mdl_ingredients_ls,
                                   predr_ctgs_ls = NULL){
  predrs_chr <- mdl_ingredients_ls$predictors_lup$short_name_chr
  mdl_type_chr <- mdl_ingredients_ls$mdls_smry_tb$Model %>%
    unique() %>%
    purrr::map_chr(~get_mdl_type_from_nm(.x,
                                         mdl_types_lup = mdl_ingredients_ls$mdl_types_lup)) %>%
    unique()
  main_mdls_ls <- predrs_chr %>% purrr::map(~paste0(paste0(.x,"_1_"),
                                                    mdl_type_chr))
  ratios_ls <- main_mdls_ls %>%
    purrr::map2(predrs_chr, ~{
      mdls_chr <- .x
      predr_1L_chr <- .y
      mdls_chr %>% purrr::map_dbl(~{
        coefs_dbl <- mdl_ingredients_ls$mdls_smry_tb %>% dplyr::filter(Model %in% .x) %>%
          dplyr::filter(Parameter %in% paste0(predr_1L_chr,(c(" baseline"," change")))) %>%
          dplyr::pull(Estimate)
        coefs_dbl[2]/coefs_dbl[1]
      })
    }
    ) %>% stats::setNames(predrs_chr)
  mdl_coef_ratios_ls = list(mean_ratios_dbl = ratios_ls %>% purrr::map_dbl(~mean(.x)))
  if(!is.null(predr_ctgs_ls)){
    append_ls <- purrr::map(predr_ctgs_ls,
                            ~ mdl_coef_ratios_ls$mean_ratios_dbl[predrs_chr %in% .x]) %>%
      stats::setNames(names(predr_ctgs_ls))
    mdl_coef_ratios_ls <- append(mdl_coef_ratios_ls, append_ls)
  }
  return(mdl_coef_ratios_ls)
}
make_mdl_coef_range_text <- function(coef_ratios_dbl,
                                     nbr_of_digits_1L_int = 2L){
  if(length(coef_ratios_dbl) == 1){
    coef_range_text_chr <- as.character(round(coef_ratios_dbl, nbr_of_digits_1L_int))
  }else{
    min_1L_dbl <- round(min(coef_ratios_dbl), nbr_of_digits_1L_int)
    max_1L_dbl <- round(max(coef_ratios_dbl), nbr_of_digits_1L_int)
    coef_range_text_chr <- paste0("between ", min_1L_dbl," and ", max_1L_dbl)
  }
  return(coef_range_text_chr)
}
make_mdl_nms_ls <- function (predr_vars_nms_ls, mdl_types_chr)
{
    mdl_nms_ls <- purrr::map2(predr_vars_nms_ls, make_unique_ls_elmt_idx_int(predr_vars_nms_ls),
        ~paste0(.x[1], "_", ifelse(is.na(.x[2]), "", paste0(.x[2],
            "_")), .y, "_", mdl_types_chr))
    return(mdl_nms_ls)
}
make_mdl_smry_elmt_tbl <- function (mat, ctg_chr)
{
    tb <- mat %>% tibble::as_tibble() %>% dplyr::mutate(Parameter = rownames(mat)) %>%
        dplyr::select(Parameter, dplyr::everything())
    mdl_elmt_sum_tb <- tb %>% dplyr::filter(F) %>% tibble::add_case(Parameter = ctg_chr) %>%
        dplyr::bind_rows(tb)
    return(mdl_elmt_sum_tb)
}
make_mdl_smry_ls <- function(mdl_types_lup = get_cndts_for_mxd_mdls(),
                             mdl_types_chr = NULL,
                             choose_from_pfx_chr = NULL,
                             folds_1L_int = 10L,
                             max_nbr_of_boruta_mdl_runs_int = 300L){
  if(is.null(mdl_types_lup))
    data("mdl_types_lup",package = "TTU", envir = environment())
  if(is.null(mdl_types_chr))
    mdl_types_chr <- mdl_types_lup$short_name_chr
  if(is.null(choose_from_pfx_chr))
    choose_from_pfx_chr <- stringr::word(mdl_types_lup$short_name_chr,1,sep="\\_") %>% unique()
  mdl_smry_ls <- list(mdl_types_lup = mdl_types_lup,
                      mdl_types_chr = mdl_types_chr,
                      choose_from_pfx_chr = choose_from_pfx_chr,
                      folds_1L_int = folds_1L_int,
                      max_nbr_of_boruta_mdl_runs_int = max_nbr_of_boruta_mdl_runs_int)
  return(mdl_smry_ls)
}
make_mdls_ls <- function(outp_smry_ls,
                         mdls_tb){
  mdls_chr <- mdls_tb$Model %>% unique()
  mdls_ls <- outp_smry_ls$prefd_mdl_types_chr %>%
    purrr::map(~mdls_chr[mdls_chr %>% endsWith(.x)])
  return(mdls_ls)
}
make_mdls_smry_tbls_ls <- function(outp_smry_ls,
                                   nbr_of_digits_1L_int = 2L){
  mdls_smry_tb <- outp_smry_ls$mdls_smry_tb  %>%
    dplyr::mutate(dplyr::across(c("Estimate","SE"),
                                ~ round(.x,nbr_of_digits_1L_int) %>%
                                  format(nsmall = nbr_of_digits_1L_int))) %>%
    dplyr::mutate(`95% CI` = `95% CI` %>% transform_chr_digit_pairs())
  rownames(mdls_smry_tb) <- NULL
  indpt_predrs_mdls_tb <- mdls_smry_tb %>% dplyr::filter(Model %in% (paste0(outp_smry_ls$predr_cmprsn_tb$predr_chr,"_1_") %>% purrr::map(~paste0(.x,
                                                                                                                                                 outp_smry_ls$prefd_mdl_types_chr)) %>% purrr::flatten_chr()))
  covar_mdls_tb <- mdls_smry_tb %>% dplyr::filter(!Model %in% indpt_predrs_mdls_tb$Model)
  mdl_types_chr <- indpt_predrs_mdls_tb$Model %>% purrr::map_chr(~get_mdl_type_from_nm(.x,
                                                                      mdl_types_lup = outp_smry_ls$mdl_types_lup)) %>%
    unique()
  prefd_predr_mdl_smry_tb <- mdl_types_chr %>%
    purrr::map_dfr(~{
      mdl_type_1L_chr <- .x
     mdl_type_smry_tb <- indpt_predrs_mdls_tb %>%
        dplyr::filter(Model %>% purrr::map_lgl(~endsWith(.x, mdl_type_1L_chr)))
     max_r2_dbl <- mdl_type_smry_tb %>%
       dplyr::filter(Parameter == "R2") %>%
       dplyr::pull(Estimate) %>%
       as.numeric() %>%
       abs() %>%
       max()
     prefd_mdl_1L_chr <- mdl_type_smry_tb %>%
       dplyr::filter(Parameter == "R2") %>%
       dplyr::filter(as.numeric(Estimate) == max_r2_dbl) %>%
       dplyr::pull(Model) %>%
       purrr::pluck(1)
     mdl_type_smry_tb %>%
       dplyr::filter(Model == prefd_mdl_1L_chr)
    })
  mdls_smry_tbls_ls <- list(indpt_predrs_mdls_tb = indpt_predrs_mdls_tb,
                            covar_mdls_tb = covar_mdls_tb,
                            prefd_predr_mdl_smry_tb = prefd_predr_mdl_smry_tb)
  return(mdls_smry_tbls_ls)
}
make_mdl_type_smry_tbl <- function(mdls_tb,
                                   mdl_nms_chr,
                                   mdl_type_1L_chr,
                                   add_mdl_nm_sfx_1L_lgl = T){
  mdl_type_smry_tbl_tb <- mdl_nms_chr %>%
    purrr::map_dfr(~make_sngl_mdl_smry_tb(mdls_tb,
                                          mdl_nm_1L_chr = .x,
                                          mdl_type_1L_chr = mdl_type_1L_chr,
                                          add_mdl_nm_sfx_1L_lgl = add_mdl_nm_sfx_1L_lgl))
  return(mdl_type_smry_tbl_tb)
}
make_mdl_desc_lines <- function(outp_smry_ls,
                                mdl_nm_1L_chr,
                                output_type_1L_chr = "PDF"){
  mdl_smry_tb <- outp_smry_ls$mdls_smry_tb %>%
    dplyr::filter(Model == mdl_nm_1L_chr)
  predictors_chr <- mdl_smry_tb$Parameter[!mdl_smry_tb$Parameter %in% c("SD (Intercept)","Intercept","R2","RMSE","Sigma")] %>%
    purrr::map_chr(~stringr::str_remove(.x," baseline") %>% stringr::str_remove(" change")) %>% unique()
  predictors_desc_chr <- predictors_chr %>%
    purrr::map_chr(~{
      scaling_1L_dbl <- ready4fun::get_from_lup_obj(outp_smry_ls$predictors_lup,
                                                    match_value_xx = .x,
                                                    match_var_nm_1L_chr = "short_name_chr",
                                                    target_var_nm_1L_chr = "mdl_scaling_dbl",
                                                    evaluate_lgl = F)
      paste0(.x,
             " (",
             ready4fun::get_from_lup_obj(outp_smry_ls$dictionary_tb,
                                         match_value_xx = .x,
                                         match_var_nm_1L_chr = "var_nm_chr",
                                         target_var_nm_1L_chr = "var_desc_chr",
                                         evaluate_lgl = F),
             ifelse(scaling_1L_dbl == 1,
                    "",
                    paste0(" (multiplied by ", scaling_1L_dbl,")")),
             ")")
    })
  if(length(predictors_desc_chr) > 1)
    predictors_desc_chr <- paste0(c(paste0("\n - ",
                                           predictors_desc_chr[-length(predictors_desc_chr)],
                                           collapse = ";"),
                                    paste0("\n - ",predictors_desc_chr[length(predictors_desc_chr)])),
                                  collapse = "; and")


  mdl_desc_lines_chr <- paste0(paste0("This model predicts values at two timepoints for ",
                                      ready4fun::get_from_lup_obj(outp_smry_ls$dictionary_tb,
                                                                  match_value_xx = outp_smry_ls$depnt_var_nm_1L_chr,
                                                                  match_var_nm_1L_chr = "var_nm_chr",
                                                                  target_var_nm_1L_chr = "var_desc_chr",
                                                                  evaluate_lgl = F),
                                      ". The predictor variables are ",
                                      "baseline values and subsequent changes in ",
                                      collapse = ""), predictors_desc_chr,". ",
                               "The catalogue reference for this model is ",
                               ifelse(output_type_1L_chr == "PDF",
                                      paste0("\\texttt{\\detokenize{",mdl_nm_1L_chr,"}}"),
                                      mdl_nm_1L_chr),
                               ".")
  return(mdl_desc_lines_chr)
}
make_nbr_at_fup_text <- function(results_ls){
  nbr_at_fup_1L_chr <- paste0("There were ",
                              results_ls$cohort_ls$n_fup_1L_dbl,
                              " participants (",
                              (results_ls$cohort_ls$n_fup_1L_dbl / results_ls$cohort_ls$n_inc_1L_dbl * 100) %>% round(1),
                              "%) who completed ",
                              results_ls$study_descs_ls$health_utl_nm_1L_chr,
                              " questions at the follow-up survey ",
                              results_ls$study_descs_ls$time_btwn_bl_and_fup_1L_chr,
                              " after baseline assessment."
  )
  return(nbr_at_fup_1L_chr)
}
make_nbr_included_text <- function(results_ls){
  paste0(ifelse(results_ls$cohort_ls$n_inc_1L_dbl == results_ls$cohort_ls$n_all_1l_dbl,
                "all ",
                paste0(results_ls$cohort_ls$n_inc_1L_dbl,
                       " out of the ")),
         results_ls$cohort_ls$n_all_1l_dbl,
         " participants with complete ",
         results_ls$study_descs_ls$health_utl_nm_1L_chr,
         " data")
}
make_output_format_ls <- function(manuscript_outp_1L_chr = "Word",
                                  manuscript_digits_1L_int = 2L,
                                  supplementary_outp_1L_chr = "PDF",
                                  supplementary_digits_1L_int = 2L){
  output_format_ls <- list(manuscript_outp_1L_chr = manuscript_outp_1L_chr,
                           manuscript_digits_1L_int = manuscript_digits_1L_int,
                           supplementary_outp_1L_chr = supplementary_outp_1L_chr,
                           supplementary_digits_1L_int = supplementary_digits_1L_int)
  return(output_format_ls)
}
make_path_params_ls <- function(path_to_data_from_top_level_chr = NULL,
                                path_from_top_level_1L_chr = NULL,
                                path_to_current_1L_chr = NULL,
                                dv_ds_nm_and_url_chr = NULL,
                                write_new_dir_1L_lgl = F,
                                use_fake_data_1L_lgl = F,
                                R_fl_nm_1L_chr = 'aaaaaaaaaa.txt'){
  if(is.null(path_to_data_from_top_level_chr))
    path_to_data_from_top_level_chr <- ifelse(use_fake_data_1L_lgl,
                                              "fake_data.rds",
                                              "data.rds")
  if(is.null(path_from_top_level_1L_chr)){
    path_from_top_level_1L_chr <- normalizePath("../") %>% strsplit("\\\\") %>% purrr::pluck(1) %>% tail(1)
  }
  if(is.null(path_to_current_1L_chr)){
    path_to_current_1L_chr <- normalizePath(".") %>% strsplit("\\\\") %>% purrr::pluck(1) %>% tail(1)
  }
  path_params_ls <- list(path_from_top_level_1L_chr = path_from_top_level_1L_chr,
                         path_to_data_from_top_level_chr = path_to_data_from_top_level_chr,
                         path_to_current_1L_chr = path_to_current_1L_chr,
                         dv_ds_nm_and_url_chr = dv_ds_nm_and_url_chr)
  if(write_new_dir_1L_lgl)
  path_params_ls$paths_ls <- write_main_oupt_dir(path_params_ls,
                                                 use_fake_data_1L_lgl = use_fake_data_1L_lgl,
                                                 R_fl_nm_1L_chr = R_fl_nm_1L_chr)

  return(path_params_ls)
}
make_paths_to_ss_plts_ls <- function(output_data_dir_1L_chr,
                                     outp_smry_ls,
                                     additional_paths_chr = "/dens_and_sctr.png"){
  paths_to_ss_plts_ls = list(combined_utl = paste0(output_data_dir_1L_chr,"/_Descriptives/combined_utl.png"),
                             composite = paste0(output_data_dir_1L_chr,additional_paths_chr[1]),
                             items = paste0(output_data_dir_1L_chr,"/_Descriptives/qstn_rspns.png"),
                             density = paste0(output_data_dir_1L_chr,"/",outp_smry_ls$file_paths_chr[outp_smry_ls$file_paths_chr %>% purrr::map_lgl(~stringr::str_detect(.x,"A_TFMN_CMPRSN_DNSTY"))]),#paste0(output_data_dir_1L_chr,"/A_TFMN_CMPRSN_DNSTY.png"), #Update
                             importance = paste0(output_data_dir_1L_chr,"/",outp_smry_ls$file_paths_chr[outp_smry_ls$file_paths_chr %>% purrr::map_lgl(~stringr::str_detect(.x,"B_PRED_CMPRSN_BORUTA_VAR_IMP"))]))
  return(paths_to_ss_plts_ls)
}
make_predr_ctgs_ls <- function(outp_smry_ls,
                               include_idx_int = NULL){
  predictors_chr <- outp_smry_ls$predr_vars_nms_ls %>% purrr::flatten_chr() %>% unique()
  categories_chr <- predictors_chr %>%
    purrr::map_chr(~outp_smry_ls$dictionary_tb %>%
                     ready4fun::get_from_lup_obj(match_value_xx = .x,
                                                 match_var_nm_1L_chr = "var_nm_chr",
                                                 target_var_nm_1L_chr = "var_ctg_chr",
                                                 evaluate_lgl = F)) %>% unique()
  predr_ctgs_ls <- categories_chr %>% purrr::map(~outp_smry_ls$dictionary_tb %>%
                                                   ready4use::remove_labels_from_ds() %>%
                                                   dplyr::filter(var_ctg_chr == .x) %>%
                                                   dplyr::pull(var_nm_chr)) %>%
    stats::setNames(categories_chr)
  if(!is.null(include_idx_int)){
    predr_ctgs_ls <- include_idx_int %>%
      purrr::map(~predr_ctgs_ls %>% purrr::pluck(.x)) %>%
      stats::setNames(categories_chr[include_idx_int])
  }
  return(predr_ctgs_ls)
}
make_predn_ds_with_one_predr <- function(model_mdl, depnt_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF",
    predr_var_nm_1L_chr, predr_vals_dbl, predn_type_1L_chr = NULL)
{
    predn_ds_tb <- tibble::tibble(`:=`(!!rlang::sym(predr_var_nm_1L_chr),
        predr_vals_dbl))
    predn_ds_tb <- predn_ds_tb %>% dplyr::mutate(`:=`(!!rlang::sym(depnt_var_nm_1L_chr),
        stats::predict(model_mdl, newdata = predn_ds_tb, type = predn_type_1L_chr) %>%
            calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
                tfmn_is_outp_1L_lgl = T)))
    return(predn_ds_tb)
}
make_predrs_for_best_mdls <- function(outp_smry_ls,
                                      old_nms_chr = NULL,
                                      new_nms_chr = NULL){
  ordered_mdl_nms_chr <- outp_smry_ls$mdls_smry_tb %>%
    dplyr::filter(Parameter=="R2") %>%
    dplyr::arrange(dplyr::desc(Estimate)) %>% dplyr::pull(Model)
  ind_predr_mdls_by_mdl_type_ls <- outp_smry_ls$prefd_mdl_types_chr %>%
    purrr::map(~{
    mdl_type_1L_chr <- .x
    paste0(outp_smry_ls$predr_cmprsn_tb$predr_chr,"_1_") %>% purrr::map(~paste0(.x,
                                                                                mdl_type_1L_chr)) %>% purrr::flatten_chr()
  })
  ordered_mdls_by_type_ls <- ind_predr_mdls_by_mdl_type_ls %>%
    purrr::map(~{
      ind_predr_mdls_chr <- .x
      ordered_mdl_nms_chr[ordered_mdl_nms_chr %>% purrr::map_lgl(~.x %in% ind_predr_mdls_chr)]
    })
  predrs_for_best_mdls_chr <- ordered_mdls_by_type_ls %>% purrr::map2_chr(outp_smry_ls$prefd_mdl_types_chr,
                                                                          ~{
                                                                            .x[1] %>% stringr::str_remove(paste0("_1_",.y))
                                                                          })
  if(!is.null(old_nms_chr)){
    predrs_for_best_mdls_chr <- transform_predr_nm_part_of_phrases(predrs_for_best_mdls_chr,
                                                                   old_nms_chr = old_nms_chr,
                                                                   new_nms_chr = new_nms_chr)
  }
  return(predrs_for_best_mdls_chr)
}
make_predr_vals <- function (predr_var_nm_1L_chr, candidate_predrs_lup = NULL)
{
    if (is.null(candidate_predrs_lup)) {
        utils::data("candidate_predrs_lup", envir = environment())
    }
    args_ls <- purrr::map_dbl(names(candidate_predrs_lup)[3:4],
                              ~candidate_predrs_lup %>% ready4fun::get_from_lup_obj(match_value_xx = predr_var_nm_1L_chr,
                                                                                    match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = .x,
                                                                                    evaluate_lgl = F)) %>% as.list()
    predr_vals_dbl <- rlang::exec(.fn = seq, !!!args_ls)
    return(predr_vals_dbl)
}
make_predrs_for_best_mdls <- function(outp_smry_ls,
                                      old_nms_chr = NULL,
                                      new_nms_chr = NULL){
  ordered_mdl_nms_chr <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Parameter=="R2") %>% dplyr::arrange(dplyr::desc(Estimate)) %>% dplyr::pull(Model)
  ind_predr_mdls_by_mdl_type_ls <- outp_smry_ls$prefd_mdl_types_chr %>% purrr::map(~{
    mdl_type_1L_chr <- .x
    paste0(outp_smry_ls$predr_cmprsn_tb$predr_chr,"_1_") %>% purrr::map(~paste0(.x,
                                                                                mdl_type_1L_chr)) %>% purrr::flatten_chr()
  })
  ordered_mdls_by_type_ls <- ind_predr_mdls_by_mdl_type_ls %>%
    purrr::map(~{
      ind_predr_mdls_chr <- .x
      ordered_mdl_nms_chr[ordered_mdl_nms_chr %>% purrr::map_lgl(~.x %in% ind_predr_mdls_chr)]
    })
  predrs_for_best_mdls_chr <- ordered_mdls_by_type_ls %>% purrr::map2_chr(outp_smry_ls$prefd_mdl_types_chr,
                                                                          ~{
                                                                            .x[1] %>% stringr::str_remove(paste0("_1_",.y))
                                                                          })
  if(!is.null(old_nms_chr)){
    predrs_for_best_mdls_chr <- transform_predr_nm_part_of_phrases(predrs_for_best_mdls_chr,
                                                                   old_nms_chr = old_nms_chr,
                                                                   new_nms_chr = new_nms_chr)
  }
  return(predrs_for_best_mdls_chr)
}
make_predr_vars_nms_ls <- function (main_predrs_chr, covars_ls, existing_predrs_ls = NULL)
{
    predr_vars_nms_ls <- covars_ls %>% purrr::map(~{
        covars_chr <- .x
        purrr::map(main_predrs_chr, ~list(c(.x), c(.x, covars_chr) %>% purrr::discard(is.na))) %>%
            purrr::flatten()
    }) %>% purrr::flatten() %>% unique()
    predr_vars_nms_ls <- predr_vars_nms_ls[order(sapply(predr_vars_nms_ls,
                                                        length))]
    if(!is.null(existing_predrs_ls)){
      predr_vars_nms_ls <- predr_vars_nms_ls[predr_vars_nms_ls %>% purrr::map_lgl(~{
        test_chr <- .x
        !any(existing_predrs_ls %>% purrr::map_lgl(~identical(.x,test_chr))
        )})]
    }

    return(predr_vars_nms_ls)
}
make_prefd_mdls_vec <- function (smry_of_sngl_predr_mdls_tb,
                                 choose_from_pfx_chr = c("BET","GLM", "OLS"),
                                 mdl_types_lup = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment(), package =  "TTU")
    ordered_mdl_types_chr <- dplyr::inner_join(smry_of_sngl_predr_mdls_tb %>%
                                                   dplyr::select(Model) %>%
                                                 dplyr::rename(long_name_chr = Model),
                                               mdl_types_lup) %>%
      dplyr::pull(short_name_chr)
    prefd_mdls_chr <- purrr::map_chr(choose_from_pfx_chr, ~ordered_mdl_types_chr[startsWith(ordered_mdl_types_chr,
                                                                                            .x)][1])
    prefd_mdls_chr <- prefd_mdls_chr[order(prefd_mdls_chr %>% purrr::map_int(~which(ordered_mdl_types_chr==.x)))]
    prefd_mdls_chr <- prefd_mdls_chr[-max((1:length(prefd_mdls_chr))[prefd_mdls_chr %>% purrr::map_lgl(~(startsWith(.x,"BET") | startsWith(.x,"GLM")))])]
    return(prefd_mdls_chr)
}
make_prmry_analysis_params_ls <- function(analysis_core_params_ls,
                                          candidate_covar_nms_chr = NA_character_,
                                          ds_tb,
                                          path_params_ls,
                                          maui_params_ls,
                                          prefd_covars_chr = NULL,
                                          prefd_mdl_types_chr = NULL,
                                          raw_ds_tfmn_fn = NULL,
                                          subtitle_1L_chr = "Methods Report 1: Analysis Program (Primary Analysis)",
                                          utl_class_fn_1L_chr = "as.numeric"){
  if(is.na(analysis_core_params_ls$candidate_covar_nms_chr[1]) & !is.na(candidate_covar_nms_chr[1])){
    analysis_core_params_ls$candidate_covar_nms_chr <- candidate_covar_nms_chr
  } # This logic is required to ensure soft deprecated arguments used in manuscript continue to work.
  if(is.null(analysis_core_params_ls$prefd_covars_chr) & !is.null(prefd_covars_chr)){
    analysis_core_params_ls$prefd_covars_chr <- prefd_covars_chr
  }
  if(is.null(analysis_core_params_ls$prefd_mdl_types_chr) & !is.null(prefd_mdl_types_chr)){
    analysis_core_params_ls$prefd_mdl_types_chr <- prefd_mdl_types_chr
  }
  prmry_analysis_params_ls <- list(ds_tb = ds_tb,
                                   raw_ds_tfmn_fn = raw_ds_tfmn_fn,
                                   subtitle_1L_chr = subtitle_1L_chr,
                                   utl_class_fn_1L_chr = utl_class_fn_1L_chr) %>%
    append(analysis_core_params_ls) %>%
    append(path_params_ls[1:2]) %>%
    append(maui_params_ls)
  return(prmry_analysis_params_ls)
}
make_psych_predrs_lup <- function(){
  predictors_lup <- TTU_predictors_lup(make_pt_TTU_predictors_lup(short_name_chr = c("k10_int","psych_well_int"),
                                                                  long_name_chr = c("Kessler Psychological Distress - 10 Item Total Score",
                                                                                    "Overall Wellbeing Measure (Winefield et al. 2012)"),
                                                                  min_val_dbl = c(10,18),
                                                                  max_val_dbl = c(50,90),
                                                                  class_chr = "integer",
                                                                  increment_dbl = 1,
                                                                  class_fn_chr = "integer",
                                                                  mdl_scaling_dbl = 0.01,
                                                                  covariate_lgl = F))
  return(predictors_lup)
}
make_random_forest_text <- function(results_ls,
                                    for_abstract_1L_lgl = F){
  if(for_abstract_1L_lgl){
    text_1L_chr <- ifelse(length(results_ls$ttu_cs_ls$rf_seq_dscdng_chr) > 1,
                          paste0(ifelse(results_ls$ttu_cs_ls$mdl_predrs_and_rf_seqs_cmprsn_1L_chr == "is consistent",
                                        " and the random forest model",
                                        paste0(", while ",
                                               results_ls$ttu_cs_ls$rf_seq_dscdng_chr[1],
                                               " was the most \'important\' predictor in the random forest model"))),
                          "")
  }else{
    text_1L_chr <- ifelse(length(results_ls$ttu_cs_ls$rf_seq_dscdng_chr) > 1,
                          paste0("This ",
                                 results_ls$ttu_cs_ls$mdl_predrs_and_rf_seqs_cmprsn_1L_chr,
                                 " with the random forest model in which ",
                                 results_ls$ttu_cs_ls$rf_seq_dscdng_chr[1],
                                 " was found to be the most \'important\' predictor"),
                          "")
  }

  return(text_1L_chr)
}
make_ranked_predrs_ls <- function(descv_tbls_ls,
                                  old_nms_chr = NULL,
                                  new_nms_chr = NULL){
  unranked_predrs_chr <- rownames(descv_tbls_ls[["bl_cors_tb"]])[-1]
  if(!is.null(old_nms_chr))
    unranked_predrs_chr <- unranked_predrs_chr %>%
    transform_predr_nm_part_of_phrases(old_nms_chr = old_nms_chr,
                                       new_nms_chr = new_nms_chr)
  ranks_dbl <- descv_tbls_ls[["bl_cors_tb"]][2:nrow(descv_tbls_ls[["bl_cors_tb"]]),1] %>%
    purrr::map_dbl(~{
      vec_dbl <- regmatches(.x, gregexpr("[[:digit:]]+", .x)) %>% unlist() %>% as.numeric()
      vec_dbl[1]+vec_dbl[2]/100
    }) %>% rank()
  unique_ranks_dbl <- unique(ranks_dbl) %>% sort(decreasing = T)
  ranked_predrs_chr <- purrr::map(unique_ranks_dbl, ~ unranked_predrs_chr[ranks_dbl == .x]) %>% purrr::flatten_chr()
  ranked_predrs_ls <- list(unranked_predrs_chr = unranked_predrs_chr,
                           ranked_predrs_chr = ranked_predrs_chr)
  return(ranked_predrs_ls)

}
make_results_ls <- function(spine_of_results_ls = NULL,
                            abstract_args_ls = NULL,
                            dv_ds_nm_and_url_chr = NULL,
                            output_format_ls = NULL,
                            params_ls_ls = NULL,
                            path_params_ls = NULL,
                            study_descs_ls = NULL,
                            fn_ls = NULL,
                            include_idx_int = NULL,
                            var_nm_change_lup = NULL,
                            ctgl_vars_regrouping_ls = NULL,
                            sig_covars_some_predrs_mdls_tb = NULL,
                            sig_thresh_covars_1L_chr = NULL,
                            version_1L_chr = NULL){
  if(is.null(spine_of_results_ls)){
    spine_of_results_ls <- make_results_ls_spine(output_format_ls = output_format_ls,
                                                 params_ls_ls = params_ls_ls,
                                                 path_params_ls = path_params_ls,
                                                 study_descs_ls = study_descs_ls,
                                                 fn_ls = fn_ls,
                                                 include_idx_int = include_idx_int,
                                                 var_nm_change_lup = var_nm_change_lup)
  }
  mdls_smry_tbls_ls <- make_mdls_smry_tbls_ls(spine_of_results_ls$outp_smry_ls,
                                              nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int)
  covars_mdls_ls <- make_mdls_ls(spine_of_results_ls$outp_smry_ls,
                                 mdls_tb = mdls_smry_tbls_ls$covar_mdls_tb)
  descv_tbls_ls <- paste0(spine_of_results_ls$output_data_dir_1L_chr,"/",spine_of_results_ls$outp_smry_ls$file_paths_chr[spine_of_results_ls$outp_smry_ls$file_paths_chr %>% purrr::map_lgl(~stringr::str_detect(.x,"descv_tbls_ls.RDS"))]) %>% readRDS()
  composite_plt <- make_cmpst_sctr_and_dnsty_plt(spine_of_results_ls$outp_smry_ls,
                                                 output_data_dir_1L_chr = spine_of_results_ls$output_data_dir_1L_chr,
                                                 predr_var_nms_chr = spine_of_results_ls$outp_smry_ls$predr_vars_nms_ls[[1]])
  cowplot::save_plot(paste0(spine_of_results_ls$output_data_dir_1L_chr,"/dens_and_sctr.png"), composite_plt, base_height = 20)
  ttu_cs_ls <- make_ttu_cs_ls(spine_of_results_ls$outp_smry_ls,
                             sig_covars_some_predrs_mdls_tb = sig_covars_some_predrs_mdls_tb,
                             sig_thresh_covars_1L_chr = sig_thresh_covars_1L_chr)
  mdl_type_descs_chr <- mdls_smry_tbls_ls$prefd_predr_mdl_smry_tb$Model %>%
    purrr::map_chr(~get_mdl_type_from_nm(.x)) %>%
    unique() %>%
    purrr::map_chr(~ready4fun::get_from_lup_obj(spine_of_results_ls$outp_smry_ls$mdl_types_lup,
                                                match_value_xx = .x,
                                                match_var_nm_1L_chr = "short_name_chr",
                                                target_var_nm_1L_chr = "long_name_chr",
                                                evaluate_lgl = F))
  ttu_cs_ls$rf_seq_dscdng_chr <- ttu_cs_ls$rf_seq_dscdng_chr %>%
    purrr::map_chr(~ifelse(.x %in% spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                                                         .x %>%
                                                   ready4fun::get_from_lup_obj(data_lookup_tb = spine_of_results_ls$var_nm_change_lup,
                                                                               match_var_nm_1L_chr = "old_nms_chr",
                                                                               target_var_nm_1L_chr = "new_nms_chr",
                                                                               evaluate_lgl = F),
                                                   .x))
  ttu_cs_ls$cs_mdls_predrs_seq_dscdng_chr <- ttu_cs_ls$cs_mdls_predrs_seq_dscdng_chr  %>%
    purrr::map_chr(~ifelse(.x %in% spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                           .x %>%
                             ready4fun::get_from_lup_obj(data_lookup_tb = spine_of_results_ls$var_nm_change_lup,
                                                         match_var_nm_1L_chr = "old_nms_chr",
                                                         target_var_nm_1L_chr = "new_nms_chr",
                                                         evaluate_lgl = F),
                           .x))
  ttu_lngl_ls = list(best_mdls_tb = tibble::tibble(model_type = mdl_type_descs_chr %>%
                                                     purrr::map_chr(~ifelse(startsWith(.x,"Ordinary Least Squares"),"LLM","GLMM")),
                                                   link_and_tfmn_chr = mdl_type_descs_chr %>%
                                                     purrr::map_chr(~ifelse(startsWith(.x,"Ordinary Least Squares"),
                                                                            stringr::str_remove(.x,"Ordinary Least Squares ") %>%
                                                                              stringr::str_sub(start = 2, end = -2) %>%
                                                                              tolower(),
                                                                            stringr::str_remove(.x,"Generalised Linear Mixed Model with ") %>%
                                                                              stringr::str_remove("Beta Regression Model with "))),
                                                   name_chr = make_predrs_for_best_mdls(spine_of_results_ls$outp_smry_ls,
                                                                                        old_nms_chr = spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                                                                                        new_nms_chr = spine_of_results_ls$var_nm_change_lup$new_nms_chr),
                                                   r2_dbl = mdls_smry_tbls_ls$prefd_predr_mdl_smry_tb %>%
                                                     dplyr::filter(Parameter == "R2") %>%
                                                     dplyr::pull(Estimate)),
                     cs_ts_ratios_tb = spine_of_results_ls$cs_ts_ratios_tb,
                     incld_covars_chr = spine_of_results_ls$outp_smry_ls$prefd_covars_chr)
  results_ls <- list(abstract_args_ls = abstract_args_ls,
                     candidate_covars_ls = spine_of_results_ls$candidate_covars_ls,
                     candidate_predrs_chr = spine_of_results_ls$candidate_predrs_chr,
                     cohort_ls = make_cohort_ls(descv_tbls_ls,
                                                ctgl_vars_regrouping_ls = ctgl_vars_regrouping_ls,
                                                nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int),
                     dv_ds_nm_and_url_chr = dv_ds_nm_and_url_chr,
                     header_yaml_args_ls = input_params_ls$header_yaml_args_ls,
                     hlth_utl_and_predrs_ls = make_hlth_utl_and_predrs_ls(spine_of_results_ls$outp_smry_ls,
                                                                          descv_tbls_ls = descv_tbls_ls,
                                                                          nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int,
                                                                          old_nms_chr = spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                                                                          new_nms_chr = spine_of_results_ls$var_nm_change_lup$new_nms_chr),
                     mdl_coef_ratios_ls = spine_of_results_ls$mdl_coef_ratios_ls,
                     mdl_ingredients_ls = spine_of_results_ls$mdl_ingredients_ls,
                     mdls_with_signft_covars_ls = spine_of_results_ls$mdls_with_signft_covars_ls,
                     output_format_ls = input_params_ls$output_format_ls,
                     path_params_ls = input_params_ls$path_params_ls,
                     paths_to_figs_ls = make_paths_to_ss_plts_ls(spine_of_results_ls$output_data_dir_1L_chr,
                                                                 outp_smry_ls = spine_of_results_ls$outp_smry_ls),
                     predr_var_nms_chr = spine_of_results_ls$outp_smry_ls$predr_vars_nms_ls[[1]] %>%
                       purrr::map_chr(~ifelse(.x %in% spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                                              .x %>%
                                                ready4fun::get_from_lup_obj(data_lookup_tb = spine_of_results_ls$var_nm_change_lup,
                                                                            match_var_nm_1L_chr = "old_nms_chr",
                                                                            target_var_nm_1L_chr = "new_nms_chr",
                                                                            evaluate_lgl = F),
                                              .x)),
                     r_version_1L_chr = paste0(spine_of_results_ls$outp_smry_ls$session_data_ls$R.version$major,
                                               ".",
                                               spine_of_results_ls$outp_smry_ls$session_data_ls$R.version$minor),
                     study_descs_ls = spine_of_results_ls$study_descs_ls,
                     tables_ls = make_ss_tbls_ls(spine_of_results_ls$outp_smry_ls,
                                                 mdls_smry_tbls_ls = mdls_smry_tbls_ls,
                                                 covars_mdls_ls = covars_mdls_ls,
                                                 descv_tbls_ls = descv_tbls_ls,
                                                 nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int),
                     ttu_cs_ls = ttu_cs_ls,
                     ttu_lngl_ls = ttu_lngl_ls,
                     ttu_version_1L_chr = spine_of_results_ls$outp_smry_ls$session_data_ls$otherPkgs$TTU$Version,
                     var_nm_change_lup = spine_of_results_ls$var_nm_change_lup,
                     version_1L_chr = version_1L_chr)
  return(results_ls)
}
make_results_ls_spine <-  function(output_format_ls = NULL,
                                   params_ls_ls = NULL,
                                   path_params_ls = NULL,
                                   study_descs_ls,
                                   fn_ls = NULL,
                                   include_idx_int = NULL,
                                   nbr_of_digits_1L_int = NULL,
                                   output_data_dir_1L_chr = NULL,
                                   var_nm_change_lup = NULL){
  if(!missing("output_data_dir_1L_chr")){
    warning("The output_data_dir_1L_chr argument is soft deprecated.")
  }else{
    output_data_dir_1L_chr <- path_params_ls$paths_ls$output_data_dir_1L_chr
  }
  if(!missing("nbr_of_digits_1L_int")){
    warning("The nbr_of_digits_1L_int argument is soft deprecated.")
  }else{
    nbr_of_digits_1L_int <- output_format_ls$manuscript_digits_1L_int
  }
  if(is.null(var_nm_change_lup)){
    var_nm_change_lup <- list(old_nms_chr = NULL,
                              new_nms_chr = NULL)
  }
  outp_smry_ls <- readRDS(paste0(output_data_dir_1L_chr,"/I_ALL_OUTPUT_.RDS"))
  mdl_ingredients_ls <- readRDS(paste0(output_data_dir_1L_chr,
                                       "/G_Shareable/Ingredients/mdl_ingredients.RDS"))
  if(!is.null(params_ls_ls)){
    covar_ctgs_chr <- params_ls_ls$params_ls$candidate_covar_nms_chr %>%
      purrr::map_chr(~ready4fun::get_from_lup_obj(params_ls_ls$params_ls$ds_descvs_ls$dictionary_tb,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "var_nm_chr",
                                                  target_var_nm_1L_chr = "var_ctg_chr",
                                                  evaluate_lgl = F)) %>%
      unique()
    candidate_covars_ls <- covar_ctgs_chr  %>%
      purrr::map(~{
        var_desc_chr <- ready4fun::get_from_lup_obj(params_ls_ls$params_ls$ds_descvs_ls$dictionary_tb,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "var_ctg_chr",
                                                  target_var_nm_1L_chr = "var_desc_chr",
                                                  evaluate_lgl = F)
        class(var_desc_chr) <- setdiff(class(var_desc_chr), "labelled")
        attr(var_desc_chr, "label") <- NULL
        var_nm_chr <- ready4fun::get_from_lup_obj(params_ls_ls$params_ls$ds_descvs_ls$dictionary_tb,
                                                    match_value_xx = .x,
                                                    match_var_nm_1L_chr = "var_ctg_chr",
                                                    target_var_nm_1L_chr = "var_nm_chr",
                                                    evaluate_lgl = F)
        var_desc_chr[var_nm_chr %in% params_ls_ls$params_ls$candidate_covar_nms_chr]
      }
        ) %>%
      stats::setNames(covar_ctgs_chr)
  }else{
    candidate_covars_ls <- NULL
  }
  study_descs_ls$predr_ctgs_ls <- make_predr_ctgs_ls(outp_smry_ls,
                                                     include_idx_int = include_idx_int)
  mdl_coef_ratios_ls <- make_mdl_coef_ratio_ls(mdl_ingredients_ls,
                                               predr_ctgs_ls = study_descs_ls$predr_ctgs_ls)
  mdls_smry_tbls_ls <- make_mdls_smry_tbls_ls(outp_smry_ls,
                                                   nbr_of_digits_1L_int = nbr_of_digits_1L_int)
  covars_mdls_ls <- make_mdls_ls(outp_smry_ls,
                                      mdls_tb = mdls_smry_tbls_ls$covar_mdls_tb)
  cs_ts_ratios_tb <- make_cs_ts_ratios_tb(predr_ctgs_ls = study_descs_ls$predr_ctgs_ls,
                                          mdl_coef_ratios_ls = mdl_coef_ratios_ls,
                                          fn_ls = fn_ls,
                                          nbr_of_digits_1L_int = nbr_of_digits_1L_int)
  spine_of_results_ls <- list(candidate_covars_ls = candidate_covars_ls,
                              candidate_predrs_chr = params_ls_ls$params_ls$ds_descvs_ls$candidate_predrs_chr,
                              cs_ts_ratios_tb = cs_ts_ratios_tb,
                              mdls_with_signft_covars_ls = get_mdls_with_signft_covars(outp_smry_ls,
                                                                                        params_ls_ls = params_ls_ls),
                              outp_smry_ls = outp_smry_ls,
                              output_data_dir_1L_chr = output_data_dir_1L_chr,
                              mdl_coef_ratios_ls = mdl_coef_ratios_ls,
                              mdl_ingredients_ls = mdl_ingredients_ls,
                              nbr_of_digits_1L_int = nbr_of_digits_1L_int,
                              study_descs_ls = study_descs_ls,
                              var_nm_change_lup = var_nm_change_lup)
  return(spine_of_results_ls)
}
make_scaling_text <- function(results_ls,
                              table_1L_chr = "cfscl"){
  if(startsWith(table_1L_chr,"cfscl")){
    table_df <- results_ls$tables_ls$ind_preds_coefs_tbl
  }else{
    if(startsWith(table_1L_chr,"coefscovarstype")){
      table_df <- results_ls$tables_ls %>%
        purrr::pluck(paste0("mdl_type_",
             table_1L_chr %>% stringr::str_remove("coefscovarstype"),
             "_covar_mdls_tb"))
    }
  }
  predrs_chr <- table_df$Parameter %>% setdiff(c("SD (Intercept)","Intercept")) %>% stringr::str_replace_all(" model","") %>% stringr::str_replace_all(" baseline","") %>% stringr::str_replace_all(" change","") %>% unique()
  predrs_lup <- results_ls$mdl_ingredients_ls$predictors_lup %>%
    dplyr::filter(short_name_chr %in% predrs_chr)
  scaling_dbl <- predrs_lup$mdl_scaling_dbl %>% unique()
  text_1L_chr <- ifelse(all(scaling_dbl==1),
                        "",
                        paste0("Note: ",

                               scaling_dbl %>%
                                 purrr::map_chr(~ {
                                   scaled_predrs_chr <- predrs_lup %>%
                                     dplyr::filter(mdl_scaling_dbl == .x) %>%
                                     dplyr::pull(short_name_chr) %>%
                                     sort()
                                   ifelse(.x ==1,
                                          "",
                                          paste0("The ",
                                                 scaled_predrs_chr %>%
                                                   paste0(collapse = ", ") %>%
                                                   stringi::stri_replace_last(fixed = ",", " and"),
                                                 " parameter",
                                                 ifelse(length(scaled_predrs_chr) == 1,
                                                        " was",
                                                        "s were"),
                                                 " first multiplied by ",
                                                 .x,
                                                 "."))
                                 }
                                 ) %>%
                                 paste0(collapse = " ")))
  return(text_1L_chr)
}
make_scndry_anlys_params <- function(scndry_anlys_params_ls = NULL,
                                     candidate_covar_nms_chr = NULL,
                                     candidate_predrs_chr = NULL,
                                     predictors_lup = NULL,
                                     prefd_covars_chr = NA_character_
){
  new_params_ls <- list(candidate_covar_nms_chr = candidate_covar_nms_chr,
                        candidate_predrs_chr = candidate_predrs_chr,
                        predictors_lup = predictors_lup,
                        prefd_covars_chr = prefd_covars_chr)
  if(!is.null(scndry_anlys_params_ls)){
    new_params_ls <- append(scndry_anlys_params_ls, list(new_params_ls)) %>%
      stats::setNames(paste0("secondary_",1:(length(scndry_anlys_params_ls) + 1)))
  }else{
    new_params_ls <- list(secondary_1 = new_params_ls)
  }
  return(new_params_ls)
}
make_scndry_anlys_text <- function(results_ls){
  text_1L_chr <- ifelse(get_nbr_of_scndry_analyses(results_ls,
                                                   as_words_1L_lgl = F) == 0,
                        "",
                        paste0(get_nbr_of_scndry_analyses(results_ls) ,
                               " secondary analys",
                               ifelse(get_nbr_of_scndry_analyses(results_ls, as_words_1L_lgl = F) == 1,
                                      "is was ","es were "),
                               "undertaken. ",
                               get_scndry_anlys_descs(results_ls)))
  return(text_1L_chr)
}
make_selected_mdl_text <- function(results_ls,
                                   for_abstract_1L_lgl = F){
  length_1L_int <- length(results_ls$ttu_cs_ls$selected_mdls_chr)
  text_1L_chr <- paste0(ifelse((length_1L_int == 2 & !for_abstract_1L_lgl),"Both ",""),
                        results_ls$ttu_cs_ls$selected_mdls_chr %>%
                          paste0(collapse = ", ") %>%
                          stringi::stri_replace_last(fixed = ",", " and"),
                        ifelse(length_1L_int >1," were"," was"),
                        ifelse(for_abstract_1L_lgl,
                               paste0(" the best peforming model",ifelse(length_1L_int >1,"s.",".")),
                               " selected for further evaluation."))
  return(text_1L_chr)
}
make_shareable_mdl <- function (fake_ds_tb, mdl_smry_tb, depnt_var_nm_1L_chr = "utl_total_w",
    id_var_nm_1L_chr = "fkClientID", tfmn_1L_chr = "CLL", mdl_type_1L_chr = "OLS_CLL",
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NA_character_,
    seed_1L_int = 12345L)
{
    if (is.null(mdl_types_lup))
        utils::data(mdl_types_lup, envir = environment())
  if (is.na(tfmn_1L_chr))
    tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
                                               match_value_xx = mdl_type_1L_chr,
                                               match_var_nm_1L_chr = "short_name_chr",
                                               target_var_nm_1L_chr = "tfmn_chr",
                                               evaluate_lgl = F)
  predr_var_nms_chr <- mdl_smry_tb$Parameter[!mdl_smry_tb$Parameter %in% c("SD (Intercept)", "Intercept",
                                                      "R2", "RMSE","Sigma")] %>%
    stringi::stri_replace_last_fixed(" baseline","_baseline") %>%
    stringi::stri_replace_last_fixed(" change","_change")
    tfmd_depnt_var_nm_1L_chr <- transform_depnt_var_nm(depnt_var_nm_1L_chr,
                                                      tfmn_1L_chr = tfmn_1L_chr)
    if (length(predr_var_nms_chr) > 1) {
        covar_var_nms_chr <- predr_var_nms_chr[2:length(predr_var_nms_chr)]
    }else{
        covar_var_nms_chr <- NA_character_
    }
    model_mdl <- make_mdl(fake_ds_tb %>%
                            dplyr::select(tidyselect::all_of(c(id_var_nm_1L_chr,
                                                               tfmd_depnt_var_nm_1L_chr,
                                                               predr_var_nms_chr))),
                          depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        predr_var_nm_1L_chr = predr_var_nms_chr[1], covar_var_nms_chr = covar_var_nms_chr,
        tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr,
        mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr,
        start_1L_chr = start_1L_chr)
    if(ready4fun::get_from_lup_obj(mdl_types_lup,
                                   match_value_xx = mdl_type_1L_chr,
                                   match_var_nm_1L_chr = "short_name_chr",
                                   target_var_nm_1L_chr = "fn_chr",
                                   evaluate_lgl = F) == "betareg::betareg"){
      model_coeffs_dbl <- model_mdl$coefficients$mean

    }else{
      model_coeffs_dbl <- model_mdl$coefficients
    }
    par_nms_chr <- model_coeffs_dbl %>% names()
    mdl_smry_tb <- mdl_smry_tb %>% dplyr::mutate(Parameter = dplyr::case_when(Parameter ==
        "Intercept" ~ "(Intercept)", TRUE ~ purrr::map_chr(Parameter,
        ~stringr::str_replace_all(.x, " ", "_")))) %>% dplyr::filter(Parameter %in%
        par_nms_chr) %>% dplyr::slice(match(par_nms_chr, Parameter))
    assertthat::assert_that(all(par_nms_chr == mdl_smry_tb$Parameter),
        msg = "Parameter names mismatch between data and model summary table")
    model_coeffs_dbl <- mdl_smry_tb$Estimate
    names(model_coeffs_dbl) <- par_nms_chr
    if(ready4fun::get_from_lup_obj(mdl_types_lup,
                                   match_value_xx = mdl_type_1L_chr,
                                   match_var_nm_1L_chr = "short_name_chr",
                                   target_var_nm_1L_chr = "fn_chr",
                                   evaluate_lgl = F) == "betareg::betareg"){
      model_mdl$coefficients$mean <- model_coeffs_dbl

    }else{
      model_mdl$coefficients <- model_coeffs_dbl
    }
    return(model_mdl)
}
make_sngl_mdl_smry_tb <- function(mdls_tb,
                                  mdl_nm_1L_chr,
                                  mdl_type_1L_chr,
                                  add_mdl_nm_sfx_1L_lgl = T){
  new_tb <- mdls_tb %>%
    dplyr::filter(Model == mdl_nm_1L_chr) %>%
    tibble::add_case(Parameter = mdl_nm_1L_chr, .before = 1) %>%
    dplyr::select(-Model)
  mdl_nm_sfx_1L_chr <- ifelse(add_mdl_nm_sfx_1L_lgl,
                              paste0("_",mdl_type_1L_chr),
                              "")
  new_tb <- new_tb %>%
    dplyr::mutate(!!rlang::sym(paste0("R2",mdl_nm_sfx_1L_chr)) := dplyr::case_when(Parameter == mdl_nm_1L_chr ~ ready4fun::get_from_lup_obj(new_tb,
                                                                                                                                            match_value_xx = "R2",
                                                                                                                                            match_var_nm_1L_chr = "Parameter",
                                                                                                                                            target_var_nm_1L_chr = "Estimate",
                                                                                                                                            evaluate_lgl = F),
                                                                                   T ~ NA_character_),
                  !!rlang::sym(paste0("Sigma",mdl_nm_sfx_1L_chr)) := dplyr::case_when(Parameter == mdl_nm_1L_chr ~ ready4fun::get_from_lup_obj(new_tb,
                                                                                                                                               match_value_xx = "Sigma",
                                                                                                                                               match_var_nm_1L_chr = "Parameter",
                                                                                                                                               target_var_nm_1L_chr = "Estimate",
                                                                                                                                               evaluate_lgl = F),
                                                                                      T ~ NA_character_)) %>%
    dplyr::filter(!Parameter %in% c("R2","RMSE","Sigma")) %>%
    dplyr::rename_with(~paste0(.x,mdl_nm_sfx_1L_chr),.cols = c("Parameter","Estimate","SE")) %>%
    dplyr::rename_with(~paste0("CI",mdl_nm_sfx_1L_chr),.cols = c("95% CI")) %>%
    dplyr::mutate(!!rlang::sym(paste0("Parameter",mdl_nm_sfx_1L_chr)) := !!rlang::sym(paste0("Parameter",mdl_nm_sfx_1L_chr)) %>% purrr::map_chr(~stringr::str_replace(.x,
                                                                                                                                                                      paste0("_1_",mdl_type_1L_chr), " model")))
  rownames(new_tb)  <- NULL
  return(new_tb)
}
make_smry_of_brm_mdl <- function (mdl_ls,
                                  data_tb,
                                  depnt_var_nm_1L_chr = "utl_total_w",
                                  predr_vars_nms_chr,
                                  mdl_nm_1L_chr = NA_character_,
                                  seed_1L_dbl = 23456,
                                  tfmn_1L_chr) {
    if (is.na(mdl_nm_1L_chr))
        mdl_nm_1L_chr <- predr_vars_nms_chr[1]
    set.seed(seed_1L_dbl)
    predictions <- stats::predict(mdl_ls, summary = F) %>%
      calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
                              tfmn_is_outp_1L_lgl = T)
    sd_intcpt_df <- summary(mdl_ls, digits = 4)$random[[1]]
    sd_intcpt_df <- sd_intcpt_df[1:nrow(sd_intcpt_df), 1:4]
    coef <- summary(mdl_ls, digits = 4)$fixed
    coef <- coef[1:nrow(coef), 1:4]
    R2 <- brms::bayes_R2(mdl_ls)
    RMSE <- psych::describe(apply(predictions, 1, calculate_rmse, y_dbl = data_tb %>%
        dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))), quant = c(0.25,
        0.75), skew = F, ranges = F)
    RMSE <- cbind(RMSE$mean, RMSE$sd, RMSE$Q0.25, RMSE$Q0.75) %>%
        as.vector()
    Sigma <- summary(mdl_ls, digits = 4)$spec_par[1:4]
    smry_of_brm_mdl_tb <- data.frame(round(rbind(sd_intcpt_df, coef, R2, RMSE,
        Sigma), 3)) %>% dplyr::mutate(Parameter = c("SD (Intercept)","Intercept",
                                                    purrr::map(predr_vars_nms_chr,
                                                               ~paste0(.x, c(" baseline",
                                                                             " change"))) %>%
                                                      purrr::flatten_chr(),
                                                    "R2", "RMSE", "Sigma"),
                                      Model = mdl_nm_1L_chr) %>%
      dplyr::mutate(`95% CI` = paste(l.95..CI,
                                     ",",
                                     u.95..CI)) %>%
      dplyr::rename(SE = Est.Error) %>%
        dplyr::select(Model, Parameter, Estimate, SE, `95% CI`)
    return(smry_of_brm_mdl_tb)
}
make_smry_of_mdl_outp <- function (data_tb,
                                   model_mdl = NULL, # REMOVE
                                   folds_1L_int = 10,
                                   depnt_var_nm_1L_chr = "utl_total_w",
                                   start_1L_chr = NULL,
                                   tfmn_1L_chr = "NTF",
                                   predr_var_nm_1L_chr,
                                   covar_var_nms_chr = NA_character_,
                                   mdl_type_1L_chr = "OLS_NTF",
                                   mdl_types_lup = NULL,
                                   predn_type_1L_chr = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- data_tb %>% dplyr::filter(!is.na(!!rlang::sym(predr_var_nm_1L_chr)))
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    mdl_desc_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
        target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F)
    folds_ls <- make_folds_ls(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        folds_1L_int = folds_1L_int)
    control_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
                                                  match_var_nm_1L_chr = "short_name_chr",
                                                  match_value_xx = mdl_type_1L_chr,
                                                  target_var_nm_1L_chr = "control_chr",
                                                  evaluate_lgl = F)
    smry_of_one_predr_mdl_tb <- purrr::map_dfr(folds_ls, ~{
        model_mdl <- make_mdl(data_tb[-.x,], depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
            start_1L_chr = start_1L_chr, tfmn_1L_chr = tfmn_1L_chr,
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr,
            mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr)
        pred_old_dbl <- stats::predict(model_mdl, type = predn_type_1L_chr) %>%
          calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, ## TEST
                                  tfmn_is_outp_1L_lgl = T)
        pred_new_dbl <- stats::predict(model_mdl, newdata = data_tb[.x, ],
            type = predn_type_1L_chr) %>% calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
            tfmn_is_outp_1L_lgl = T)
        tibble::tibble(Rsquared = caret::R2(pred_old_dbl, data_tb[-.x,
            ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
            form = "traditional"), RMSE = caret::RMSE(pred_old_dbl,
            data_tb[-.x, ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
            MAE = caret::MAE(pred_old_dbl, data_tb[-.x, ] %>%
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
            RsquaredP = caret::R2(pred_new_dbl, data_tb[.x, ] %>%
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
                form = "traditional"), RMSEP = caret::RMSE(pred_new_dbl,
                data_tb[.x, ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
            MAEP = caret::MAE(pred_new_dbl, data_tb[.x, ] %>%
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))))
    }) %>% dplyr::summarise_all(mean) %>% dplyr::mutate(Model = mdl_desc_1L_chr) %>%
        dplyr::select(Model, dplyr::everything())
    return(smry_of_one_predr_mdl_tb)
}
make_smry_of_ts_mdl_outp <- function (data_tb,
                                      predr_vars_nms_chr, mdl_nm_1L_chr, path_to_write_to_1L_chr = NA_character_,
                                      depnt_var_nm_1L_chr = "utl_total_w", id_var_nm_1L_chr = "fkClientID",
                                      round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", predictors_lup, utl_min_val_1L_dbl = -1,
                                      backend_1L_chr = getOption("brms.backend", "rstan"), iters_1L_int = 4000L, mdl_types_lup,
                                      seed_1L_int = 1000L, prior_ls = NULL, control_ls = NULL)
{
  scaling_fctr_dbl <- predr_vars_nms_chr %>% purrr::map_dbl(~
      ifelse(.x %in% predictors_lup$short_name_chr,
             ready4fun::get_from_lup_obj(predictors_lup,
                                         target_var_nm_1L_chr = "mdl_scaling_dbl",
                                         match_value_xx = .x,
                                         match_var_nm_1L_chr = "short_name_chr",
                                         evaluate_lgl = F),
             1))
  mdl_type_1L_chr <- mdl_nm_1L_chr %>%
    stringr::str_remove(paste0(predr_vars_nms_chr[1], "_", ifelse(is.na(predr_vars_nms_chr[2]), "", paste0(predr_vars_nms_chr[2],
                                                                                                    "_"))))
  mdl_type_1L_chr <- stringr::str_sub(mdl_type_1L_chr, start = 1 + (mdl_type_1L_chr %>% stringi::stri_locate_first_fixed("_"))[1,2] %>% as.vector())
  tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
                                             target_var_nm_1L_chr = "tfmn_chr",
                                             match_value_xx = mdl_type_1L_chr,
                                             match_var_nm_1L_chr = "short_name_chr",
                                             evaluate_lgl = F)
  tfd_data_tb <- transform_tb_to_mdl_inp(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr,
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr,
        scaling_fctr_dbl = scaling_fctr_dbl,
        tfmn_1L_chr = tfmn_1L_chr)
    tfd_depnt_var_nm_1L_chr <- transform_depnt_var_nm(depnt_var_nm_1L_chr,
                                                      tfmn_1L_chr = tfmn_1L_chr)
    family_fn_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
                                                   match_var_nm_1L_chr = "short_name_chr",
                                                   match_value_xx = mdl_type_1L_chr,
                                                   target_var_nm_1L_chr = "family_chr",
                                                   evaluate_lgl = F)
    family_fn_1L_chr <- ifelse(is.na(family_fn_1L_chr),
                               ifelse(startsWith(mdl_type_1L_chr,"BET"),
                                      paste0("brms::Beta(link = \"",
                                             get_link_from_tfmn(stringr::str_sub(mdl_type_1L_chr,
                                                                          start = -3)),
                                             "\")"),
                                      "gaussian(identity)"),
                               family_fn_1L_chr)
    args_ls <- list(data_tb = tfd_data_tb, depnt_var_nm_1L_chr = tfd_depnt_var_nm_1L_chr,
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, iters_1L_int = iters_1L_int,
        backend_1L_chr = backend_1L_chr,
        family_fn_1L_chr = family_fn_1L_chr,
        seed_1L_int = seed_1L_int,
        prior_ls = prior_ls, control_ls = control_ls)
    # if(startsWith(mdl_type_1L_chr, "GLM_BNL")){
    # WRITE FN
    # }else{
      mdl_ls <- rlang::exec(fit_ts_model_with_brm, !!!args_ls)
    # }
    smry_of_ts_mdl_ls <- list(smry_of_ts_mdl_tb = make_smry_of_brm_mdl(mdl_ls,
        data_tb = tfd_data_tb, depnt_var_nm_1L_chr = tfd_depnt_var_nm_1L_chr,
        predr_vars_nms_chr = predr_vars_nms_chr,
        mdl_nm_1L_chr = mdl_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr))
    if (!is.na(path_to_write_to_1L_chr)) {
        smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr <- paste0(path_to_write_to_1L_chr,
            "/", mdl_nm_1L_chr, ".RDS")
        if (file.exists(smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr))
            file.remove(smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)
        saveRDS(mdl_ls, smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)
        smry_of_ts_mdl_ls$paths_to_mdl_plts_chr <- write_ts_mdl_plts(mdl_ls,
                                                                        tfd_data_tb = tfd_data_tb,
                                                                        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                                        mdl_nm_1L_chr = mdl_nm_1L_chr,
                                                                        path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                                                        round_var_nm_1L_chr = round_var_nm_1L_chr,
                                                                        tfmn_1L_chr = tfmn_1L_chr,
                                                                        utl_min_val_1L_dbl = utl_min_val_1L_dbl)
    }
    return(smry_of_ts_mdl_ls)
}
make_ss_tbls_ls <- function(outp_smry_ls,
                            mdls_smry_tbls_ls,
                            covars_mdls_ls,
                            descv_tbls_ls,
                            nbr_of_digits_1L_int = 2L){
  mdl_types_tables_ls <- purrr::map(1:length(outp_smry_ls$prefd_mdl_types_chr),
                                    ~ make_mdl_type_smry_tbl(mdls_tb = mdls_smry_tbls_ls$covar_mdls_tb,
                                                             mdl_nms_chr = covars_mdls_ls[[.x]],
                                                             mdl_type_1L_chr = outp_smry_ls$prefd_mdl_types_chr[.x],
                                                             add_mdl_nm_sfx_1L_lgl = F)) %>%
    stats::setNames(1:length(outp_smry_ls$prefd_mdl_types_chr) %>%
                      purrr::map_chr(~paste0("mdl_type_",.x,"_covar_mdls_tb")))
  ss_tbls_ls <- append(mdl_types_tables_ls,
                       list(ind_preds_coefs_tbl = make_all_mdl_types_smry_tbl(outp_smry_ls,
                                                                              mdls_tb = mdls_smry_tbls_ls$indpt_predrs_mdls_tb),
                     participant_descs = descv_tbls_ls$cohort_desc_tb,
                     pred_dist_and_cors = descv_tbls_ls$predr_pars_and_cors_tb,
                     tenf_prefd_mdl_tb =  outp_smry_ls[["smry_of_mdl_sngl_predrs_tb"]] %>%
                       tibble::as_tibble() %>%
                       dplyr::mutate(dplyr::across(where(is.numeric), ~ .x %>% purrr::map_dbl(~min(max(.x,-1.1),1.1)))) %>%
                       transform_tbl_to_rnd_vars(nbr_of_digits_1L_int = nbr_of_digits_1L_int) %>%
                       dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ .x %>%
                                                     stringr::str_replace_all("-1.10","< -1.00") %>%
                                                     stringr::str_replace_all("1.10","> 1.00"))),
                     tenf_sngl_predr_tb = make_tfd_sngl_predr_mdls_tb(outp_smry_ls,
                                                             nbr_of_digits_1L_int = nbr_of_digits_1L_int)))
  return(ss_tbls_ls)
}
make_study_descs_ls <- function(input_params_ls = NULL,
                                time_btwn_bl_and_fup_1L_chr,
                                background_1L_chr = "",
                                coi_1L_chr = "None declared.",
                                conclusion_1L_chr = "",
                                ethics_1L_chr = NULL,
                                funding_1L_chr = NULL,
                                health_utl_nm_1L_chr = NULL,
                                params_ls_ls = NULL,
                                predr_ctgs_ls = NULL,
                                sample_desc_1L_chr = NULL,
                                var_nm_change_lup = NULL){
  if(!missing("health_utl_nm_1L_chr")){
    warning("The health_utl_nm_1L_chr argument is now soft deprecated. Passing a valid value to params_ls_ls will ensure the names of the MAUI instrument are included in function operations.")
  }else{
    health_utl_nm_1L_chr <- input_params_ls$short_and_long_nm[1]
  }
  if(!missing("params_ls_ls")){
    warning("The params_ls_ls argument is now soft deprecated. Use input_params_ls instead.")
  }else{
    params_ls_ls <- input_params_ls
  }
  health_utl_long_nm_1L_chr <- ifelse(!is.null(params_ls_ls$short_and_long_nm),
                                      params_ls_ls$short_and_long_nm[2],
                                      NA_character_)
  study_descs_ls <- list(background_1L_chr = background_1L_chr,
                         coi_1L_chr = coi_1L_chr,
                         conclusion_1L_chr = conclusion_1L_chr,
                         ethics_1L_chr = ethics_1L_chr,
                         funding_1L_chr = funding_1L_chr,
                         health_utl_nm_1L_chr = health_utl_nm_1L_chr,
                         health_utl_long_nm_1L_chr = health_utl_long_nm_1L_chr,
                         time_btwn_bl_and_fup_1L_chr = time_btwn_bl_and_fup_1L_chr,
                         predr_ctgs_ls = predr_ctgs_ls,
                         sample_desc_1L_chr = sample_desc_1L_chr,
                         var_nm_change_lup = var_nm_change_lup)
  if(!is.null(input_params_ls)){
    input_params_ls$study_descs_ls <- study_descs_ls
  }else{
    input_params_ls <- study_descs_ls
  }
  return(input_params_ls)
}
make_ten_fold_text <- function(results_ls,
                               for_abstract_1L_lgl = F){
  mdls_chr <- get_ordered_sngl_csnl_mdls(results_ls)
  if(for_abstract_1L_lgl){
    text_1L_chr <- ifelse(length(mdls_chr)>1,
                          paste0(mdls_chr[1],
                                 " had the highest predictive ability in ten fold cross validation"),
                          "")
  }else{
    text_1L_chr <- ifelse(length(mdls_chr)>1,
                          paste0(mdls_chr[1],
                                 " had the highest predictive ability followed by ",
                                 get_ordered_sngl_csnl_mdls(results_ls,
                                                            select_int = -1,
                                                            collapse_1L_lgl = T),
                                 ". ",
                                 ifelse(length(mdls_chr)>2,
                                        paste0(mdls_chr[length(mdls_chr)],
                                               " had the least predictive capability."),
                                        "")),
                          paste0("the predictive ability of ", mdls_chr[1]))
  }
  return(text_1L_chr)
}
make_ten_folds_tbl_title <- function(results_ls,
                                     ref_1L_int = 1){
  title_1L_chr <- ifelse(ref_1L_int == 1,
                         paste0('10-fold cross-validated model fitting index for different ',
                                results_ls$ttu_cs_ls$best_mdl_types_ls %>%
                                  names() %>%
                                  paste0(collapse = ", ") %>%
                                  stringi::stri_replace_last(fixed = ",", " and"),
                                ' models using ',
                                results_ls$ttu_cs_ls$cs_mdls_predrs_seq_dscdng_chr[1],
                                ' as predictor with the baseline data'),
                         paste0('10-fold cross-validated model fitting index for different candidate predictors estimated using ',
                                results_ls$ttu_cs_ls$selected_mdls_chr[1],
                                ' with the baseline data'))
  return(title_1L_chr)
}
make_tfd_sngl_predr_mdls_tb <- function(outp_smry_ls,
                                        nbr_of_digits_1L_int = 2L,
                                        mdl_pfx_ls = list(OLS = "Ordinary Least Squares ",
                                                          GLM = c("Generalised Linear Mixed Model with ",
                                                                  "Beta Regression Model with Binomial "))){
  tfd_sngl_predr_mdls_tb <- mdl_pfx_ls %>%
    purrr::map2(names(mdl_pfx_ls),
                ~{
                  pfx_chr <- .x
                  mdls_tb <- outp_smry_ls$smry_of_sngl_predr_mdls_tb %>%
                    dplyr::filter(Model %>%
                                    purrr::map_lgl(~{
                                      term_1L_chr <- .x
                                      pfx_chr %>%
                                        purrr::map_lgl(~startsWith(term_1L_chr,.x)) %>%
                                        any()
                                    }))
                  pfx_chr %>%
                    purrr::map_dfr(~ {
                      pfx_1L_chr <- .x
                      mdls_tb %>%
                        dplyr::filter(startsWith(Model, pfx_1L_chr)) %>%
                        dplyr::mutate(Model = dplyr::case_when(Model %>% startsWith(mdl_pfx_ls[[2]][2]) ~ stringr::str_replace_all(Model,pfx_1L_chr,"Beta "),
                                                               T ~ stringr::str_remove_all(Model,pfx_1L_chr))
                        ) %>%
                        dplyr::mutate(Model = Model %>% stringr::str_remove_all("\\(|\\)"))
                    }) %>%
                    tibble::add_case(Model = .y, .before = 1)
                }) %>%
    purrr::map_dfr(~.x) %>%
    transform_tbl_to_rnd_vars(nbr_of_digits_1L_int = nbr_of_digits_1L_int)
  return(tfd_sngl_predr_mdls_tb)
}
make_tfmn_cmprsn_plt <-  function(data_tb,
                                  depnt_var_nm_1L_chr,
                                  dictionary_tb){
  tfmn_cmprsn_plt <- tidyr::gather(data_tb %>%
                                     dplyr::mutate(!!rlang::sym(paste0(depnt_var_nm_1L_chr,"_log")) := log(!!rlang::sym(depnt_var_nm_1L_chr)),
                                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_logit")) := psych::logit(!!rlang::sym(depnt_var_nm_1L_chr)),
                                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_loglog")) := -log(-log(!!rlang::sym(depnt_var_nm_1L_chr))),
                                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_cloglog")) := log(-log(1-!!rlang::sym(depnt_var_nm_1L_chr)))),
                                   variable,
                                   value,
                                   !!rlang::sym(depnt_var_nm_1L_chr),
                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_log")),
                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_logit")),
                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_loglog")),
                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_cloglog"))) %>%
    dplyr::mutate(variable = factor(variable, levels = paste0(depnt_var_nm_1L_chr,c("",
                                                                                    "_log", "_logit", "_loglog", "_cloglog")),
                                    labels= c("No transformation",
                                              "Log", "Logit","Log-log", "Complementary log-log"))) %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = variable)) +
    ggalt::geom_bkde() +
    ggplot2::geom_rug() +
    viridis::scale_fill_viridis(guide = "none", discrete = TRUE) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::theme_bw() + ggplot2::labs(x = paste0("Transformed ",
                                                   dictionary_tb %>%
                                                     ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "var_nm_chr",
                                                                                 match_value_xx = depnt_var_nm_1L_chr,
                                                                                 target_var_nm_1L_chr = "var_desc_chr",
                                                                                 evaluate_lgl = F)))
  return(tfmn_cmprsn_plt)
}
make_ttu_cs_ls <-  function(outp_smry_ls,
                            sig_covars_some_predrs_mdls_tb,
                            sig_thresh_covars_1L_chr){
  mdl_type_descs_chr <- outp_smry_ls$prefd_mdl_types_chr %>%
    purrr::map_chr(~ready4fun::get_from_lup_obj(outp_smry_ls$mdl_types_lup,
                                                match_var_nm_1L_chr = "short_name_chr",
                                                match_value_xx = .x,
                                                target_var_nm_1L_chr = "long_name_chr",
                                                evaluate_lgl = F))
  ttu_cs_ls <- list(best_mdl_types_ls = list(GLM = c("Gaussian distribution and log link"), # LEGACY ISSUE FROM MANUSCRIPT - OTHERWISE UNUSED
                                             OLS = c("no transformation","log transformation", "clog-log transformation")),
                    selected_mdls_chr = mdl_type_descs_chr %>%
                      purrr::map_chr(~paste0(ifelse(startsWith(.x,"Ordinary Least Squares"),"OLS","GLM"),
                                             " with ",
                                             ifelse(startsWith(.x,"Ordinary Least Squares"),
                                                    stringr::str_remove(.x,"Ordinary Least Squares ") %>%
                                                      stringr::str_sub(start = 2, end = -2) %>%
                                                      tolower(),
                                                    stringr::str_remove(.x,"Generalised Linear Mixed Model with ") %>%
                                                      stringr::str_remove("Beta Regression Model with ")))) ,
                    cs_mdls_predrs_seq_dscdng_chr = outp_smry_ls$smry_of_mdl_sngl_predrs_tb$Predictor,
                    sig_covars_all_predrs_mdls_chr = outp_smry_ls$signt_covars_chr,
                    sig_thresh_covars_1L_chr = sig_thresh_covars_1L_chr,
                    sig_covars_some_predrs_mdls_tb = sig_covars_some_predrs_mdls_tb,
                    rf_seq_dscdng_chr = outp_smry_ls$predr_cmprsn_tb$predr_chr,
                    mdl_predrs_and_rf_seqs_cmprsn_1L_chr = ifelse(outp_smry_ls$smry_of_mdl_sngl_predrs_tb$Predictor[1] == outp_smry_ls$predr_cmprsn_tb$predr_chr[1],"is consistent","contrasts") # COMPUTE
  )
  return(ttu_cs_ls)
}
make_uid_rename_lup <- function(data_tb,
                                id_var_nm_1L_chr = "UID"){
  uid_rename_lup_tb <- tibble::tibble(old_id_xx = data_tb %>%
                                        dplyr::pull(id_var_nm_1L_chr) %>%
                                        unique(),
                                      new_id_int = 1:length(data_tb %>%
                                                              dplyr::pull(id_var_nm_1L_chr) %>%
                                                              unique()))
  return(uid_rename_lup_tb)
}
make_unique_ls_elmt_idx_int <- function (data_ls)
{
  combos_tb <- tibble::as_tibble(data_ls, .name_repair = ~paste0("r_",
                                                                 1:length(data_ls))) %>% t() %>% as.data.frame()
  combos_tb <- combos_tb %>% tibble::as_tibble()
  if(ncol(combos_tb)>1){
    combos_tb <- combos_tb %>% dplyr::mutate(V2 = dplyr::case_when(V1 == V2 ~ NA_character_,
                                                                   T ~ V2)) %>%
      dplyr::mutate(combo_chr = purrr::map2_chr(V1,
                                                V2, ~ifelse(ncol(combos_tb) == 1 | is.na(.y),
                                                            .x,
                                                            paste0(.x,"_", .y))))
    combos_tb <- combos_tb %>% dplyr::group_by(combo_chr) %>%
      dplyr::mutate(combo_id = dplyr::row_number())
    unique_ls_elmt_idx_int <- purrr::map(data_ls %>% unique(),
                                         ~ready4fun::get_from_lup_obj(combos_tb %>% dplyr::ungroup(),
                                                                      match_var_nm_1L_chr = "combo_chr", match_value_xx = paste0(.x[1],
                                                                                                                                 ifelse(is.na(.x[2]), "", paste0("_", .x[2]))),
                                                                      target_var_nm_1L_chr = "combo_id", evaluate_lgl = F)) %>%
      purrr::flatten_int()
  }else{
    unique_ls_elmt_idx_int <- 1
  }
  return(unique_ls_elmt_idx_int)
}
make_valid_params_ls_ls <- function(analysis_core_params_ls,
                                    ds_tb,
                                    path_params_ls,
                                    maui_params_ls,
                                    candidate_covar_nms_chr = NA_character_,
                                    prefd_covars_chr = NULL,
                                    prefd_mdl_types_chr = NULL,
                                    raw_ds_tfmn_fn = NULL,
                                    scndry_analysis_extra_vars_chr = NA_character_,
                                    subtitle_1L_chr = "Methods Report 1: Analysis Program (Primary Analysis)",
                                    utl_class_fn_1L_chr = "as.numeric"){
  if(!missing("candidate_covar_nms_chr"))
    warning("The candidate_covar_nms_chr argument is soft deprecated. We recommend declaring the value for this argument when using the make_analysis_core_params_ls function.")
  if(!missing("prefd_covars_chr"))
    warning("The prefd_covars_chr argument is soft deprecated. We recommend declaring the value for this argument when using the make_analysis_core_params_ls function.")
  if(!missing("prefd_mdl_types_chr"))
    warning("The prefd_mdl_types_chr argument is soft deprecated. We recommend declaring the value for this argument when using the make_analysis_core_params_ls function.")
  valid_params_ls_ls <- make_prmry_analysis_params_ls(analysis_core_params_ls = analysis_core_params_ls,
                                                      candidate_covar_nms_chr = candidate_covar_nms_chr,
                                                      ds_tb = ds_tb,
                                                      path_params_ls = path_params_ls,
                                                      maui_params_ls = maui_params_ls,
                                                      prefd_covars_chr = prefd_covars_chr,
                                                      prefd_mdl_types_chr = prefd_mdl_types_chr,
                                                      raw_ds_tfmn_fn = raw_ds_tfmn_fn,
                                                      subtitle_1L_chr = subtitle_1L_chr,
                                                      utl_class_fn_1L_chr = utl_class_fn_1L_chr) %>%
    transform_params_ls_to_valid(scndry_analysis_extra_vars_chr = scndry_analysis_extra_vars_chr)
  valid_params_ls_ls$params_ls$short_and_long_nm <- NULL
  valid_params_ls_ls$short_and_long_nm <- maui_params_ls$short_and_long_nm
  valid_params_ls_ls$path_params_ls <- path_params_ls
  return(valid_params_ls_ls)
}
make_within_between_ratios_text <- function(results_ls){
  text_1L_chr <- results_ls$ttu_lngl_ls$cs_ts_ratios_tb %>%
    purrr::pmap_chr(~paste0(..2, " for ", ..1)) %>%
    paste0(collapse = ", ") %>%
    stringi::stri_replace_last(fixed = ",", " and")
  return(text_1L_chr)
}
