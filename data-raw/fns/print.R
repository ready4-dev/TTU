print_all_plts_for_mdl_set <- function (output_ls, start_from_1L_int = 0L)
{
    label_refs_mat <- paste0("lab", 1:(output_ls %>% purrr::flatten() %>%
        length()) + start_from_1L_int) %>% matrix(nrow = length(output_ls),
        byrow = T)
    purrr::walk(1:length(output_ls), ~print_ts_mdl_plts(output_ls[[.x]],
        title_1L_chr = names(output_ls)[.x], label_refs_chr = label_refs_mat[.x,
            ]))
}
print_cohort_table <- function(params_ls,
                               caption_1L_chr,
                               mkdn_tbl_ref_1L_chr){
  results_ls <- params_ls$results_ls
  df <- results_ls$tables_ls$participant_descs
  df$variable <- gsub("\\s*\\([^\\)]+\\)","",df$variable)
  df <- df %>%
    dplyr::mutate(variable = variable %>%
                    purrr::map_chr(~Hmisc::capitalize(.x)))
  if(params_ls$output_type_1L_chr == "PDF"){
    df <- df %>%
      dplyr::mutate_all(~ stringr::str_replace(.x,"%","\\\\%") %>%
                          stringr::str_replace(",","\\\\,"))
  }
  if(params_ls$output_type_1L_chr == "PDF"){
    names(df) <- c("","",
                   "(N =",paste0(results_ls$cohort_ls$n_inc_1L_dbl,")"),
                   "(N =",paste0(results_ls$cohort_ls$n_fup_1L_dbl,")"))
    df %>%
      kableExtra::kbl(booktabs = T,
                      caption = knitr::opts_current$get("tab.cap"),
                      escape = F) %>%
      kableExtra::kable_styling() %>%
      kableExtra::column_spec(3:6, width = "3em") %>%
      kableExtra::column_spec(1, bold = T, width = "14em") %>%
      kableExtra::add_header_above(c(" ", " ", "Baseline" = 2, "Follow-Up" = 2)) %>%
      kableExtra::collapse_rows(columns = 1)
  }else{
    df <- df %>% youthvars::transform_tb_for_merged_col_1(output_type_1L_chr = params_ls$output_type_1L_chr)
    add_to_row_ls <- make_bl_fup_add_to_row_ls(df,
                                               n_at_bl_1L_int = results_ls$cohort_ls$n_inc_1L_dbl,
                                               n_at_fup_1L_int = results_ls$cohort_ls$n_fup_1L_dbl)
    df %>%
      ready4show::print_table(output_type_1L_chr = params_ls$output_type_1L_chr,
                              caption_1L_chr = caption_1L_chr,
                              mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                              use_rdocx_1L_lgl = ifelse(params_ls$output_type_1L_chr=="Word",T,F),
                              add_to_row_ls = add_to_row_ls,
                              sanitize_fn = force)
  }
}
print_corls_tbl <- function(params_ls,
                            caption_1L_chr,
                            mkdn_tbl_ref_1L_chr){
  results_ls <- params_ls$results_ls
  tb <- results_ls$tables_ls$pred_dist_and_cors
  tb <- tb %>%
    dplyr::mutate(label = label %>%
                    purrr::map_chr(~stringr::str_remove_all(.x," \\(weighted total\\)")))
  if(params_ls$output_type_1L_chr == "PDF"){
    names(tb) <- c("","",
                   "(N =",
                   paste0(results_ls$cohort_ls$n_inc_1L_dbl,")"),
                   "(N =",paste0(results_ls$cohort_ls$n_fup_1L_dbl,")"),
                   "\\textit{p}")
    tb %>%
      kableExtra::kbl(booktabs = T,
                      caption = knitr::opts_current$get("tab.cap"),
                      escape = F) %>%
      kableExtra::kable_styling() %>%
      kableExtra::column_spec(3:6, width = "3em") %>%
      kableExtra::column_spec(1, bold = T, width = "14em") %>%
      kableExtra::add_header_above(c(" ", " ", "Baseline" = 2, "Follow-Up" = 2, " ")) %>%
      kableExtra::collapse_rows(columns = 1)
  }else{
    tb <- tb %>%
      youthvars::transform_tb_for_merged_col_1(output_type_1L_chr = params_ls$output_type_1L_chr)
    tb %>%
      ready4show::print_table(output_type_1L_chr = params_ls$output_type_1L_chr,
                              caption_1L_chr = caption_1L_chr,
                              mkdn_tbl_ref_1L_chr = mkdn_tbl_ref_1L_chr,
                              use_rdocx_1L_lgl = ifelse(params_ls$output_type_1L_chr=="Word",T,F),
                              add_to_row_ls = add_to_row_ls,
                              sanitize_fn = force)
  }
}
print_covar_ttu_tbls <- function(params_ls,
                                 caption_1L_chr,
                                 table_1L_chr,
                                 ref_1L_int = 1){
  results_ls <- params_ls$results_ls
  df <- results_ls$tables_ls %>%
    purrr::pluck(paste0("mdl_type_",ref_1L_int,"_covar_mdls_tb"))
  if(!is.null(df)){
    df <- df %>%
      dplyr::mutate(Parameter = Parameter %>% purrr::map_chr(~stringr::str_replace_all(.x,"_"," ")))
    df %>%
      print_lngl_ttu_tbls(caption_1L_chr = caption_1L_chr,
                          params_ls = params_ls,
                          ref_1L_int = ref_1L_int,
                          table_1L_chr = table_1L_chr)
  }
}
print_indpnt_predrs_coefs_tbl <- function(params_ls,
                                          caption_1L_chr,
                                          mkdn_tbl_ref_1L_chr){
  results_ls <- params_ls$results_ls
  tb <- results_ls$tables_ls$ind_preds_coefs_tbl
  add_to_row_ls <- list()
  add_to_row_ls$pos <- list(-1, 0, nrow(tb))
  add_to_row_ls$command <- c(paste0("\\toprule \n",
                                    "&\\multicolumn{5}{c}{",
                                    paste0(results_ls$ttu_lngl_ls$best_mdls_tb[[1,"model_type"]],
                                           " - ",
                                           results_ls$ttu_lngl_ls$best_mdls_tb[[1,"link_and_tfmn_chr"]]),
                                    "}&\\multicolumn{5}{c}{",
                                    paste0(results_ls$ttu_lngl_ls$best_mdls_tb[[2,"model_type"]],
                                           " - ",                                             results_ls$ttu_lngl_ls$best_mdls_tb[[2,"link_and_tfmn_chr"]]),
                                    "}\\\\\n"),
                             paste0("Parameter & Estimate	& SE	& 95CI & R2	& Sigma & Estimate & SE	& 95CI & R2 & Sigma \\\\\n",
                                    "\\midrule \n"),
                             paste0("\\hline\n","{\\footnotesize ",
                                    make_scaling_text(results_ls),
                                    "}\n")
  )
  if(params_ls$output_type_1L_chr =="Word"){
    tb$Parameter <- stringr::str_replace_all(stringr::str_replace_all(stringr::str_replace_all(tb$Parameter, '\\\\textbf', ''), '\\{', ''), '\\}', '')
  }
  if(params_ls$output_type_1L_chr == "PDF"){
    tb$Parameter <- tb$Parameter %>% purrr::map_chr(~ifelse(endsWith(.x,
                                                                     " model"),
                                                            paste0("\\textbf{",.x,"}"),.x))
  }
  tb %>%
    ready4show::print_table(output_type_1L_chr = params_ls$output_type_1L_chr,
                            caption_1L_chr = caption_1L_chr,
                            mkdn_tbl_ref_1L_chr = mkdn_tbl_ref_1L_chr,
                            use_rdocx_1L_lgl = ifelse(params_ls$output_type_1L_chr=="Word",T,F),
                            add_to_row_ls = add_to_row_ls,
                            footnotes_chr = make_scaling_text(results_ls),
                            sanitize_fn = force)

}
print_indpnt_predrs_lngl_mdl_coefs <- function(params_ls,
                                               caption_1L_chr,
                                               ref_1L_int = 1,
                                               table_1L_chr){
  results_ls <- params_ls$results_ls
  mdl_type_1L_chr <- results_ls$mdl_ingredients_ls$mdls_lup$mdl_type_chr %>%
    unique() %>%
    purrr::pluck(ref_1L_int)
  tb <- results_ls$tables_ls$ind_preds_coefs_tbl %>%
    dplyr::select(Parameter, names(.) %>%
                    purrr::keep(~{
                      endsWith(.x,mdl_type_1L_chr)})) %>%
    dplyr::rename_with(.cols = names(.) %>%
                         purrr::keep(~{
                           endsWith(.x,mdl_type_1L_chr)}),
                       ~ stringr::str_remove(.x,paste0("_",mdl_type_1L_chr)))
  tb %>%
    print_lngl_ttu_tbls(caption_1L_chr = caption_1L_chr,
                        params_ls = params_ls,
                        ref_1L_int = ref_1L_int,
                        table_1L_chr = table_1L_chr)
}
print_lngl_ttu_tbls <- function(table_df,
                                params_ls,
                                caption_1L_chr,
                                table_1L_chr,
                                ref_1L_int = 1){
  results_ls <- params_ls$results_ls
  if(params_ls$output_type_1L_chr == "PDF"){
    add_to_row_ls <- list()
    add_to_row_ls$pos <- list(0, nrow(table_df))
    add_to_row_ls$command <- c("Parameter & Estimate	& SE	& 95CI & R2	& Sigma\\\\\n",
                               paste0("\\hline\n",
                                      "{\\footnotesize ",
                                      make_scaling_text(results_ls,
                                                        table_1L_chr = knitr::opts_current$get("tab.id")),
                                      "}\n"))
    table_df$Parameter <- table_df$Parameter %>% purrr::map_chr(~ifelse(endsWith(.x,
                                                                                 " model"),
                                                                        paste0("\\textbf{",.x,"}"),.x))
  }else{
    add_to_row_ls <- NULL
  }
  table_df %>%
    ready4show::print_table(output_type_1L_chr = params_ls$output_type_1L_chr,
                            caption_1L_chr = caption_1L_chr,
                            mkdn_tbl_ref_1L_chr = paste0("tab:",table_1L_chr),
                            use_rdocx_1L_lgl = ifelse(params_ls$output_type_1L_chr=="Word",T,F),
                            add_to_row_ls = add_to_row_ls,
                            footnotes_chr = make_scaling_text(results_ls,
                                                              table_1L_chr = table_1L_chr),
                            hline_after_ls = c(-1,0),
                            sanitize_fn = force)

}
print_ten_folds_tbl <- function(params_ls,
                                caption_1L_chr,
                                mkdn_tbl_ref_1L_chr,
                                ref_1L_int = 1){
  results_ls <- params_ls$results_ls
  if(ref_1L_int ==1){
    df <- results_ls$tables_ls$tenf_sngl_predr_tb %>%
      dplyr::mutate(Model = gsub('"', '', Model)) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ .x %>%
                                    stringr::str_replace_all("  NA", NA_character_)))
  }else{
    df <- results_ls$tables_ls$tenf_prefd_mdl_tb
  }
  #df <- df
  if(params_ls$output_type_1L_chr == "PDF"){
    add_to_row_ls <- list()
    add_to_row_ls$pos <- list(0,0,0,nrow(df))
    add_to_row_ls$command <- c("&\\multicolumn{3}{c}{Training model fit}&\\multicolumn{3}{c}{Testing model fit}\\\\\n",
                               "&\\multicolumn{3}{c}{(averaged over 10 folds)}&\\multicolumn{3}{c}{(averaged over 10 folds)}\\\\\n",
                               "Model & R2	& RMSE & MAE &  R2	& RMSE & MAE  \\\\\n",
                               paste0("\\hline\n","{\\footnotesize  RMSE: Root Mean Squared Error; MAE: Mean Absolute Error}\n"))
    if(ref_1L_int ==1){
      df$Model <- df$Model %>% purrr::map_chr(~ifelse(.x %in% c("GLM","OLS"),
                                                      paste0("\\textbf{",
                                                             .x,
                                                             "}"),.x))
    }else{
      df$Predictor <- df$Predictor %>% purrr::map_chr(~paste0("\\textbf{",.x,"}"))
    }
  }else{
    add_to_row_ls <- NULL
  }
  df %>%
    ready4show::print_table(output_type_1L_chr = params_ls$output_type_1L_chr,
                            caption_1L_chr = caption_1L_chr,
                            mkdn_tbl_ref_1L_chr = mkdn_tbl_ref_1L_chr,
                            use_rdocx_1L_lgl = ifelse(params_ls$output_type_1L_chr=="Word",T,F),
                            add_to_row_ls = add_to_row_ls,
                            hline_after_ls = c(-1,0),
                            sanitize_fn = force)

}
print_ts_mdl_plts <- function (paths_to_plts_chr, title_1L_chr, label_refs_chr, mdl_smry_ls)
{
    cat("\n")
    cat("## ", title_1L_chr, "\n")
    cat("\n")
    purrr::pwalk(list(paths_to_plts_chr, paste0("Caption", 1:length(paths_to_plts_chr)),
        label_refs_chr), ~{
        cat(paste0("\n![", "](", ..1, ")\n\n"))
        cat(fig_nums(..3, ..2))
        cat("\n\n")
    })
    cat("\n\n")
}
