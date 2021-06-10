add_prefd_predr_var_to_mdl_smry_ls <- function(mdl_smry_ls,
                                               ds_smry_ls){
  mdl_smry_ls$predr_var_nm_1L_chr <- ds_smry_ls$candidate_predrs_chr[1]
  mdl_smry_ls$predr_var_desc_1L_chr <- ds_smry_ls$predictors_lup %>% ready4fun::get_from_lup_obj(match_value_xx = mdl_smry_ls$predr_var_nm_1L_chr, match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = "long_name_chr",
                                                                                                 evaluate_lgl = F)
  mdl_smry_ls$predr_vals_dbl <- make_predr_vals(mdl_smry_ls$predr_var_nm_1L_chr,
                                                candidate_predrs_lup = ds_smry_ls$predictors_lup)
  return(mdl_smry_ls)
}
add_tfmd_var_to_ds <- function(data_tb,
                               depnt_var_nm_1L_chr,
                               tfmn_1L_chr){
  data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_depnt_var_nm(depnt_var_nm_1L_chr,
                                                                                tfmn_1L_chr = tfmn_1L_chr)), !!rlang::sym(depnt_var_nm_1L_chr) %>%
                                              calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr)))
  return(data_tb)
}
add_uids_to_tbs_ls <- function (tbs_ls, prefix_1L_chr, id_var_nm_1L_chr = "fkClientID")
{
  participant_ids <- paste0(prefix_1L_chr, 1:nrow(tbs_ls[[1]])) %>%
    sample(nrow(tbs_ls[[1]]))
  tbs_ls <- purrr::map(tbs_ls, ~{
    .x %>% dplyr::mutate(`:=`(!!rlang::sym(id_var_nm_1L_chr),
                              tidyselect::all_of(participant_ids[1:nrow(.x)]))) %>%
      dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr) %>%
                       purrr::map_chr(~stringr::str_replace(.x, prefix_1L_chr,
                                                            "")) %>% as.numeric())
  }) %>% stats::setNames(names(tbs_ls))
  return(tbs_ls)
}
add_utility_predn_to_ds <- function (data_tb, model_mdl, tfmn_1L_chr, depnt_var_nm_1L_chr, predr_vars_nms_chr = NULL,
                                     force_min_max_1L_lgl = T, utl_min_val_1L_dbl = 0.03, impute_1L_lgl = T, utl_cls_fn = NULL,
                                     rmv_tfd_depnt_var_1L_lgl = F)
{
    dep_vars_chr <- c(depnt_var_nm_1L_chr, transform_depnt_var_nm(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        tfmn_1L_chr = tfmn_1L_chr)) %>% unique()
    data_tb <- purrr::reduce(dep_vars_chr, .init = data_tb, ~dplyr::mutate(.x,
                                                                           !!rlang::sym(.y):= NA_real_))
    predictions_dbl <- predict_utility(data_tb = data_tb, tfmn_1L_chr = tfmn_1L_chr,
                                       model_mdl = model_mdl,
                                       force_min_max_1L_lgl = force_min_max_1L_lgl,
                                       utl_min_val_1L_dbl = utl_min_val_1L_dbl,
                                       impute_1L_lgl = impute_1L_lgl,
                                       utl_cls_fn = utl_cls_fn)
    data_tb <- data_tb %>% dplyr::mutate(!!rlang::sym(depnt_var_nm_1L_chr):=predictions_dbl)
    if(!is.null(predr_vars_nms_chr)){
      data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(purrr::map(predr_vars_nms_chr, ~ paste0(.x,c("_baseline","_change"))) %>%
                                             purrr::flatten_chr()))
    }
    if(rmv_tfd_depnt_var_1L_lgl){
      data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(dep_vars_chr[dep_vars_chr!=depnt_var_nm_1L_chr]))
    }
    return(data_tb)
}

