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
add_utl_predn_to_new_ds <- function(data_tb,
                                    ingredients_ls,
                                    mdl_nm_1L_chr,
                                    deterministic_1L_lgl = T,
                                    force_min_max_1L_lgl = T,
                                    model_mdl = NULL,
                                    #impute_1L_lgl = T,
                                    new_data_is_1L_chr = "Simulated",
                                    predr_vars_nms_chr = NULL,
                                    round_var_nm_1L_chr = "Timepoint",
                                    round_bl_val_1L_chr = "BL",
                                    utl_cls_fn = NULL,
                                    utl_var_nm_1L_chr = NULL){
  if(is.null(model_mdl))
    model_mdl <- get_table_predn_mdl(mdl_nm_1L_chr,
                                     ingredients_ls = ingredients_ls)
  tfmn_1L_chr <- ready4fun::get_from_lup_obj(ingredients_ls$mdl_types_lup,
                                             match_value_xx = mdl_type_1L_chr,
                                             match_var_nm_1L_chr = "short_name_chr",
                                             target_var_nm_1L_chr = "tfmn_chr",
                                             evaluate_lgl = F)
  predn_type_1L_chr <- ready4fun::get_from_lup_obj(ingredients_ls$mdl_types_lup,
                                                   match_value_xx = mdl_type_1L_chr,
                                                   match_var_nm_1L_chr = "short_name_chr",
                                                   target_var_nm_1L_chr = "predn_type_chr",
                                                   evaluate_lgl = F)
  if(is.na(predn_type_1L_chr))
    predn_type_1L_chr <- NULL
  if(!is.null(predr_vars_nms_chr)){
    data_tb <- rename_from_nmd_vec(data_tb,
                                   nmd_vec_chr = predr_vars_nms_chr,
                                   vec_nms_as_new_1L_lgl = T)
  }
  mdl_predr_terms_chr <- ingredients_ls$mdls_lup %>%
    dplyr::filter(mdl_nms_chr == mdl_nm_1L_chr) %>%
    dplyr::pull(predrs_ls) %>%
    purrr::flatten_chr()
  original_ds_vars_chr <- names(data_tb)[!names(data_tb) %in% c(mdl_predr_terms_chr,
                                                                ifelse(!is.null(utl_var_nm_1L_chr),
                                                                       utl_var_nm_1L_chr,
                                                                       ingredients_ls$depnt_var_nm_1L_chr))]
  updated_tb <- data_tb %>%
    transform_ds_to_predn_ds(predr_vars_nms_chr = mdl_predr_terms_chr,
                             tfmn_1L_chr = tfmn_1L_chr,
                             depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr,
                             id_var_nm_1L_chr = ingredients_ls$id_var_nm_1L_chr,
                             round_var_nm_1L_chr = round_var_nm_1L_chr,
                             round_bl_val_1L_chr = round_bl_val_1L_chr,
                             predictors_lup = ingredients_ls$predictors_lup) %>%
    add_utility_predn_to_ds(model_mdl = model_mdl,
                            tfmn_1L_chr = tfmn_1L_chr,
                            depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr,
                            predr_vars_nms_chr = mdl_predr_terms_chr,
                            force_min_max_1L_lgl = force_min_max_1L_lgl,
                            force_new_data_1L_lgl = T,
                            impute_1L_lgl = T, # Redundant?
                            is_brms_mdl_1L_lgl = inherits(model_mdl,"brmsfit"),
                            new_data_is_1L_chr = new_data_is_1L_chr,
                            predn_type_1L_chr = NULL,
                            rmv_tfd_depnt_var_1L_lgl = T,
                            utl_cls_fn = utl_cls_fn,
                            utl_min_val_1L_dbl = ingredients_ls$utl_min_val_1L_dbl,
                            sd_dbl = get_random_intercept(ingredients_ls$mdls_smry_tb,
                                                          mdl_nm_1L_chr = mdl_nm_1L_chr,
                                                          deterministic_1L_lgl = deterministic_1L_lgl))
  if(!is.null(utl_var_nm_1L_chr)){
    updated_tb <- updated_tb %>%
      dplyr::rename(!!rlang::sym(utl_var_nm_1L_chr):=tidyselect::all_of(mdl_dep_var_1L_chr))
  }
  if(!is.null(names(predr_vars_nms_chr))){
    updated_tb <- rename_from_nmd_vec(updated_tb,
                                      nmd_vec_chr = predr_vars_nms_chr,
                                      vec_nms_as_new_1L_lgl = F)
  }
  names_to_incl_chr <- c(names(updated_tb),
                         setdiff(names(data_tb),
                                 names(updated_tb)))
  updated_tb <- dplyr::left_join(data_tb %>% dplyr::select(tidyselect::all_of(original_ds_vars_chr)),
                                 updated_tb)
  updated_tb <- updated_tb %>%
    dplyr::select(tidyselect::all_of(names_to_incl_chr[names_to_incl_chr %in% names(updated_tb)]))
  return(updated_tb)
}
add_utility_predn_to_ds <- function (data_tb,
                                     model_mdl,
                                     tfmn_1L_chr,
                                     depnt_var_nm_1L_chr,
                                     force_min_max_1L_lgl = T,
                                     force_new_data_1L_lgl = F,
                                     impute_1L_lgl = T,
                                     is_brms_mdl_1L_lgl = T,
                                     new_data_is_1L_chr = "Predicted",
                                     predn_type_1L_chr = NULL,
                                     predr_vars_nms_chr = NULL,
                                     rmv_tfd_depnt_var_1L_lgl = F,
                                     sd_dbl = NA_real_,
                                     utl_cls_fn = NULL,
                                     utl_min_val_1L_dbl = -1)
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
                                       new_data_is_1L_chr = new_data_is_1L_chr,
                                       utl_cls_fn = utl_cls_fn,
                                       is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl,
                                       force_new_data_1L_lgl = force_new_data_1L_lgl,
                                       predn_type_1L_chr = predn_type_1L_chr,
                                       sd_dbl = sd_dbl)
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

