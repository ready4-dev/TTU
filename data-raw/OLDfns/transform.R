transform_chr_digit_pairs <- function(digit_pairs_chr,
                                      nbr_of_digits_1L_int = 2L){
  tfd_digit_pairs_chr <- digit_pairs_chr %>%
    purrr::map_chr(~{
      # abs_vals_elmnts_chr <- .x  %>% regmatches(gregexpr("[[:digit:]]+", .)) %>% purrr::pluck(1)
      # abs_vals_chr <- c(paste0(abs_vals_elmnts_chr[1:2], collapse = "."),
      #                   paste0(abs_vals_elmnts_chr[3:4], collapse = "."))
      abs_vals_chr <- .x %>% strsplit(",") %>% purrr::pluck(1) %>% stringr::str_squish()
      abs_vals_chr[1] <- ifelse(startsWith(.x,paste0("-",abs_vals_chr[1])),paste0("-",abs_vals_chr[1]),abs_vals_chr[1])
      abs_vals_chr[2] <- ifelse(endsWith(.x,paste0("-",abs_vals_chr[2])),paste0("-",abs_vals_chr[2]),abs_vals_chr[2])
      as.numeric(abs_vals_chr) %>% round(digits=nbr_of_digits_1L_int) %>% format(nsmall=nbr_of_digits_1L_int) %>% paste0(collapse = ", ")
    })
  return(tfd_digit_pairs_chr)
}
transform_data_tb_for_cmprsn <- function (data_tb, model_mdl, depnt_var_nm_1L_chr = "utl_total_w",
    source_data_nm_1L_chr = "Original", new_data_is_1L_chr = "Predicted",
    predn_type_1L_chr = NULL, family_1L_chr = NA_character_, impute_1L_lgl = F, is_brms_mdl_1L_lgl = F,
    sd_dbl = NA_real_, sfx_1L_chr = "", tfmn_for_bnml_1L_lgl = F, tfmn_1L_chr = "NTF", utl_cls_fn = NULL, utl_min_val_1L_dbl = NA_real_)
{
  new_data_dbl <- predict_utility(data_tb = data_tb,
                                  tfmn_1L_chr = tfmn_1L_chr,
                                  model_mdl = model_mdl,
                                  force_min_max_1L_lgl = !is.na(utl_min_val_1L_dbl),
                                  force_new_data_1L_lgl = T,
                                  utl_min_val_1L_dbl = utl_min_val_1L_dbl,
                                  impute_1L_lgl = impute_1L_lgl,
                                  utl_cls_fn = utl_cls_fn,
                                  new_data_is_1L_chr = new_data_is_1L_chr,
                                  predn_type_1L_chr = predn_type_1L_chr,
                                  sd_dbl = sd_dbl,
                                  tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl,
                                  family_1L_chr = family_1L_chr,
                                  is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl)
    tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_predd_var_nm(new_data_is_1L_chr,
                                                                                      sfx_1L_chr = sfx_1L_chr,
                                                                                      utl_min_val_1L_dbl = utl_min_val_1L_dbl)),
        new_data_dbl), `:=`(!!rlang::sym(source_data_nm_1L_chr),
        !!rlang::sym(depnt_var_nm_1L_chr)))
    return(tfd_data_tb)
}

transform_depnt_var_nm <- function (depnt_var_nm_1L_chr, tfmn_1L_chr = "NTF")
{
    tfd_depnt_var_nm_1L_chr <- paste0(depnt_var_nm_1L_chr,
                                      ifelse(tfmn_1L_chr == "NTF",
                                             "",
                                             paste0("_", tfmn_1L_chr)))
    return(tfd_depnt_var_nm_1L_chr)
}
# transform_depnt_var_nm_for_cll <- function (depnt_var_nm_1L_chr)
# {
#     tfd_depnt_var_nm_1L_chr <- paste0(depnt_var_nm_1L_chr, "_cloglog")
#     return(tfd_depnt_var_nm_1L_chr)
# }
transform_dict_with_rename_lup <- function(dictionary_tb,
                                           rename_lup){
  var_lbl_1L_chr <- Hmisc::label(dictionary_tb$var_nm_chr)
  tfd_dictionary_tb <- dictionary_tb %>%
    dplyr::mutate(var_nm_chr = var_nm_chr %>%
                    purrr::map_chr(~ifelse(.x %in% rename_lup$old_nms_chr,
                                           ready4::get_from_lup_obj(rename_lup,
                                                                       match_value_xx = .x,
                                                                       match_var_nm_1L_chr = "old_nms_chr",
                                                                       target_var_nm_1L_chr = "new_nms_chr",
                                                                       evaluate_1L_lgl = F),
                                           .x)))
  Hmisc::label(tfd_dictionary_tb[["var_nm_chr"]]) <- var_lbl_1L_chr
  return(tfd_dictionary_tb)
}
transform_ds_for_all_cmprsn_plts <- function(tfd_data_tb,
                                             model_mdl,
                                             depnt_var_nm_1L_chr,
                                             is_brms_mdl_1L_lgl,
                                             predn_type_1L_chr,
                                             sd_dbl,
                                             sfx_1L_chr = "",
                                             tfmn_1L_chr,
                                             utl_min_val_1L_dbl = -1){
  tfd_data_tb <- transform_data_tb_for_cmprsn(tfd_data_tb %>% dplyr::ungroup(),
                                              model_mdl = model_mdl,
                                              depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                              predn_type_1L_chr = predn_type_1L_chr,
                                              sfx_1L_chr = sfx_1L_chr,
                                              tfmn_1L_chr = tfmn_1L_chr) %>%
    transform_data_tb_for_cmprsn(model_mdl = model_mdl,
                                 depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                 family_1L_chr = NA_character_,
                                 is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl,
                                 new_data_is_1L_chr = "Simulated",
                                 predn_type_1L_chr = predn_type_1L_chr,
                                 sd_dbl = sd_dbl,
                                 sfx_1L_chr = sfx_1L_chr,
                                 tfmn_1L_chr = tfmn_1L_chr,
                                 tfmn_for_bnml_1L_lgl = FALSE) %>%
    transform_data_tb_for_cmprsn(model_mdl = model_mdl,
                                 depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                 predn_type_1L_chr = predn_type_1L_chr,
                                 sfx_1L_chr = sfx_1L_chr,
                                 tfmn_1L_chr = tfmn_1L_chr,
                                 utl_min_val_1L_dbl = utl_min_val_1L_dbl) %>%
    transform_data_tb_for_cmprsn(model_mdl = model_mdl,
                                 depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                 family_1L_chr = NA_character_,
                                 is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl,
                                 new_data_is_1L_chr = "Simulated",
                                 predn_type_1L_chr = predn_type_1L_chr,
                                 sfx_1L_chr = sfx_1L_chr,
                                 sd_dbl = sd_dbl,
                                 tfmn_1L_chr = tfmn_1L_chr,
                                 tfmn_for_bnml_1L_lgl = FALSE,
                                 utl_min_val_1L_dbl = utl_min_val_1L_dbl)
  return(tfd_data_tb)
}
transform_ds_for_mdlng <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_var_nm_1L_chr,
    covar_var_nms_chr = NA_character_)
{
    mdl_vars_chr <- c(names(data_tb)[names(data_tb) %>% startsWith(depnt_var_nm_1L_chr)],
        predr_var_nm_1L_chr, covar_var_nms_chr) %>% purrr::discard(is.na)
    tfd_data_tb <- data_tb %>% tidyr::drop_na(!!!rlang::syms(mdl_vars_chr)) %>%
        dplyr::select(!!!rlang::syms(mdl_vars_chr))
    return(tfd_data_tb)
}
transform_ds_to_predn_ds <- function(data_tb,
                                     predr_vars_nms_chr,
                                     tfmn_1L_chr,
                                     depnt_var_nm_1L_chr,
                                     id_var_nm_1L_chr,
                                     round_var_nm_1L_chr,
                                     round_bl_val_1L_chr,
                                     predictors_lup){#=NULL
  # if(is.null(predictors_lup))
  #   utils::data("predictors_lup", package = "youthvars", envir = environment())
  data_tb <- data_tb %>%
    dplyr::mutate(!!rlang::sym(depnt_var_nm_1L_chr):= NA_real_)
  data_tb <- purrr::reduce(predr_vars_nms_chr,
                           .init = data_tb,
                           ~ {
                             predr_cls_fn <- eval(parse(text=ready4::get_from_lup_obj(predictors_lup,
                                                                                         match_var_nm_1L_chr = "short_name_chr",
                                                                                         match_value_xx = .y,
                                                                                         target_var_nm_1L_chr = "class_fn_chr",
                                                                                         evaluate_1L_lgl = F)))
                             dplyr::mutate(.x,
                                           !!rlang::sym(.y) := !!rlang::sym(.y) %>% rlang::exec(.fn = predr_cls_fn))
                           })
  data_tb <- data_tb %>% transform_tb_to_mdl_inp(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                      predr_vars_nms_chr = predr_vars_nms_chr,
                                                      id_var_nm_1L_chr = id_var_nm_1L_chr,
                                                      round_var_nm_1L_chr = round_var_nm_1L_chr,
                                                      round_bl_val_1L_chr = round_bl_val_1L_chr,
                                                      drop_all_msng_1L_lgl = F,
                                                      scaling_fctr_dbl = purrr::map_dbl(predr_vars_nms_chr,
                                                                                        ~ ready4::get_from_lup_obj(predictors_lup,
                                                                                                                      target_var_nm_1L_chr = "mdl_scaling_dbl",
                                                                                                                      match_var_nm_1L_chr = "short_name_chr",
                                                                                                                      match_value_xx = .x,
                                                                                                                      evaluate_1L_lgl = F)),
                                                      ungroup_1L_lgl = T,
                                                      tfmn_1L_chr = tfmn_1L_chr)
  return(data_tb)
}
# transform_ds_with_rename_lup <- function(ds_tb,
#                                          rename_lup,
#                                          target_var_nms_chr = NULL){
#   if(is.null(target_var_nms_chr))
#     target_var_nms_chr <- intersect(names(ds_tb),rename_lup$old_nms_chr)
#   tfd_ds_tb <- dplyr::rename_with(ds_tb,
#                                    .cols = target_var_nms_chr,
#                                    ~ ready4::get_from_lup_obj(rename_lup,
#                                                                  match_value_xx = .x,
#                                                                  match_var_nm_1L_chr = "old_nms_chr",
#                                                                  target_var_nm_1L_chr = "new_nms_chr",
#                                                                  evaluate_1L_lgl = F))
#   return(tfd_ds_tb)
#
# }
transform_mdl_vars_with_clss <- function(ds_tb,
                                         predictors_lup = NULL,
                                         prototype_lup = NULL,
                                         depnt_var_nm_1L_chr = "utl_total_w", # Remove default
                                         class_fn_1L_chr = "as.numeric"){
  if(is.null(predictors_lup))
    data("predictors_lup", package = "youthvars", envir = environment())
  if(is.null(prototype_lup))
    prototype_lup <- ready4use::Ready4useRepos(gh_repo_1L_chr = "ready4-dev/ready4",
                                               gh_tag_1L_chr = "Documentation_0.0") %>%
      ready4::ingest(fls_to_ingest_chr = "prototype_lup",
                     metadata_1L_lgl = F)
  if(!is.null(depnt_var_nm_1L_chr)){
    predictors_lup <- tibble::add_case(predictors_lup,
                                       short_name_chr = depnt_var_nm_1L_chr,
                                       class_chr = "numeric",
                                       class_fn_chr = class_fn_1L_chr)
  }
  tfd_ds_tb <- purrr::reduce(predictors_lup$short_name_chr,
                             .init = ds_tb,
                             ~ if(.y %in% names(.x)){
                               label_1L_chr <- Hmisc::label(.x[[.y]])
                               class_1L_chr <- ready4::get_from_lup_obj(predictors_lup,
                                                                        match_var_nm_1L_chr = "short_name_chr",
                                                                        match_value_xx = .y,
                                                                        target_var_nm_1L_chr = "class_chr",
                                                                        evaluate_1L_lgl = F)
                               ns_1L_chr <- ready4::get_from_lup_obj(prototype_lup,
                                                                     match_var_nm_1L_chr = "type_chr",
                                                                     match_value_xx = class_1L_chr,
                                                                     target_var_nm_1L_chr = "pt_ns_chr",
                                                                     evaluate_1L_lgl = F)
                               ns_and_ext_1L_chr <- ifelse(ns_1L_chr == "base",
                                                           "",
                                                           paste0(ns_1L_chr,"::"))
                               fn <- ifelse(exists(paste0("as.",class_1L_chr), where = paste0("package:",ns_1L_chr)),
                                            eval(parse(text = paste0(ns_and_ext_1L_chr,"as.",class_1L_chr))),
                                            ifelse(exists(paste0("as_",class_1L_chr), where = paste0("package:",ns_1L_chr)),
                                                   eval(parse(text = paste0(ns_and_ext_1L_chr,"as_",class_1L_chr))),
                                                   eval(parse(text = paste0(ns_and_ext_1L_chr,class_1L_chr)))))
                               tb <- .x %>% dplyr::mutate(!!rlang::sym(.y) := rlang::exec(ready4::get_from_lup_obj(predictors_lup,
                                                                                                                   match_var_nm_1L_chr = "short_name_chr",
                                                                                                                   match_value_xx = .y,
                                                                                                                   target_var_nm_1L_chr = "class_fn_chr",
                                                                                                                   evaluate_1L_lgl = T),
                                                                                          !!rlang::sym(.y) %>% fn))
                               if(label_1L_chr != ""){
                                 Hmisc::label(tb[[.y]]) <- label_1L_chr
                               }
                               tb
                             }else{
                               .x
                             })
  return(tfd_ds_tb)
}
transform_names <- function(names_chr,
                            rename_lup,
                            invert_1L_lgl = F){
  new_names_chr <- names_chr %>%
    purrr::map_chr(~ifelse((!invert_1L_lgl & .x %in% rename_lup$old_nms_chr) | (invert_1L_lgl & .x %in% rename_lup$new_nms_chr),
                           .x %>%
                             ready4::get_from_lup_obj(data_lookup_tb = rename_lup,
                                                         match_var_nm_1L_chr = ifelse(invert_1L_lgl,"new_nms_chr","old_nms_chr"),
                                                         target_var_nm_1L_chr = ifelse(invert_1L_lgl,"old_nms_chr","new_nms_chr"),
                                                         evaluate_1L_lgl = F),
                           .x))
  return(new_names_chr)
}
transform_nms_in_mdl_tbl <- function(mdl_tbl_tb,
                                     col_nm_1L_chr = "Parameter",
                                     var_nm_change_lup = NULL){
  if(is.null(var_nm_change_lup)){
    tfd_mdl_tbl_tb <-  mdl_tbl_tb
  }else{
    tfd_mdl_tbl_tb <-  mdl_tbl_tb %>%
      dplyr::mutate(!!rlang::sym(col_nm_1L_chr) := dplyr::case_when(!!rlang::sym(col_nm_1L_chr) %>%
                                                                      purrr::map_lgl(~(endsWith(.x," model") | endsWith(.x," baseline") | endsWith(.x," change"))) ~ !!rlang::sym(col_nm_1L_chr) %>% purrr::map_chr(~{
                                                                        sfx_starts_1L_int <- stringi::stri_locate_first_fixed(.x," ")[[1,1]]
                                                                        paste0(stringr::str_sub(.x,end=(sfx_starts_1L_int-1)) %>%
                                                                                 strsplit("_") %>%
                                                                                 purrr::pluck(1) %>%
                                                                                 transform_names(rename_lup = var_nm_change_lup) %>%
                                                                                 paste0(collapse = "_"),
                                                                               stringr::str_sub(.x,start=sfx_starts_1L_int))}),
                                                                    T ~ !!rlang::sym(col_nm_1L_chr)))
  }
 return(tfd_mdl_tbl_tb)
}
transform_params_ls_to_valid <- function(params_ls,
                                         scndry_analysis_extra_vars_chr = NA_character_){
  target_var_nms_chr <- c(params_ls$ds_descvs_ls$candidate_predrs_chr,
                          params_ls$candidate_covar_nms_chr,
                          scndry_analysis_extra_vars_chr) %>%
    purrr::discard(is.na) %>%
    unique()
  valid_var_nms_chr <- target_var_nms_chr %>%
    stringi::stri_replace_last_fixed("_dbl","") %>%
    stringi::stri_replace_last_fixed("_int","") %>%
    stringi::stri_replace_all_fixed("_","")
  unchanged_var_nms_chr <- setdiff(params_ls$ds_descvs_ls$dictionary_tb$var_nm_chr,
                                   target_var_nms_chr)
  rename_lup <- tibble::tibble(old_nms_chr = c(unchanged_var_nms_chr,target_var_nms_chr),
                               new_nms_chr = make.unique(c(unchanged_var_nms_chr,
                                                           valid_var_nms_chr), sep="V")) %>%
    dplyr::filter(!old_nms_chr %in% unchanged_var_nms_chr)
  params_ls$ds_tb <- youthvars::transform_ds_with_rename_lup(params_ls$ds_tb,
                                                  rename_lup = rename_lup,
                                                  target_var_nms_chr = target_var_nms_chr)
    # dplyr::rename_with(params_ls$ds_tb,
    #                                     .cols = target_var_nms_chr,
    #                                     ~ ready4::get_from_lup_obj(rename_lup,
    #                                                                   match_value_xx = .x,
    #                                                                   match_var_nm_1L_chr = "old_nms_chr",
    #                                                                   target_var_nm_1L_chr = "new_nms_chr",
    #                                                                   evaluate_1L_lgl = F))
  params_ls$ds_descvs_ls$dictionary_tb <- params_ls$ds_descvs_ls$dictionary_tb %>%
    transform_dict_with_rename_lup(rename_lup = rename_lup)
  # var_lbl_1L_chr <- Hmisc::label(params_ls$ds_descvs_ls$dictionary_tb$var_nm_chr)
  # params_ls$ds_descvs_ls$dictionary_tb <- params_ls$ds_descvs_ls$dictionary_tb %>%
  #   dplyr::mutate(var_nm_chr = var_nm_chr %>%
  #                   purrr::map_chr(~ifelse(.x %in% rename_lup$old_nms_chr,
  #                                          ready4::get_from_lup_obj(rename_lup,
  #                                                                      match_value_xx = .x,
  #                                                                      match_var_nm_1L_chr = "old_nms_chr",
  #                                                                      target_var_nm_1L_chr = "new_nms_chr",
  #                                                                      evaluate_1L_lgl = F),
  #                                          .x)))
  # Hmisc::label(params_ls$ds_descvs_ls$dictionary_tb[["var_nm_chr"]]) <- var_lbl_1L_chr#
  rename_lup <- rename_lup %>%
    dplyr::filter(old_nms_chr != new_nms_chr)
  valid_params_ls_ls <- list(params_ls = params_ls %>%
                               transform_params_ls_from_lup(rename_lup = rename_lup),
                             rename_lup = rename_lup)
  return(valid_params_ls_ls)
}
transform_params_ls_from_lup <- function(params_ls,
                                         rename_lup){
  if(!is.null(params_ls$ds_descvs_ls)){
    params_ls$ds_descvs_ls$candidate_predrs_chr <- params_ls$ds_descvs_ls$candidate_predrs_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_1L_lgl = F)))
    params_ls$ds_descvs_ls$cohort_descv_var_nms_chr <- params_ls$ds_descvs_ls$cohort_descv_var_nms_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_1L_lgl = F)))
  }
  if(!is.null(params_ls$predictors_lup)){
    params_ls$predictors_lup$short_name_chr <-  params_ls$predictors_lup$short_name_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_1L_lgl = F)))
  }
  params_ls$candidate_covar_nms_chr <- params_ls$candidate_covar_nms_chr %>%
    purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                           .x,
                           ready4::get_from_lup_obj(rename_lup,
                                                match_value_xx = .x,
                                                match_var_nm_1L_chr = "old_nms_chr",
                                                target_var_nm_1L_chr = "new_nms_chr",
                                                evaluate_1L_lgl = F)))
  if(!is.na(params_ls$prefd_covars_chr)){
    params_ls$prefd_covars_chr <- params_ls$prefd_covars_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_1L_lgl = F)))
  }
  if(!is.null(params_ls$candidate_predrs_chr)){
    params_ls$candidate_predrs_chr <- params_ls$candidate_predrs_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_1L_lgl = F)))
  }
  return(params_ls)
}
transform_paths_ls_for_scndry <- function(paths_ls,
                                          reference_1L_int = 1,
                                          remove_prmry_1L_lgl = F,
                                          remove_mkdn_1L_lgl = F){
  paths_ls$prmry_analysis_dir_nm_1L_chr <- paths_ls$write_to_dir_nm_1L_chr
  paths_ls$write_to_dir_nm_1L_chr <- paste0(paths_ls$write_to_dir_nm_1L_chr,
                                            "/secondary_",
                                            reference_1L_int)
  paths_ls$reports_dir_1L_chr <- paste0(paths_ls$reports_dir_1L_chr %>%
    stringr::str_sub(end = -(nchar(paths_ls$prmry_analysis_dir_nm_1L_chr)+10)),
    "/",
    paths_ls$write_to_dir_nm_1L_chr,
    "/Reports")
  if(remove_prmry_1L_lgl)
    paths_ls <- paths_ls[names(paths_ls) != "prmry_analysis_dir_nm_1L_chr"]
  if(remove_mkdn_1L_lgl)
    paths_ls <- paths_ls[names(paths_ls) != "reports_dir_1L_chr"]
  return(paths_ls)
}
transform_predd_var_nm <- function (new_data_is_1L_chr,
                                    sfx_1L_chr = "",
                                    utl_min_val_1L_dbl = NA_real_)
{
  tfd_predd_var_nm_1L_chr <- paste0(new_data_is_1L_chr,
                                     sfx_1L_chr,
                                     ifelse(!is.na(utl_min_val_1L_dbl),
                                            " (constrained)",
                                            ""))
  return(tfd_predd_var_nm_1L_chr)
}
transform_predr_nm_part_of_phrases <- function(phrases_chr,
                                               old_nms_chr = NULL,
                                               new_nms_chr = NULL){
  if(is.null(old_nms_chr)){
    tfd_phrases_chr <- phrases_chr
  }else{
    nm_changes_lup_tb = tibble::tibble(old_nms_chr = old_nms_chr,
                                       new_nms_chr = new_nms_chr)
    tfd_phrases_chr <- phrases_chr %>%
      purrr::map_chr(~{
        phrase_1L_chr <- .x
        match_lgl <- nm_changes_lup_tb$old_nms_chr %>%
          purrr::map_lgl(~stringr::str_detect(phrase_1L_chr, .x))
        if(any(match_lgl)){
          stringr::str_replace(phrase_1L_chr,
                               nm_changes_lup_tb$old_nms_chr[match_lgl],
                               nm_changes_lup_tb$new_nms_chr[match_lgl])
        }else{
          phrase_1L_chr
        }
      })
  }
  return(tfd_phrases_chr)
}
transform_rprt_lup <- function(rprt_lup,
                               add_suplry_rprt_1L_lgl = T,
                               add_sharing_rprt_1L_lgl = F,
                               start_at_int = NULL,
                               reference_1L_int = NULL){
  if(add_suplry_rprt_1L_lgl){
    rprt_lup <- rprt_lup  %>%
      tibble::add_case(rprt_nms_chr = "AAA_SUPLRY_ANLYS_MTH",
                       title_chr = "Report outlining the algorithm to run the supplemenatary analysis.",
                       paths_to_rmd_dir_1L_chr = NA_character_,
                       pkg_dirs_chr = "Markdown",
                       packages_chr = "TTU",
                       nms_of_rmd_chr = "Supplement.Rmd") %>%
      dplyr::filter(rprt_nms_chr != "AAA_PMRY_ANLYS_MTH")
  }
  if(add_sharing_rprt_1L_lgl){
    rprt_lup <- rprt_lup  %>%
      tibble::add_case(rprt_nms_chr = "AAA_SHARING_MTH",
                       title_chr = "Supplementary report outlining the algorithm to create and disseminate shareable study output.",
                       paths_to_rmd_dir_1L_chr = NA_character_,
                       pkg_dirs_chr = "Markdown",
                       packages_chr = "TTU",
                       nms_of_rmd_chr = "Share.Rmd")
  }
  if(!is.null(start_at_int[1])){
    rprt_lup <- dplyr::mutate(rprt_lup,
                              title_chr = dplyr::case_when(rprt_nms_chr %in% c("AAA_PMRY_ANLYS_MTH") ~ paste0("Methods Report ",
                                                                                                              start_at_int[1],
                                                                                                              ": Analysis Program (",
                                                                                                              "Primary Analysis",
                                                                                                              ")"),
                                                           rprt_nms_chr %in% c("AAA_SUPLRY_ANLYS_MTH") ~ paste0("Methods Report ",
                                                                                     start_at_int[1]+3,
                                                                                     ": Analysis Program (",
                                                                                     "Secondary Analysis",
                                                                                     ")"),
                                                           rprt_nms_chr %in% c("AAA_RPRT_WRTNG_MTH") ~ paste0("Methods Report ",
                                                                                                           start_at_int[1] + 1,
                                                                                                            ": Reporting Program"),
                                                           rprt_nms_chr %in% c("AAA_SHARING_MTH") ~ paste0("Methods Report ",
                                                                                                           start_at_int[1] + 2,
                                                                                                           ": Sharing Program"),
                                                           rprt_nms_chr %in% c("AAA_TTU_MDL_CTG") ~ paste0("Results Report ",
                                                                                                            ifelse(is.null(reference_1L_int),
                                                                                                                   start_at_int[2],
                                                                                                                   start_at_int[2]+reference_1L_int),
                                                                                                            ": Catalogue of longitudinal models (",
                                                                                                            ifelse(is.null(reference_1L_int),
                                                                                                                   "Primary Analysis",
                                                                                                                   paste0("Secondary Analysis ",LETTERS[reference_1L_int])),
                                                                                                            ")"),
                                                           T ~ title_chr))
  }
  if(!is.null(reference_1L_int)){
    rprt_lup <- dplyr::mutate(rprt_lup,
                              rprt_nms_chr = dplyr::case_when(rprt_nms_chr %in% c("AAA_TTU_MDL_CTG") ~
                                                                paste0("AAA_TTU_MDL_CTG",
                                                                       ifelse(is.null(reference_1L_int),
                                                                              "",
                                                                              ifelse(reference_1L_int == 0,
                                                                                     "",
                                                                                     paste0("-",reference_1L_int)))),
                                                              T ~ rprt_nms_chr))
  }
  return(rprt_lup)
}
# transform_tb_for_merged_col_1 <- function(df,
#                                           output_type_1L_chr = "PDF"){
#   df[[1]] <- as.character(df[[1]])
#   rle.lengths <- rle(df[[1]])$lengths
#   first <- !duplicated(df[[1]])
#   df[[1]][!first] <- ""
#   if(output_type_1L_chr == "PDF")
#     df[[1]][first] <- paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", df[[1]][first], "}}")
#   return(df)
# }
transform_tb_to_mdl_inp <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr,
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round",
    round_bl_val_1L_chr = "Baseline", drop_all_msng_1L_lgl = T, scaling_fctr_dbl = 1,
    tfmn_1L_chr = "NTF", ungroup_1L_lgl = F)
{
    if(length(scaling_fctr_dbl)!=length(predr_vars_nms_chr)){
      scaling_fctr_dbl <- rep(scaling_fctr_dbl[1],length(predr_vars_nms_chr))
    }
  data_tb <- data.frame(data_tb) %>%
    ready4use::remove_labels_from_ds()
    tfd_for_mdl_inp_tb <- data_tb %>% dplyr::select(dplyr::all_of(id_var_nm_1L_chr), dplyr::all_of(round_var_nm_1L_chr),
        dplyr::all_of(predr_vars_nms_chr),
        dplyr::all_of(depnt_var_nm_1L_chr) # Moved from first var in tb
        ) %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>%
        dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr), !!rlang::sym(round_var_nm_1L_chr))
    tfd_for_mdl_inp_tb <- purrr::reduce(1:length(predr_vars_nms_chr),
                                        .init = tfd_for_mdl_inp_tb,
                                        ~ {
                                          idx_1L_int <- as.integer(.y)
                                          .x %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr[idx_1L_int]),
                                                                             .fns = list(baseline = ~dplyr::first(.)*scaling_fctr_dbl[idx_1L_int],
                                                                                         change = ~ifelse(!!rlang::sym(round_var_nm_1L_chr) == round_bl_val_1L_chr,
                                                                                                          0,
                                                                                                          (. - dplyr::lag(.))*scaling_fctr_dbl[idx_1L_int]))))
                                          })
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
        add_tfd_var_to_ds(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                           tfmn_1L_chr = tfmn_1L_chr,
                           depnt_var_max_val_1L_dbl = 0.999)
        if(drop_all_msng_1L_lgl){
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
        stats::na.omit()
    }
    if(ungroup_1L_lgl){
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
        dplyr::ungroup()
    }
      # rename_tb <- tibble::tibble(old_id_xx = tfd_for_mdl_inp_tb %>% dplyr::pull(id_var_nm_1L_chr) %>% unique(),
      #                             new_id_int = 1:length(tfd_for_mdl_inp_tb %>% dplyr::pull(id_var_nm_1L_chr) %>% unique()))
      tfd_for_mdl_inp_tb  <- tfd_for_mdl_inp_tb %>%
        transform_uid_var(id_var_nm_1L_chr = id_var_nm_1L_chr)
        # dplyr::mutate(!!rlang::sym(id_var_nm_1L_chr) := !!rlang::sym(id_var_nm_1L_chr) %>% purrr::map_int(~ ready4::get_from_lup_obj(rename_tb,
        #                                                                                                                   match_value_xx = .x,
        #                                                                                                                   match_var_nm_1L_chr = "old_id_xx",
        #                                                                                                                   target_var_nm_1L_chr = "new_id_int",
        #                                                                                                                   evaluate_1L_lgl = F)))
    return(tfd_for_mdl_inp_tb)
}
transform_tbl_to_rnd_vars <- function(ds_tb,
                                      nbr_of_digits_1L_int = 2L){
  numeric_vars_chr <- ds_tb %>% dplyr::select(where(is.numeric)) %>% names()
  tfd_ds_tb <- ds_tb %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x,nbr_of_digits_1L_int) %>%
                                  format(nsmall = nbr_of_digits_1L_int)))
  return(tfd_ds_tb)
}
transform_timepoint_vals <- function(timepoint_vals_chr,
                                     timepoint_levels_chr,
                                     bl_val_1L_chr){
  if(length(timepoint_vals_chr)==1){
    timepoint_vals_chr <- bl_val_1L_chr
  }else{
    unique_vals_chr <- unique(timepoint_vals_chr)
    if(length(timepoint_vals_chr) >  length(unique_vals_chr))
      timepoint_vals_chr <- c(unique_vals_chr,
                              setdiff(c(bl_val_1L_chr,
                                        setdiff(timepoint_levels_chr,
                                                bl_val_1L_chr)),
                                      unique_vals_chr)[1:(length(timepoint_vals_chr) - length(unique_vals_chr))])
  }
  return(timepoint_vals_chr)
}
transform_ts_mdl_data <- function (mdl_ls, data_tb, depnt_var_nm_1L_chr = "utl_total_w",
                                   predr_vars_nms_chr, id_var_nm_1L_chr = "fkClientID", mdl_nm_1L_chr)
{
  old_data_tb <- data_tb %>% dplyr::select(c(dplyr::all_of(id_var_nm_1L_chr),
                                             dplyr::all_of(depnt_var_nm_1L_chr), predr_vars_nms_chr %>%
                                               purrr::map(~paste0(.x, c("", "_baseline", "_change"))) %>%
                                               purrr::flatten_chr()))
  cnfdl_mdl_ls <- mdl_ls
  cnfdl_mdl_ls$data <- old_data_tb %>% as.data.frame() %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~sample(.x,
                                                                1)))
  return(cnfdl_mdl_ls)
}
transform_uid_var <- function(data_tb,
                              id_var_nm_1L_chr,
                              rename_tb = NULL,
                              old_new_chr = c("old_id_xx","new_id_int")){
  if(is.null(rename_tb)){
    rename_tb <- make_uid_rename_lup(data_tb,
                                     id_var_nm_1L_chr = id_var_nm_1L_chr)
  }
  if(!identical(rename_tb$old_id_xx,rename_tb$new_id_int)){
    fn <- ifelse("character" %in% class(rename_tb %>% dplyr::pull(old_new_chr[2])),
                 purrr::flatten_chr,
                 ifelse("integer" %in% class(rename_tb %>% dplyr::pull(old_new_chr[2])),
                        purrr::flatten_int,
                        purrr::flatten_dbl))
    # first_fn <- ifelse("character" %in% class(rename_tb %>% dplyr::pull(old_new_chr[2])),
    #                    as.character,
    #                    as.numeric)
    tfd_data_tb <- data_tb %>%
      dplyr::mutate(`:=`(!!rlang::sym(id_var_nm_1L_chr),
                         !!rlang::sym(id_var_nm_1L_chr) %>%
                           purrr::map(~ready4::get_from_lup_obj(rename_tb,
                                                                   match_value_xx = .x,
                                                                   match_var_nm_1L_chr = old_new_chr[1],
                                                                   target_var_nm_1L_chr = old_new_chr[2],
                                                                   evaluate_1L_lgl = F)) %>% fn()))
  }else{
    tfd_data_tb <- data_tb
  }
  return(tfd_data_tb)
}
