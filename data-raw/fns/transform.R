transform_chr_digit_pairs <- function(digit_pairs_chr,
                                      nbr_of_digits_1L_int = 2L){
  tfd_digit_pairs_chr <- digit_pairs_chr %>%
    purrr::map_chr(~{
      abs_vals_elmnts_chr <- .x  %>% regmatches(gregexpr("[[:digit:]]+", .)) %>% purrr::pluck(1)
      abs_vals_chr <- c(paste0(abs_vals_elmnts_chr[1:2], collapse = "."),
                        paste0(abs_vals_elmnts_chr[3:4], collapse = "."))
      abs_vals_chr[1] <- ifelse(startsWith(.x,paste0("-",abs_vals_chr[1])),paste0("-",abs_vals_chr[1]),abs_vals_chr[1])
      abs_vals_chr[2] <- ifelse(endsWith(.x,paste0("-",abs_vals_chr[2])),paste0("-",abs_vals_chr[2]),abs_vals_chr[2])
      as.numeric(abs_vals_chr) %>% round(digits=nbr_of_digits_1L_int) %>% format(nsmall=nbr_of_digits_1L_int) %>% paste0(collapse = ", ")
    })
  return(tfd_digit_pairs_chr)
}
transform_data_tb_for_cmprsn <- function (data_tb, model_mdl, depnt_var_nm_1L_chr = "utl_total_w",
    source_data_nm_1L_chr = "Original", tf_type_1L_chr = "Predicted",
    predn_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_, tfmn_1L_chr = "NTF")
{
    if(tf_type_1L_chr == "Predicted")
        new_data_dbl <- stats::predict(model_mdl, type = predn_type_1L_chr)
    if(tf_type_1L_chr == "Simulated"){
      if("betareg" %in% class(model_mdl)){
        new_data_dbl <- rlang::exec(enrichwith::get_simulate_function(model_mdl),
                                    coef(enrichwith::enrich(model_mdl, with = "auxiliary functions")))
      }else{
        if (!tfmn_for_bnml_1L_lgl){
          new_data_dbl <- stats::simulate(model_mdl)$sim_1
        }else{
          new_data_dbl <- (stats::predict(model_mdl) + stats::rnorm(nrow(data_tb),
                                                                    0, stats::sigma(model_mdl)))
        }
      }
    }
        new_data_dbl <- new_data_dbl %>%
          calculate_dpnt_var_tfmn(tfmn_1L_chr = ifelse(tfmn_for_bnml_1L_lgl & tf_type_1L_chr == "Simulated",
                                                       ifelse(family_1L_chr == "quasibinomial(log)",
                                                              "LOG",
                                                              ifelse(family_1L_chr == "quasibinomial(logit)",
                                                                     "LOGIT",
                                                                     ifelse(family_1L_chr == "quasibinomial(cloglog)",
                                                                            "CLL", "NTF"))),
                                                       tfmn_1L_chr),
                                  tfmn_is_outp_1L_lgl = T)
    tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(tf_type_1L_chr),
        new_data_dbl), `:=`(!!rlang::sym(source_data_nm_1L_chr),
        !!rlang::sym(depnt_var_nm_1L_chr)))
    return(tfd_data_tb)
}
transform_depnt_var_nm <- function (depnt_var_nm_1L_chr, tfmn_1L_chr = "NTF")
{
    tfd_depnt_var_nm_1L_chr <- paste0(depnt_var_nm_1L_chr, ifelse(tfmn_1L_chr ==
        "NTF", "", paste0("_", tfmn_1L_chr)))
    return(tfd_depnt_var_nm_1L_chr)
}
transform_depnt_var_nm_for_cll <- function (depnt_var_nm_1L_chr)
{
    tfd_depnt_var_nm_1L_chr <- paste0(depnt_var_nm_1L_chr, "_cloglog")
    return(tfd_depnt_var_nm_1L_chr)
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
transform_mdl_vars_with_clss <- function(ds_tb,
                                         predictors_lup = NULL,
                                         prototype_lup = NULL,
                                         depnt_var_nm_1L_chr = "utl_total_w",
                                         class_fn_1L_chr = "as.numeric"){
  if(is.null(predictors_lup))
    data("predictors_lup", package = "youthvars", envir = environment())
  if(is.null(prototype_lup))
    data("prototype_lup", package = "TTU", envir = environment())
  predictors_lup <- tibble::add_case(predictors_lup,
                                     short_name_chr = depnt_var_nm_1L_chr,
                                     class_chr = "numeric",
                                     class_fn_chr = class_fn_1L_chr)
  tfd_ds_tb <- purrr::reduce(predictors_lup$short_name_chr,
                             .init = ds_tb,
                             ~ if(.y %in% names(.x)){
                               label_1L_chr <- Hmisc::label(.x[[.y]])
                               class_1L_chr <- ready4fun::get_from_lup_obj(predictors_lup,
                                                                           match_var_nm_1L_chr = "short_name_chr",
                                                                           match_value_xx = .y,
                                                                           target_var_nm_1L_chr = "class_chr",
                                                                           evaluate_lgl = F)
                               ns_1L_chr <- ready4fun::get_from_lup_obj(prototype_lup,
                                                                        match_var_nm_1L_chr = "type_chr",
                                                                        match_value_xx = class_1L_chr,
                                                                        target_var_nm_1L_chr = "pt_ns_chr",
                                                                        evaluate_lgl = F)
                               ns_and_ext_1L_chr <- ifelse(ns_1L_chr == "base",
                                                           "",
                                                           paste0(ns_1L_chr,"::"))
                               fn <- ifelse(exists(paste0("as.",class_1L_chr), where = paste0("package:",ns_1L_chr)),
                                            eval(parse(text = paste0(ns_and_ext_1L_chr,"as.",class_1L_chr))),
                                            ifelse(exists(paste0("as_",class_1L_chr), where = paste0("package:",ns_1L_chr)),
                                                   eval(parse(text = paste0(ns_and_ext_1L_chr,"as_",class_1L_chr))),
                                                   eval(parse(text = paste0(ns_and_ext_1L_chr,class_1L_chr)))))
                               tb <- .x %>% dplyr::mutate(!!rlang::sym(.y) := rlang::exec(ready4fun::get_from_lup_obj(predictors_lup,
                                                                                                                      match_var_nm_1L_chr = "short_name_chr",
                                                                                                                      match_value_xx = .y,
                                                                                                                      target_var_nm_1L_chr = "class_fn_chr",
                                                                                                                      evaluate_lgl = T),
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
  params_ls$ds_tb <- dplyr::rename_with(params_ls$ds_tb,
                                        .cols = target_var_nms_chr,
                                        ~ ready4fun::get_from_lup_obj(rename_lup,
                                                                      match_value_xx = .x,
                                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                                      evaluate_lgl = F))
  var_lbl_1L_chr <- Hmisc::label(params_ls$ds_descvs_ls$dictionary_tb$var_nm_chr)
  params_ls$ds_descvs_ls$dictionary_tb <- params_ls$ds_descvs_ls$dictionary_tb %>%
    dplyr::mutate(var_nm_chr = var_nm_chr %>%
                    purrr::map_chr(~ifelse(.x %in% rename_lup$old_nms_chr,
                                           ready4fun::get_from_lup_obj(rename_lup,
                                                                       match_value_xx = .x,
                                                                       match_var_nm_1L_chr = "old_nms_chr",
                                                                       target_var_nm_1L_chr = "new_nms_chr",
                                                                       evaluate_lgl = F),
                                           .x)))
  Hmisc::label(params_ls$ds_descvs_ls$dictionary_tb[["var_nm_chr"]]) <- var_lbl_1L_chr#
  valid_params_ls_ls <- list(params_ls = params_ls %>%
                               transform_params_ls_from_lup(rename_lup = rename_lup),
                             rename_lup = rename_lup)
  return(valid_params_ls_ls)
}
transform_params_ls_from_lup <- function(params_ls,
                                         rename_lup){
  if(!is.null(params_ls$ds_descvs_ls)){
    params_ls$ds_descvs_ls$candidate_predrs_chr <- params_ls$ds_descvs_ls$candidate_predrs_chr %>%
      purrr::map_chr(~ready4fun::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_lgl = F))
    valid_params_ls_ls$params_ls$ds_descvs_ls$cohort_descv_var_nms_chr <- valid_params_ls_ls$params_ls$ds_descvs_ls$cohort_descv_var_nms_chr %>%
      purrr::map_chr(~ready4fun::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_lgl = F))
  }
  if(!is.null(params_ls$predictors_lup)){
    params_ls$predictors_lup$short_name_chr <-  params_ls$predictors_lup$short_name_chr %>%
      purrr::map_chr(~ready4fun::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_lgl = F))
  }
  params_ls$candidate_covar_nms_chr <- params_ls$candidate_covar_nms_chr %>%
    purrr::map_chr(~ready4fun::get_from_lup_obj(rename_lup,
                                                match_value_xx = .x,
                                                match_var_nm_1L_chr = "old_nms_chr",
                                                target_var_nm_1L_chr = "new_nms_chr",
                                                evaluate_lgl = F))
  if(!is.na(params_ls$prefd_covars_chr)){
    params_ls$prefd_covars_chr <- params_ls$prefd_covars_chr %>%
      purrr::map_chr(~ready4fun::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_lgl = F))
  }
  if(!is.null(params_ls$candidate_predrs_chr)){
    params_ls$candidate_predrs_chr <- params_ls$candidate_predrs_chr %>%
      purrr::map_chr(~ready4fun::get_from_lup_obj(rename_lup,
                                                  match_value_xx = .x,
                                                  match_var_nm_1L_chr = "old_nms_chr",
                                                  target_var_nm_1L_chr = "new_nms_chr",
                                                  evaluate_lgl = F))
  }
  return(params_ls)
}
transform_paths_ls_for_scndry <- function(paths_ls,
                                          reference_1L_int = 1,
                                          remove_prmry_1L_lgl = F){
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
  return(paths_ls)
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
                               add_sharing_rprt_1L_lgl = F){
  if(add_suplry_rprt_1L_lgl){
    rprt_lup <- rprt_lup  %>%
      tibble::add_case(rprt_nms_chr = "Suplry_Analysis_Rprt",
                       title_chr = "Report outlining the algorithm to run the supplemenatary analysis.",
                       paths_to_rmd_dir_1L_chr = NA_character_,
                       pkg_dirs_chr = "Markdown",
                       packages_chr = "TTU",
                       nms_of_rmd_chr = "Supplement.Rmd") %>%
      dplyr::filter(rprt_nms_chr != "Main_Analysis_Rprt")
  }
  if(add_sharing_rprt_1L_lgl){
    rprt_lup <- rprt_lup  %>%
      tibble::add_case(rprt_nms_chr = "Share_Outp_Rprt",
                       title_chr = "Supplementary report outlining the algorithm to create and disseminate shareable study output.",
                       paths_to_rmd_dir_1L_chr = NA_character_,
                       pkg_dirs_chr = "Markdown",
                       packages_chr = "TTU",
                       nms_of_rmd_chr = "Share.Rmd")
  }
  return(rprt_lup)
}
transform_tb_to_mdl_inp <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr,
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round",
    round_bl_val_1L_chr = "Baseline", drop_all_msng_1L_lgl = T, scaling_fctr_dbl = 0.01, ungroup_1L_lgl = F,
    add_cll_tfmn_1L_lgl = T)
{
    if(length(scaling_fctr_dbl)!=length(predr_vars_nms_chr)){
      scaling_fctr_dbl <- rep(scaling_fctr_dbl[1],length(predr_vars_nms_chr))
    }
  data_tb <- data.frame(data_tb)
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
    if(add_cll_tfmn_1L_lgl){
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
      dplyr::mutate(`:=`(!!rlang::sym(transform_depnt_var_nm_for_cll(depnt_var_nm_1L_chr)),
        log(-log(1 - ifelse(!!rlang::sym(depnt_var_nm_1L_chr) ==
            1, 0.999, !!rlang::sym(depnt_var_nm_1L_chr))))))
    }
    if(drop_all_msng_1L_lgl){
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
        stats::na.omit()
    }
    if(ungroup_1L_lgl){
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
        dplyr::ungroup()
    }
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
