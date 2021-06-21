#' Transform character vector digit pairs
#' @description transform_chr_digit_pairs() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform character vector digit pairs. Function argument digit_pairs_chr specifies the object to be updated. Argument nbr_of_digits_1L_int provides the object to be updated. The function returns Transformed digit pairs (a character vector).
#' @param digit_pairs_chr Digit pairs (a character vector)
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @return Transformed digit pairs (a character vector)
#' @rdname transform_chr_digit_pairs
#' @export 
#' @importFrom purrr map_chr pluck
#' @keywords internal
transform_chr_digit_pairs <- function (digit_pairs_chr, nbr_of_digits_1L_int = 2L) 
{
    tfd_digit_pairs_chr <- digit_pairs_chr %>% purrr::map_chr(~{
        abs_vals_elmnts_chr <- .x %>% regmatches(gregexpr("[[:digit:]]+", 
            .)) %>% purrr::pluck(1)
        abs_vals_chr <- c(paste0(abs_vals_elmnts_chr[1:2], collapse = "."), 
            paste0(abs_vals_elmnts_chr[3:4], collapse = "."))
        abs_vals_chr[1] <- ifelse(startsWith(.x, paste0("-", 
            abs_vals_chr[1])), paste0("-", abs_vals_chr[1]), 
            abs_vals_chr[1])
        abs_vals_chr[2] <- ifelse(endsWith(.x, paste0("-", abs_vals_chr[2])), 
            paste0("-", abs_vals_chr[2]), abs_vals_chr[2])
        as.numeric(abs_vals_chr) %>% round(digits = nbr_of_digits_1L_int) %>% 
            format(nsmall = nbr_of_digits_1L_int) %>% paste0(collapse = ", ")
    })
    return(tfd_digit_pairs_chr)
}
#' Transform data tibble for comparison
#' @description transform_data_tb_for_cmprsn() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform data tibble for comparison. Function argument data_tb specifies the object to be updated. Argument model_mdl provides the object to be updated. The function returns Transformed data (a tibble).
#' @param data_tb Data (a tibble)
#' @param model_mdl Model (a model)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param source_data_nm_1L_chr Source data name (a character vector of length one), Default: 'Original'
#' @param new_data_is_1L_chr New data is (a character vector of length one), Default: 'Predicted'
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @param family_1L_chr Family (a character vector of length one), Default: 'NA'
#' @param impute_1L_lgl Impute (a logical vector of length one), Default: F
#' @param is_brms_mdl_1L_lgl Is bayesian regression models model (a logical vector of length one), Default: F
#' @param sd_dbl Standard deviation (a double vector), Default: NA
#' @param tfmn_for_bnml_1L_lgl Transformation for binomial (a logical vector of length one), Default: F
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param utl_cls_fn Utility class (a function), Default: NULL
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: NA
#' @return Transformed data (a tibble)
#' @rdname transform_data_tb_for_cmprsn
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @keywords internal
transform_data_tb_for_cmprsn <- function (data_tb, model_mdl, depnt_var_nm_1L_chr = "utl_total_w", 
    source_data_nm_1L_chr = "Original", new_data_is_1L_chr = "Predicted", 
    predn_type_1L_chr = NULL, family_1L_chr = NA_character_, 
    impute_1L_lgl = F, is_brms_mdl_1L_lgl = F, sd_dbl = NA_real_, 
    tfmn_for_bnml_1L_lgl = F, tfmn_1L_chr = "NTF", utl_cls_fn = NULL, 
    utl_min_val_1L_dbl = NA_real_) 
{
    new_data_dbl <- predict_utility(data_tb = data_tb, tfmn_1L_chr = tfmn_1L_chr, 
        model_mdl = model_mdl, force_min_max_1L_lgl = !is.na(utl_min_val_1L_dbl), 
        force_new_data_1L_lgl = T, utl_min_val_1L_dbl = utl_min_val_1L_dbl, 
        impute_1L_lgl = impute_1L_lgl, utl_cls_fn = utl_cls_fn, 
        new_data_is_1L_chr = new_data_is_1L_chr, predn_type_1L_chr = predn_type_1L_chr, 
        sd_dbl = sd_dbl, tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, 
        family_1L_chr = family_1L_chr, is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl)
    tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_predd_var_nm(new_data_is_1L_chr, 
        utl_min_val_1L_dbl = utl_min_val_1L_dbl)), new_data_dbl), 
        `:=`(!!rlang::sym(source_data_nm_1L_chr), !!rlang::sym(depnt_var_nm_1L_chr)))
    return(tfd_data_tb)
}
#' Transform dependent variable name
#' @description transform_depnt_var_nm() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dependent variable name. Function argument depnt_var_nm_1L_chr specifies the object to be updated. Argument tfmn_1L_chr provides the object to be updated. The function returns Transformed dependent variable name (a character vector of length one).
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one)
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @return Transformed dependent variable name (a character vector of length one)
#' @rdname transform_depnt_var_nm
#' @export 

#' @keywords internal
transform_depnt_var_nm <- function (depnt_var_nm_1L_chr, tfmn_1L_chr = "NTF") 
{
    tfd_depnt_var_nm_1L_chr <- paste0(depnt_var_nm_1L_chr, ifelse(tfmn_1L_chr == 
        "NTF", "", paste0("_", tfmn_1L_chr)))
    return(tfd_depnt_var_nm_1L_chr)
}
#' Transform dataset for modelling
#' @description transform_ds_for_mdlng() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dataset for modelling. Function argument data_tb specifies the object to be updated. Argument depnt_var_nm_1L_chr provides the object to be updated. The function returns Transformed data (a tibble).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param covar_var_nms_chr Covariate variable names (a character vector), Default: 'NA'
#' @return Transformed data (a tibble)
#' @rdname transform_ds_for_mdlng
#' @export 
#' @importFrom purrr discard
#' @importFrom tidyr drop_na
#' @importFrom rlang syms
#' @importFrom dplyr select
#' @keywords internal
transform_ds_for_mdlng <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_var_nm_1L_chr, 
    covar_var_nms_chr = NA_character_) 
{
    mdl_vars_chr <- c(names(data_tb)[names(data_tb) %>% startsWith(depnt_var_nm_1L_chr)], 
        predr_var_nm_1L_chr, covar_var_nms_chr) %>% purrr::discard(is.na)
    tfd_data_tb <- data_tb %>% tidyr::drop_na(!!!rlang::syms(mdl_vars_chr)) %>% 
        dplyr::select(!!!rlang::syms(mdl_vars_chr))
    return(tfd_data_tb)
}
#' Transform model variables with classes
#' @description transform_mdl_vars_with_clss() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform model variables with classes. Function argument ds_tb specifies the object to be updated. Argument predictors_lup provides the object to be updated. The function returns Transformed dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param predictors_lup Predictors (a lookup table), Default: NULL
#' @param prototype_lup Prototype (a lookup table), Default: NULL
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param class_fn_1L_chr Class function (a character vector of length one), Default: 'as.numeric'
#' @return Transformed dataset (a tibble)
#' @rdname transform_mdl_vars_with_clss
#' @export 
#' @importFrom tibble add_case
#' @importFrom purrr reduce
#' @importFrom Hmisc label
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom dplyr mutate
#' @importFrom rlang sym exec
transform_mdl_vars_with_clss <- function (ds_tb, predictors_lup = NULL, prototype_lup = NULL, 
    depnt_var_nm_1L_chr = "utl_total_w", class_fn_1L_chr = "as.numeric") 
{
    if (is.null(predictors_lup)) 
        data("predictors_lup", package = "youthvars", envir = environment())
    if (is.null(prototype_lup)) 
        data("prototype_lup", package = "TTU", envir = environment())
    predictors_lup <- tibble::add_case(predictors_lup, short_name_chr = depnt_var_nm_1L_chr, 
        class_chr = "numeric", class_fn_chr = class_fn_1L_chr)
    tfd_ds_tb <- purrr::reduce(predictors_lup$short_name_chr, 
        .init = ds_tb, ~if (.y %in% names(.x)) {
            label_1L_chr <- Hmisc::label(.x[[.y]])
            class_1L_chr <- ready4fun::get_from_lup_obj(predictors_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .y, 
                target_var_nm_1L_chr = "class_chr", evaluate_lgl = F)
            ns_1L_chr <- ready4fun::get_from_lup_obj(prototype_lup, 
                match_var_nm_1L_chr = "type_chr", match_value_xx = class_1L_chr, 
                target_var_nm_1L_chr = "pt_ns_chr", evaluate_lgl = F)
            ns_and_ext_1L_chr <- ifelse(ns_1L_chr == "base", 
                "", paste0(ns_1L_chr, "::"))
            fn <- ifelse(exists(paste0("as.", class_1L_chr), 
                where = paste0("package:", ns_1L_chr)), eval(parse(text = paste0(ns_and_ext_1L_chr, 
                "as.", class_1L_chr))), ifelse(exists(paste0("as_", 
                class_1L_chr), where = paste0("package:", ns_1L_chr)), 
                eval(parse(text = paste0(ns_and_ext_1L_chr, "as_", 
                  class_1L_chr))), eval(parse(text = paste0(ns_and_ext_1L_chr, 
                  class_1L_chr)))))
            tb <- .x %>% dplyr::mutate(`:=`(!!rlang::sym(.y), 
                rlang::exec(ready4fun::get_from_lup_obj(predictors_lup, 
                  match_var_nm_1L_chr = "short_name_chr", match_value_xx = .y, 
                  target_var_nm_1L_chr = "class_fn_chr", evaluate_lgl = T), 
                  !!rlang::sym(.y) %>% fn)))
            if (label_1L_chr != "") {
                Hmisc::label(tb[[.y]]) <- label_1L_chr
            }
            tb
        }
        else {
            .x
        })
    return(tfd_ds_tb)
}
#' Transform params list from
#' @description transform_params_ls_from_lup() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform params list from lookup table. Function argument params_ls specifies the object to be updated. Argument rename_lup provides the object to be updated. The function returns Params (a list).
#' @param params_ls Params (a list)
#' @param rename_lup Rename (a lookup table)
#' @return Params (a list)
#' @rdname transform_params_ls_from_lup
#' @export 
#' @importFrom purrr map_chr
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
transform_params_ls_from_lup <- function (params_ls, rename_lup) 
{
    if (!is.null(params_ls$ds_descvs_ls)) {
        params_ls$ds_descvs_ls$candidate_predrs_chr <- params_ls$ds_descvs_ls$candidate_predrs_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4fun::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_lgl = F)))
        params_ls$ds_descvs_ls$cohort_descv_var_nms_chr <- params_ls$ds_descvs_ls$cohort_descv_var_nms_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4fun::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_lgl = F)))
    }
    if (!is.null(params_ls$predictors_lup)) {
        params_ls$predictors_lup$short_name_chr <- params_ls$predictors_lup$short_name_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4fun::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_lgl = F)))
    }
    params_ls$candidate_covar_nms_chr <- params_ls$candidate_covar_nms_chr %>% 
        purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
            .x, ready4fun::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                evaluate_lgl = F)))
    if (!is.na(params_ls$prefd_covars_chr)) {
        params_ls$prefd_covars_chr <- params_ls$prefd_covars_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4fun::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_lgl = F)))
    }
    if (!is.null(params_ls$candidate_predrs_chr)) {
        params_ls$candidate_predrs_chr <- params_ls$candidate_predrs_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4fun::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_lgl = F)))
    }
    return(params_ls)
}
#' Transform params list to valid
#' @description transform_params_ls_to_valid() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform params list to valid. Function argument params_ls specifies the object to be updated. Argument scndry_analysis_extra_vars_chr provides the object to be updated. The function returns Valid params (a list of lists).
#' @param params_ls Params (a list)
#' @param scndry_analysis_extra_vars_chr Scndry analysis extra variables (a character vector), Default: 'NA'
#' @return Valid params (a list of lists)
#' @rdname transform_params_ls_to_valid
#' @export 
#' @importFrom purrr discard map_chr
#' @importFrom stringi stri_replace_last_fixed stri_replace_all_fixed
#' @importFrom tibble tibble
#' @importFrom dplyr filter rename_with mutate
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom Hmisc label
#' @keywords internal
transform_params_ls_to_valid <- function (params_ls, scndry_analysis_extra_vars_chr = NA_character_) 
{
    target_var_nms_chr <- c(params_ls$ds_descvs_ls$candidate_predrs_chr, 
        params_ls$candidate_covar_nms_chr, scndry_analysis_extra_vars_chr) %>% 
        purrr::discard(is.na) %>% unique()
    valid_var_nms_chr <- target_var_nms_chr %>% stringi::stri_replace_last_fixed("_dbl", 
        "") %>% stringi::stri_replace_last_fixed("_int", "") %>% 
        stringi::stri_replace_all_fixed("_", "")
    unchanged_var_nms_chr <- setdiff(params_ls$ds_descvs_ls$dictionary_tb$var_nm_chr, 
        target_var_nms_chr)
    rename_lup <- tibble::tibble(old_nms_chr = c(unchanged_var_nms_chr, 
        target_var_nms_chr), new_nms_chr = make.unique(c(unchanged_var_nms_chr, 
        valid_var_nms_chr), sep = "V")) %>% dplyr::filter(!old_nms_chr %in% 
        unchanged_var_nms_chr)
    params_ls$ds_tb <- dplyr::rename_with(params_ls$ds_tb, .cols = target_var_nms_chr, 
        ~ready4fun::get_from_lup_obj(rename_lup, match_value_xx = .x, 
            match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
            evaluate_lgl = F))
    var_lbl_1L_chr <- Hmisc::label(params_ls$ds_descvs_ls$dictionary_tb$var_nm_chr)
    params_ls$ds_descvs_ls$dictionary_tb <- params_ls$ds_descvs_ls$dictionary_tb %>% 
        dplyr::mutate(var_nm_chr = var_nm_chr %>% purrr::map_chr(~ifelse(.x %in% 
            rename_lup$old_nms_chr, ready4fun::get_from_lup_obj(rename_lup, 
            match_value_xx = .x, match_var_nm_1L_chr = "old_nms_chr", 
            target_var_nm_1L_chr = "new_nms_chr", evaluate_lgl = F), 
            .x)))
    Hmisc::label(params_ls$ds_descvs_ls$dictionary_tb[["var_nm_chr"]]) <- var_lbl_1L_chr
    valid_params_ls_ls <- list(params_ls = params_ls %>% transform_params_ls_from_lup(rename_lup = rename_lup), 
        rename_lup = rename_lup)
    return(valid_params_ls_ls)
}
#' Transform paths list for scndry
#' @description transform_paths_ls_for_scndry() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform paths list for scndry. Function argument paths_ls specifies the object to be updated. Argument reference_1L_int provides the object to be updated. The function returns Paths (a list).
#' @param paths_ls Paths (a list)
#' @param reference_1L_int Reference (an integer vector of length one), Default: 1
#' @param remove_prmry_1L_lgl Remove prmry (a logical vector of length one), Default: F
#' @return Paths (a list)
#' @rdname transform_paths_ls_for_scndry
#' @export 
#' @importFrom stringr str_sub
transform_paths_ls_for_scndry <- function (paths_ls, reference_1L_int = 1, remove_prmry_1L_lgl = F) 
{
    paths_ls$prmry_analysis_dir_nm_1L_chr <- paths_ls$write_to_dir_nm_1L_chr
    paths_ls$write_to_dir_nm_1L_chr <- paste0(paths_ls$write_to_dir_nm_1L_chr, 
        "/secondary_", reference_1L_int)
    paths_ls$reports_dir_1L_chr <- paste0(paths_ls$reports_dir_1L_chr %>% 
        stringr::str_sub(end = -(nchar(paths_ls$prmry_analysis_dir_nm_1L_chr) + 
            10)), "/", paths_ls$write_to_dir_nm_1L_chr, "/Reports")
    if (remove_prmry_1L_lgl) 
        paths_ls <- paths_ls[names(paths_ls) != "prmry_analysis_dir_nm_1L_chr"]
    return(paths_ls)
}
#' Transform predicted variable name
#' @description transform_predd_var_nm() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform predicted variable name. Function argument new_data_is_1L_chr specifies the object to be updated. Argument utl_min_val_1L_dbl provides the object to be updated. The function returns Tfmd predicted variable name (a character vector of length one).
#' @param new_data_is_1L_chr New data is (a character vector of length one)
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: NA
#' @return Tfmd predicted variable name (a character vector of length one)
#' @rdname transform_predd_var_nm
#' @export 

#' @keywords internal
transform_predd_var_nm <- function (new_data_is_1L_chr, utl_min_val_1L_dbl = NA_real_) 
{
    tfmd_predd_var_nm_1L_chr <- paste0(new_data_is_1L_chr, ifelse(!is.na(utl_min_val_1L_dbl), 
        " (constrained)", ""))
    return(tfmd_predd_var_nm_1L_chr)
}
#' Transform predictor name part of phrases
#' @description transform_predr_nm_part_of_phrases() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform predictor name part of phrases. Function argument phrases_chr specifies the object to be updated. Argument old_nms_chr provides the object to be updated. The function returns Transformed phrases (a character vector).
#' @param phrases_chr Phrases (a character vector)
#' @param old_nms_chr Old names (a character vector), Default: NULL
#' @param new_nms_chr New names (a character vector), Default: NULL
#' @return Transformed phrases (a character vector)
#' @rdname transform_predr_nm_part_of_phrases
#' @export 
#' @importFrom tibble tibble
#' @importFrom purrr map_chr map_lgl
#' @importFrom stringr str_detect str_replace
#' @keywords internal
transform_predr_nm_part_of_phrases <- function (phrases_chr, old_nms_chr = NULL, new_nms_chr = NULL) 
{
    if (is.null(old_nms_chr)) {
        tfd_phrases_chr <- phrases_chr
    }
    else {
        nm_changes_lup_tb = tibble::tibble(old_nms_chr = old_nms_chr, 
            new_nms_chr = new_nms_chr)
        tfd_phrases_chr <- phrases_chr %>% purrr::map_chr(~{
            phrase_1L_chr <- .x
            match_lgl <- nm_changes_lup_tb$old_nms_chr %>% purrr::map_lgl(~stringr::str_detect(phrase_1L_chr, 
                .x))
            if (any(match_lgl)) {
                stringr::str_replace(phrase_1L_chr, nm_changes_lup_tb$old_nms_chr[match_lgl], 
                  nm_changes_lup_tb$new_nms_chr[match_lgl])
            }
            else {
                phrase_1L_chr
            }
        })
    }
    return(tfd_phrases_chr)
}
#' Transform report
#' @description transform_rprt_lup() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform report lookup table. Function argument rprt_lup specifies the object to be updated. Argument add_suplry_rprt_1L_lgl provides the object to be updated. The function returns Report (a lookup table).
#' @param rprt_lup Report (a lookup table)
#' @param add_suplry_rprt_1L_lgl Add suplry report (a logical vector of length one), Default: T
#' @param add_sharing_rprt_1L_lgl Add sharing report (a logical vector of length one), Default: F
#' @return Report (a lookup table)
#' @rdname transform_rprt_lup
#' @export 
#' @importFrom tibble add_case
#' @importFrom dplyr filter
transform_rprt_lup <- function (rprt_lup, add_suplry_rprt_1L_lgl = T, add_sharing_rprt_1L_lgl = F) 
{
    if (add_suplry_rprt_1L_lgl) {
        rprt_lup <- rprt_lup %>% tibble::add_case(rprt_nms_chr = "Suplry_Analysis_Rprt", 
            title_chr = "Report outlining the algorithm to run the supplemenatary analysis.", 
            paths_to_rmd_dir_1L_chr = NA_character_, pkg_dirs_chr = "Markdown", 
            packages_chr = "TTU", nms_of_rmd_chr = "Supplement.Rmd") %>% 
            dplyr::filter(rprt_nms_chr != "Main_Analysis_Rprt")
    }
    if (add_sharing_rprt_1L_lgl) {
        rprt_lup <- rprt_lup %>% tibble::add_case(rprt_nms_chr = "Share_Outp_Rprt", 
            title_chr = "Supplementary report outlining the algorithm to create and disseminate shareable study output.", 
            paths_to_rmd_dir_1L_chr = NA_character_, pkg_dirs_chr = "Markdown", 
            packages_chr = "TTU", nms_of_rmd_chr = "Share.Rmd")
    }
    return(rprt_lup)
}
#' Transform tibble to model input
#' @description transform_tb_to_mdl_inp() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform tibble to model input. Function argument data_tb specifies the object to be updated. Argument depnt_var_nm_1L_chr provides the object to be updated. The function returns Transformed for model input (a tibble).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one), Default: 'Baseline'
#' @param drop_all_msng_1L_lgl Drop all missing (a logical vector of length one), Default: T
#' @param scaling_fctr_dbl Scaling factor (a double vector), Default: 0.01
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param ungroup_1L_lgl Ungroup (a logical vector of length one), Default: F
#' @return Transformed for model input (a tibble)
#' @rdname transform_tb_to_mdl_inp
#' @export 
#' @importFrom dplyr select all_of group_by arrange mutate across first lag ungroup
#' @importFrom rlang sym
#' @importFrom purrr reduce
#' @importFrom stats na.omit
transform_tb_to_mdl_inp <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr, 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline", drop_all_msng_1L_lgl = T, 
    scaling_fctr_dbl = 0.01, tfmn_1L_chr = "NTF", ungroup_1L_lgl = F) 
{
    if (length(scaling_fctr_dbl) != length(predr_vars_nms_chr)) {
        scaling_fctr_dbl <- rep(scaling_fctr_dbl[1], length(predr_vars_nms_chr))
    }
    data_tb <- data.frame(data_tb)
    tfd_for_mdl_inp_tb <- data_tb %>% dplyr::select(dplyr::all_of(id_var_nm_1L_chr), 
        dplyr::all_of(round_var_nm_1L_chr), dplyr::all_of(predr_vars_nms_chr), 
        dplyr::all_of(depnt_var_nm_1L_chr)) %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>% 
        dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr), !!rlang::sym(round_var_nm_1L_chr))
    tfd_for_mdl_inp_tb <- purrr::reduce(1:length(predr_vars_nms_chr), 
        .init = tfd_for_mdl_inp_tb, ~{
            idx_1L_int <- as.integer(.y)
            .x %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr[idx_1L_int]), 
                .fns = list(baseline = ~dplyr::first(.) * scaling_fctr_dbl[idx_1L_int], 
                  change = ~ifelse(!!rlang::sym(round_var_nm_1L_chr) == 
                    round_bl_val_1L_chr, 0, (. - dplyr::lag(.)) * 
                    scaling_fctr_dbl[idx_1L_int]))))
        })
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_depnt_var_nm(depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)), !!rlang::sym(depnt_var_nm_1L_chr) %>% 
        calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, tfmn_is_outp_1L_lgl = F, 
            dep_var_max_val_1L_dbl = 0.999)))
    if (drop_all_msng_1L_lgl) {
        tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% stats::na.omit()
    }
    if (ungroup_1L_lgl) {
        tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::ungroup()
    }
    return(tfd_for_mdl_inp_tb)
}
#' Transform table to rnd variables
#' @description transform_tbl_to_rnd_vars() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform table to rnd variables. Function argument ds_tb specifies the object to be updated. Argument nbr_of_digits_1L_int provides the object to be updated. The function returns Transformed dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @return Transformed dataset (a tibble)
#' @rdname transform_tbl_to_rnd_vars
#' @export 
#' @importFrom dplyr select mutate across
#' @importFrom tibble as_tibble
#' @keywords internal
transform_tbl_to_rnd_vars <- function (ds_tb, nbr_of_digits_1L_int = 2L) 
{
    numeric_vars_chr <- ds_tb %>% dplyr::select(where(is.numeric)) %>% 
        names()
    tfd_ds_tb <- ds_tb %>% tibble::as_tibble() %>% dplyr::mutate(dplyr::across(where(is.numeric), 
        ~round(.x, nbr_of_digits_1L_int) %>% format(nsmall = nbr_of_digits_1L_int)))
    return(tfd_ds_tb)
}
#' Transform time series model data
#' @description transform_ts_mdl_data() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform time series model data. Function argument mdl_ls specifies the object to be updated. Argument data_tb provides the object to be updated. The function returns Cnfdl (a list of models).
#' @param mdl_ls Model list (a list of models)
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @return Cnfdl (a list of models)
#' @rdname transform_ts_mdl_data
#' @export 
#' @importFrom dplyr select all_of summarise across everything
#' @importFrom purrr map flatten_chr
#' @keywords internal
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
