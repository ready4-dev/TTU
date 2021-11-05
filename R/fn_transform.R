#' Transform character vector digit pairs
#' @description transform_chr_digit_pairs() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform character vector digit pairs. Function argument digit_pairs_chr specifies the object to be updated. Argument nbr_of_digits_1L_int provides the object to be updated. The function returns Transformed digit pairs (a character vector).
#' @param digit_pairs_chr Digit pairs (a character vector)
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @return Transformed digit pairs (a character vector)
#' @rdname transform_chr_digit_pairs
#' @export 
#' @importFrom purrr map_chr pluck
#' @importFrom stringr str_squish
#' @keywords internal
transform_chr_digit_pairs <- function (digit_pairs_chr, nbr_of_digits_1L_int = 2L) 
{
    tfd_digit_pairs_chr <- digit_pairs_chr %>% purrr::map_chr(~{
        abs_vals_chr <- .x %>% strsplit(",") %>% purrr::pluck(1) %>% 
            stringr::str_squish()
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
#' @param sfx_1L_chr Suffix (a character vector of length one), Default: ''
#' @param tfmn_for_bnml_1L_lgl Transformation for binomial (a logical vector of length one), Default: F
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param utl_cls_fn Utility class (a function), Default: NULL
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: NA
#' @return Transformed data (a tibble)
#' @rdname transform_data_tb_for_cmprsn
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
transform_data_tb_for_cmprsn <- function (data_tb, model_mdl, depnt_var_nm_1L_chr = "utl_total_w", 
    source_data_nm_1L_chr = "Original", new_data_is_1L_chr = "Predicted", 
    predn_type_1L_chr = NULL, family_1L_chr = NA_character_, 
    impute_1L_lgl = F, is_brms_mdl_1L_lgl = F, sd_dbl = NA_real_, 
    sfx_1L_chr = "", tfmn_for_bnml_1L_lgl = F, tfmn_1L_chr = "NTF", 
    utl_cls_fn = NULL, utl_min_val_1L_dbl = NA_real_) 
{
    new_data_dbl <- predict_utility(data_tb = data_tb, tfmn_1L_chr = tfmn_1L_chr, 
        model_mdl = model_mdl, force_min_max_1L_lgl = !is.na(utl_min_val_1L_dbl), 
        force_new_data_1L_lgl = T, utl_min_val_1L_dbl = utl_min_val_1L_dbl, 
        impute_1L_lgl = impute_1L_lgl, utl_cls_fn = utl_cls_fn, 
        new_data_is_1L_chr = new_data_is_1L_chr, predn_type_1L_chr = predn_type_1L_chr, 
        sd_dbl = sd_dbl, tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, 
        family_1L_chr = family_1L_chr, is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl)
    tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_predd_var_nm(new_data_is_1L_chr, 
        sfx_1L_chr = sfx_1L_chr, utl_min_val_1L_dbl = utl_min_val_1L_dbl)), 
        new_data_dbl), `:=`(!!rlang::sym(source_data_nm_1L_chr), 
        !!rlang::sym(depnt_var_nm_1L_chr)))
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
#' Transform dictionary with rename lookup table
#' @description transform_dict_with_rename_lup() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dictionary with rename lookup table. Function argument dictionary_tb specifies the object to be updated. Argument rename_lup provides the object to be updated. The function returns Transformed dictionary (a tibble).
#' @param dictionary_tb Dictionary (a tibble)
#' @param rename_lup Rename (a lookup table)
#' @return Transformed dictionary (a tibble)
#' @rdname transform_dict_with_rename_lup
#' @export 
#' @importFrom Hmisc label
#' @importFrom dplyr mutate
#' @importFrom purrr map_chr
#' @importFrom ready4 get_from_lup_obj
#' @keywords internal
transform_dict_with_rename_lup <- function (dictionary_tb, rename_lup) 
{
    var_lbl_1L_chr <- Hmisc::label(dictionary_tb$var_nm_chr)
    tfd_dictionary_tb <- dictionary_tb %>% dplyr::mutate(var_nm_chr = var_nm_chr %>% 
        purrr::map_chr(~ifelse(.x %in% rename_lup$old_nms_chr, 
            ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                evaluate_1L_lgl = F), .x)))
    Hmisc::label(tfd_dictionary_tb[["var_nm_chr"]]) <- var_lbl_1L_chr
    return(tfd_dictionary_tb)
}
#' Transform dataset for all comparison plots
#' @description transform_ds_for_all_cmprsn_plts() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dataset for all comparison plots. Function argument tfd_data_tb specifies the object to be updated. Argument model_mdl provides the object to be updated. The function returns Transformed data (a tibble).
#' @param tfd_data_tb Transformed data (a tibble)
#' @param model_mdl Model (a model)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one)
#' @param is_brms_mdl_1L_lgl Is bayesian regression models model (a logical vector of length one)
#' @param predn_type_1L_chr Prediction type (a character vector of length one)
#' @param sd_dbl Standard deviation (a double vector)
#' @param sfx_1L_chr Suffix (a character vector of length one), Default: ''
#' @param tfmn_1L_chr Transformation (a character vector of length one)
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: -1
#' @return Transformed data (a tibble)
#' @rdname transform_ds_for_all_cmprsn_plts
#' @export 
#' @importFrom dplyr ungroup
#' @keywords internal
transform_ds_for_all_cmprsn_plts <- function (tfd_data_tb, model_mdl, depnt_var_nm_1L_chr, is_brms_mdl_1L_lgl, 
    predn_type_1L_chr, sd_dbl, sfx_1L_chr = "", tfmn_1L_chr, 
    utl_min_val_1L_dbl = -1) 
{
    tfd_data_tb <- transform_data_tb_for_cmprsn(tfd_data_tb %>% 
        dplyr::ungroup(), model_mdl = model_mdl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predn_type_1L_chr = predn_type_1L_chr, sfx_1L_chr = sfx_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr) %>% transform_data_tb_for_cmprsn(model_mdl = model_mdl, 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, family_1L_chr = NA_character_, 
        is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl, new_data_is_1L_chr = "Simulated", 
        predn_type_1L_chr = predn_type_1L_chr, sd_dbl = sd_dbl, 
        sfx_1L_chr = sfx_1L_chr, tfmn_1L_chr = tfmn_1L_chr, tfmn_for_bnml_1L_lgl = FALSE) %>% 
        transform_data_tb_for_cmprsn(model_mdl = model_mdl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
            predn_type_1L_chr = predn_type_1L_chr, sfx_1L_chr = sfx_1L_chr, 
            tfmn_1L_chr = tfmn_1L_chr, utl_min_val_1L_dbl = utl_min_val_1L_dbl) %>% 
        transform_data_tb_for_cmprsn(model_mdl = model_mdl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
            family_1L_chr = NA_character_, is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl, 
            new_data_is_1L_chr = "Simulated", predn_type_1L_chr = predn_type_1L_chr, 
            sfx_1L_chr = sfx_1L_chr, sd_dbl = sd_dbl, tfmn_1L_chr = tfmn_1L_chr, 
            tfmn_for_bnml_1L_lgl = FALSE, utl_min_val_1L_dbl = utl_min_val_1L_dbl)
    return(tfd_data_tb)
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
transform_ds_for_mdlng <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_var_nm_1L_chr, 
    covar_var_nms_chr = NA_character_) 
{
    mdl_vars_chr <- c(names(data_tb)[names(data_tb) %>% startsWith(depnt_var_nm_1L_chr)], 
        predr_var_nm_1L_chr, covar_var_nms_chr) %>% purrr::discard(is.na)
    tfd_data_tb <- data_tb %>% tidyr::drop_na(!!!rlang::syms(mdl_vars_chr)) %>% 
        dplyr::select(!!!rlang::syms(mdl_vars_chr))
    return(tfd_data_tb)
}
#' Transform dataset to prediction dataset
#' @description transform_ds_to_predn_ds() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dataset to prediction dataset. Function argument data_tb specifies the object to be updated. Argument predr_vars_nms_chr provides the object to be updated. The function returns Data (a tibble).
#' @param data_tb Data (a tibble)
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param tfmn_1L_chr Transformation (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one)
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one)
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one)
#' @param predictors_lup Predictors (a lookup table)
#' @return Data (a tibble)
#' @rdname transform_ds_to_predn_ds
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym exec
#' @importFrom purrr reduce map_dbl
#' @importFrom ready4 get_from_lup_obj
transform_ds_to_predn_ds <- function (data_tb, predr_vars_nms_chr, tfmn_1L_chr, depnt_var_nm_1L_chr, 
    id_var_nm_1L_chr, round_var_nm_1L_chr, round_bl_val_1L_chr, 
    predictors_lup) 
{
    data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(depnt_var_nm_1L_chr), 
        NA_real_))
    data_tb <- purrr::reduce(predr_vars_nms_chr, .init = data_tb, 
        ~{
            predr_cls_fn <- eval(parse(text = ready4::get_from_lup_obj(predictors_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .y, 
                target_var_nm_1L_chr = "class_fn_chr", evaluate_1L_lgl = F)))
            dplyr::mutate(.x, `:=`(!!rlang::sym(.y), !!rlang::sym(.y) %>% 
                rlang::exec(.fn = predr_cls_fn)))
        })
    data_tb <- data_tb %>% transform_tb_to_mdl_inp(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr, 
        drop_all_msng_1L_lgl = F, scaling_fctr_dbl = purrr::map_dbl(predr_vars_nms_chr, 
            ~ready4::get_from_lup_obj(predictors_lup, target_var_nm_1L_chr = "mdl_scaling_dbl", 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .x, 
                evaluate_1L_lgl = F)), ungroup_1L_lgl = T, tfmn_1L_chr = tfmn_1L_chr)
    return(data_tb)
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
#' @importFrom ready4 get_rds_from_dv get_from_lup_obj
#' @importFrom tibble add_case
#' @importFrom purrr reduce
#' @importFrom Hmisc label
#' @importFrom dplyr mutate
#' @importFrom rlang sym exec
transform_mdl_vars_with_clss <- function (ds_tb, predictors_lup = NULL, prototype_lup = NULL, 
    depnt_var_nm_1L_chr = "utl_total_w", class_fn_1L_chr = "as.numeric") 
{
    if (is.null(predictors_lup)) 
        data("predictors_lup", package = "youthvars", envir = environment())
    if (is.null(prototype_lup)) 
        prototype_lup <- ready4::get_rds_from_dv("prototype_lup", 
            server_1L_chr = "dataverse.harvard.edu")
    if (!is.null(depnt_var_nm_1L_chr)) {
        predictors_lup <- tibble::add_case(predictors_lup, short_name_chr = depnt_var_nm_1L_chr, 
            class_chr = "numeric", class_fn_chr = class_fn_1L_chr)
    }
    tfd_ds_tb <- purrr::reduce(predictors_lup$short_name_chr, 
        .init = ds_tb, ~if (.y %in% names(.x)) {
            label_1L_chr <- Hmisc::label(.x[[.y]])
            class_1L_chr <- ready4::get_from_lup_obj(predictors_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .y, 
                target_var_nm_1L_chr = "class_chr", evaluate_1L_lgl = F)
            ns_1L_chr <- ready4::get_from_lup_obj(prototype_lup, 
                match_var_nm_1L_chr = "type_chr", match_value_xx = class_1L_chr, 
                target_var_nm_1L_chr = "pt_ns_chr", evaluate_1L_lgl = F)
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
                rlang::exec(ready4::get_from_lup_obj(predictors_lup, 
                  match_var_nm_1L_chr = "short_name_chr", match_value_xx = .y, 
                  target_var_nm_1L_chr = "class_fn_chr", evaluate_1L_lgl = T), 
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
#' Transform names
#' @description transform_names() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform names. Function argument names_chr specifies the object to be updated. Argument rename_lup provides the object to be updated. The function returns New names (a character vector).
#' @param names_chr Names (a character vector)
#' @param rename_lup Rename (a lookup table)
#' @param invert_1L_lgl Invert (a logical vector of length one), Default: F
#' @return New names (a character vector)
#' @rdname transform_names
#' @export 
#' @importFrom purrr map_chr
#' @importFrom ready4 get_from_lup_obj
#' @keywords internal
transform_names <- function (names_chr, rename_lup, invert_1L_lgl = F) 
{
    new_names_chr <- names_chr %>% purrr::map_chr(~ifelse((!invert_1L_lgl & 
        .x %in% rename_lup$old_nms_chr) | (invert_1L_lgl & .x %in% 
        rename_lup$new_nms_chr), .x %>% ready4::get_from_lup_obj(data_lookup_tb = rename_lup, 
        match_var_nm_1L_chr = ifelse(invert_1L_lgl, "new_nms_chr", 
            "old_nms_chr"), target_var_nm_1L_chr = ifelse(invert_1L_lgl, 
            "old_nms_chr", "new_nms_chr"), evaluate_1L_lgl = F), 
        .x))
    return(new_names_chr)
}
#' Transform names in model table
#' @description transform_nms_in_mdl_tbl() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform names in model table. Function argument mdl_tbl_tb specifies the object to be updated. Argument col_nm_1L_chr provides the object to be updated. The function returns Transformed model table (a tibble).
#' @param mdl_tbl_tb Model table (a tibble)
#' @param col_nm_1L_chr Column name (a character vector of length one), Default: 'Parameter'
#' @param var_nm_change_lup Variable name change (a lookup table), Default: NULL
#' @return Transformed model table (a tibble)
#' @rdname transform_nms_in_mdl_tbl
#' @export 
#' @importFrom dplyr mutate case_when
#' @importFrom rlang sym
#' @importFrom purrr map_lgl map_chr pluck
#' @importFrom stringi stri_locate_first_fixed
#' @importFrom stringr str_sub
#' @keywords internal
transform_nms_in_mdl_tbl <- function (mdl_tbl_tb, col_nm_1L_chr = "Parameter", var_nm_change_lup = NULL) 
{
    if (is.null(var_nm_change_lup)) {
        tfd_mdl_tbl_tb <- mdl_tbl_tb
    }
    else {
        tfd_mdl_tbl_tb <- mdl_tbl_tb %>% dplyr::mutate(`:=`(!!rlang::sym(col_nm_1L_chr), 
            dplyr::case_when(!!rlang::sym(col_nm_1L_chr) %>% 
                purrr::map_lgl(~(endsWith(.x, " model") | endsWith(.x, 
                  " baseline") | endsWith(.x, " change"))) ~ 
                !!rlang::sym(col_nm_1L_chr) %>% purrr::map_chr(~{
                  sfx_starts_1L_int <- stringi::stri_locate_first_fixed(.x, 
                    " ")[[1, 1]]
                  paste0(stringr::str_sub(.x, end = (sfx_starts_1L_int - 
                    1)) %>% strsplit("_") %>% purrr::pluck(1) %>% 
                    transform_names(rename_lup = var_nm_change_lup) %>% 
                    paste0(collapse = "_"), stringr::str_sub(.x, 
                    start = sfx_starts_1L_int))
                }), T ~ !!rlang::sym(col_nm_1L_chr))))
    }
    return(tfd_mdl_tbl_tb)
}
#' Transform parameters list from lookup table
#' @description transform_params_ls_from_lup() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform parameters list from lookup table. Function argument params_ls specifies the object to be updated. Argument rename_lup provides the object to be updated. The function returns Parameters (a list).
#' @param params_ls Parameters (a list)
#' @param rename_lup Rename (a lookup table)
#' @return Parameters (a list)
#' @rdname transform_params_ls_from_lup
#' @export 
#' @importFrom purrr map_chr
#' @importFrom ready4 get_from_lup_obj
#' @keywords internal
transform_params_ls_from_lup <- function (params_ls, rename_lup) 
{
    if (!is.null(params_ls$ds_descvs_ls)) {
        params_ls$ds_descvs_ls$candidate_predrs_chr <- params_ls$ds_descvs_ls$candidate_predrs_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_1L_lgl = F)))
        params_ls$ds_descvs_ls$cohort_descv_var_nms_chr <- params_ls$ds_descvs_ls$cohort_descv_var_nms_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_1L_lgl = F)))
    }
    if (!is.null(params_ls$predictors_lup)) {
        params_ls$predictors_lup$short_name_chr <- params_ls$predictors_lup$short_name_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_1L_lgl = F)))
    }
    params_ls$candidate_covar_nms_chr <- params_ls$candidate_covar_nms_chr %>% 
        purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
            .x, ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                evaluate_1L_lgl = F)))
    if (!is.na(params_ls$prefd_covars_chr)) {
        params_ls$prefd_covars_chr <- params_ls$prefd_covars_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_1L_lgl = F)))
    }
    if (!is.null(params_ls$candidate_predrs_chr)) {
        params_ls$candidate_predrs_chr <- params_ls$candidate_predrs_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_1L_lgl = F)))
    }
    return(params_ls)
}
#' Transform parameters list to valid
#' @description transform_params_ls_to_valid() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform parameters list to valid. Function argument params_ls specifies the object to be updated. Argument scndry_analysis_extra_vars_chr provides the object to be updated. The function returns Valid parameters (a list of lists).
#' @param params_ls Parameters (a list)
#' @param scndry_analysis_extra_vars_chr Secondary analysis extra variables (a character vector), Default: 'NA'
#' @return Valid parameters (a list of lists)
#' @rdname transform_params_ls_to_valid
#' @export 
#' @importFrom purrr discard
#' @importFrom stringi stri_replace_last_fixed stri_replace_all_fixed
#' @importFrom tibble tibble
#' @importFrom dplyr filter
#' @importFrom youthvars transform_ds_with_rename_lup
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
    params_ls$ds_tb <- youthvars::transform_ds_with_rename_lup(params_ls$ds_tb, 
        rename_lup = rename_lup, target_var_nms_chr = target_var_nms_chr)
    params_ls$ds_descvs_ls$dictionary_tb <- params_ls$ds_descvs_ls$dictionary_tb %>% 
        transform_dict_with_rename_lup(rename_lup = rename_lup)
    rename_lup <- rename_lup %>% dplyr::filter(old_nms_chr != 
        new_nms_chr)
    valid_params_ls_ls <- list(params_ls = params_ls %>% transform_params_ls_from_lup(rename_lup = rename_lup), 
        rename_lup = rename_lup)
    return(valid_params_ls_ls)
}
#' Transform paths list for secondary
#' @description transform_paths_ls_for_scndry() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform paths list for secondary. Function argument paths_ls specifies the object to be updated. Argument reference_1L_int provides the object to be updated. The function returns Paths (a list).
#' @param paths_ls Paths (a list)
#' @param reference_1L_int Reference (an integer vector of length one), Default: 1
#' @param remove_prmry_1L_lgl Remove primary (a logical vector of length one), Default: F
#' @param remove_mkdn_1L_lgl Remove markdown (a logical vector of length one), Default: F
#' @return Paths (a list)
#' @rdname transform_paths_ls_for_scndry
#' @export 
#' @importFrom stringr str_sub
#' @keywords internal
transform_paths_ls_for_scndry <- function (paths_ls, reference_1L_int = 1, remove_prmry_1L_lgl = F, 
    remove_mkdn_1L_lgl = F) 
{
    paths_ls$prmry_analysis_dir_nm_1L_chr <- paths_ls$write_to_dir_nm_1L_chr
    paths_ls$write_to_dir_nm_1L_chr <- paste0(paths_ls$write_to_dir_nm_1L_chr, 
        "/secondary_", reference_1L_int)
    paths_ls$reports_dir_1L_chr <- paste0(paths_ls$reports_dir_1L_chr %>% 
        stringr::str_sub(end = -(nchar(paths_ls$prmry_analysis_dir_nm_1L_chr) + 
            10)), "/", paths_ls$write_to_dir_nm_1L_chr, "/Reports")
    if (remove_prmry_1L_lgl) 
        paths_ls <- paths_ls[names(paths_ls) != "prmry_analysis_dir_nm_1L_chr"]
    if (remove_mkdn_1L_lgl) 
        paths_ls <- paths_ls[names(paths_ls) != "reports_dir_1L_chr"]
    return(paths_ls)
}
#' Transform predicted variable name
#' @description transform_predd_var_nm() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform predicted variable name. Function argument new_data_is_1L_chr specifies the object to be updated. Argument sfx_1L_chr provides the object to be updated. The function returns Transformed predicted variable name (a character vector of length one).
#' @param new_data_is_1L_chr New data is (a character vector of length one)
#' @param sfx_1L_chr Suffix (a character vector of length one), Default: ''
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: NA
#' @return Transformed predicted variable name (a character vector of length one)
#' @rdname transform_predd_var_nm
#' @export 
#' @keywords internal
transform_predd_var_nm <- function (new_data_is_1L_chr, sfx_1L_chr = "", utl_min_val_1L_dbl = NA_real_) 
{
    tfd_predd_var_nm_1L_chr <- paste0(new_data_is_1L_chr, sfx_1L_chr, 
        ifelse(!is.na(utl_min_val_1L_dbl), " (constrained)", 
            ""))
    return(tfd_predd_var_nm_1L_chr)
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
#' Transform report lookup table
#' @description transform_rprt_lup() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform report lookup table. Function argument rprt_lup specifies the object to be updated. Argument add_suplry_rprt_1L_lgl provides the object to be updated. The function returns Report (a lookup table).
#' @param rprt_lup Report (a lookup table)
#' @param add_suplry_rprt_1L_lgl Add supplementary report (a logical vector of length one), Default: T
#' @param add_sharing_rprt_1L_lgl Add sharing report (a logical vector of length one), Default: F
#' @param start_at_int Start at (an integer vector), Default: NULL
#' @param reference_1L_int Reference (an integer vector of length one), Default: NULL
#' @return Report (a lookup table)
#' @rdname transform_rprt_lup
#' @export 
#' @importFrom tibble add_case
#' @importFrom dplyr filter mutate case_when
#' @keywords internal
transform_rprt_lup <- function (rprt_lup, add_suplry_rprt_1L_lgl = T, add_sharing_rprt_1L_lgl = F, 
    start_at_int = NULL, reference_1L_int = NULL) 
{
    if (add_suplry_rprt_1L_lgl) {
        rprt_lup <- rprt_lup %>% tibble::add_case(rprt_nms_chr = "AAA_SUPLRY_ANLYS_MTH", 
            title_chr = "Report outlining the algorithm to run the supplemenatary analysis.", 
            paths_to_rmd_dir_1L_chr = NA_character_, pkg_dirs_chr = "Markdown", 
            packages_chr = "TTU", nms_of_rmd_chr = "Supplement.Rmd") %>% 
            dplyr::filter(rprt_nms_chr != "AAA_PMRY_ANLYS_MTH")
    }
    if (add_sharing_rprt_1L_lgl) {
        rprt_lup <- rprt_lup %>% tibble::add_case(rprt_nms_chr = "AAA_SHARING_MTH", 
            title_chr = "Supplementary report outlining the algorithm to create and disseminate shareable study output.", 
            paths_to_rmd_dir_1L_chr = NA_character_, pkg_dirs_chr = "Markdown", 
            packages_chr = "TTU", nms_of_rmd_chr = "Share.Rmd")
    }
    if (!is.null(start_at_int[1])) {
        rprt_lup <- dplyr::mutate(rprt_lup, title_chr = dplyr::case_when(rprt_nms_chr %in% 
            c("AAA_PMRY_ANLYS_MTH") ~ paste0("Methods Report ", 
            start_at_int[1], ": Analysis Program (", "Primary Analysis", 
            ")"), rprt_nms_chr %in% c("AAA_SUPLRY_ANLYS_MTH") ~ 
            paste0("Methods Report ", start_at_int[1] + 3, ": Analysis Program (", 
                "Secondary Analysis", ")"), rprt_nms_chr %in% 
            c("AAA_RPRT_WRTNG_MTH") ~ paste0("Methods Report ", 
            start_at_int[1] + 1, ": Reporting Program"), rprt_nms_chr %in% 
            c("AAA_SHARING_MTH") ~ paste0("Methods Report ", 
            start_at_int[1] + 2, ": Sharing Program"), rprt_nms_chr %in% 
            c("AAA_TTU_MDL_CTG") ~ paste0("Results Report ", 
            ifelse(is.null(reference_1L_int), start_at_int[2], 
                start_at_int[2] + reference_1L_int), ": Catalogue of longitudinal models (", 
            ifelse(is.null(reference_1L_int), "Primary Analysis", 
                paste0("Secondary Analysis ", LETTERS[reference_1L_int])), 
            ")"), T ~ title_chr))
    }
    if (!is.null(reference_1L_int)) {
        rprt_lup <- dplyr::mutate(rprt_lup, rprt_nms_chr = dplyr::case_when(rprt_nms_chr %in% 
            c("AAA_TTU_MDL_CTG") ~ paste0("AAA_TTU_MDL_CTG", 
            ifelse(is.null(reference_1L_int), "", ifelse(reference_1L_int == 
                0, "", paste0("-", reference_1L_int)))), T ~ 
            rprt_nms_chr))
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
#' @param scaling_fctr_dbl Scaling factor (a double vector), Default: 1
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param ungroup_1L_lgl Ungroup (a logical vector of length one), Default: F
#' @return Transformed for model input (a tibble)
#' @rdname transform_tb_to_mdl_inp
#' @export 
#' @importFrom ready4use remove_labels_from_ds
#' @importFrom dplyr select all_of group_by arrange mutate across first lag ungroup
#' @importFrom rlang sym
#' @importFrom purrr reduce
#' @importFrom stats na.omit
transform_tb_to_mdl_inp <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr, 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline", drop_all_msng_1L_lgl = T, 
    scaling_fctr_dbl = 1, tfmn_1L_chr = "NTF", ungroup_1L_lgl = F) 
{
    if (length(scaling_fctr_dbl) != length(predr_vars_nms_chr)) {
        scaling_fctr_dbl <- rep(scaling_fctr_dbl[1], length(predr_vars_nms_chr))
    }
    data_tb <- data.frame(data_tb) %>% ready4use::remove_labels_from_ds()
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
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% add_tfd_var_to_ds(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr, depnt_var_max_val_1L_dbl = 0.999)
    if (drop_all_msng_1L_lgl) {
        tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% stats::na.omit()
    }
    if (ungroup_1L_lgl) {
        tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::ungroup()
    }
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% transform_uid_var(id_var_nm_1L_chr = id_var_nm_1L_chr)
    return(tfd_for_mdl_inp_tb)
}
#' Transform table to round variables
#' @description transform_tbl_to_rnd_vars() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform table to round variables. Function argument ds_tb specifies the object to be updated. Argument nbr_of_digits_1L_int provides the object to be updated. The function returns Transformed dataset (a tibble).
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
#' Transform timepoint values
#' @description transform_timepoint_vals() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform timepoint values. Function argument timepoint_vals_chr specifies the object to be updated. Argument timepoint_levels_chr provides the object to be updated. The function returns Timepoint values (a character vector).
#' @param timepoint_vals_chr Timepoint values (a character vector)
#' @param timepoint_levels_chr Timepoint levels (a character vector)
#' @param bl_val_1L_chr Baseline value (a character vector of length one)
#' @return Timepoint values (a character vector)
#' @rdname transform_timepoint_vals
#' @export 
#' @keywords internal
transform_timepoint_vals <- function (timepoint_vals_chr, timepoint_levels_chr, bl_val_1L_chr) 
{
    if (length(timepoint_vals_chr) == 1) {
        timepoint_vals_chr <- bl_val_1L_chr
    }
    else {
        unique_vals_chr <- unique(timepoint_vals_chr)
        if (length(timepoint_vals_chr) > length(unique_vals_chr)) 
            timepoint_vals_chr <- c(unique_vals_chr, setdiff(c(bl_val_1L_chr, 
                setdiff(timepoint_levels_chr, bl_val_1L_chr)), 
                unique_vals_chr)[1:(length(timepoint_vals_chr) - 
                length(unique_vals_chr))])
    }
    return(timepoint_vals_chr)
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
#' Transform unique identifier variable
#' @description transform_uid_var() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform unique identifier variable. Function argument data_tb specifies the object to be updated. Argument id_var_nm_1L_chr provides the object to be updated. The function returns Transformed data (a tibble).
#' @param data_tb Data (a tibble)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one)
#' @param rename_tb Rename (a tibble), Default: NULL
#' @param old_new_chr Old new (a character vector), Default: c("old_id_xx", "new_id_int")
#' @return Transformed data (a tibble)
#' @rdname transform_uid_var
#' @export 
#' @importFrom dplyr pull mutate
#' @importFrom purrr flatten_chr flatten_int flatten_dbl map
#' @importFrom rlang sym
#' @importFrom ready4 get_from_lup_obj
#' @keywords internal
transform_uid_var <- function (data_tb, id_var_nm_1L_chr, rename_tb = NULL, old_new_chr = c("old_id_xx", 
    "new_id_int")) 
{
    if (is.null(rename_tb)) {
        rename_tb <- make_uid_rename_lup(data_tb, id_var_nm_1L_chr = id_var_nm_1L_chr)
    }
    if (!identical(rename_tb$old_id_xx, rename_tb$new_id_int)) {
        fn <- ifelse("character" %in% class(rename_tb %>% dplyr::pull(old_new_chr[2])), 
            purrr::flatten_chr, ifelse("integer" %in% class(rename_tb %>% 
                dplyr::pull(old_new_chr[2])), purrr::flatten_int, 
                purrr::flatten_dbl))
        tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(id_var_nm_1L_chr), 
            !!rlang::sym(id_var_nm_1L_chr) %>% purrr::map(~ready4::get_from_lup_obj(rename_tb, 
                match_value_xx = .x, match_var_nm_1L_chr = old_new_chr[1], 
                target_var_nm_1L_chr = old_new_chr[2], evaluate_1L_lgl = F)) %>% 
                fn()))
    }
    else {
        tfd_data_tb <- data_tb
    }
    return(tfd_data_tb)
}
