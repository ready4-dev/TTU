#' Transform data tibble for comparison
#' @description transform_data_tb_for_cmprsn() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform data tibble for comparison. Function argument data_tb specifies the object to be updated. Argument model_mdl provides the object to be updated. The function returns Transformed data (a tibble).
#' @param data_tb Data (a tibble)
#' @param model_mdl Model (a model)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param source_data_nm_1L_chr Source data name (a character vector of length one), Default: 'Original'
#' @param tf_type_1L_chr Transform type (a character vector of length one), Default: 'Predicted'
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @param tfmn_for_bnml_1L_lgl Transformation for binomial (a logical vector of length one), Default: F
#' @param family_1L_chr Family (a character vector of length one), Default: 'NA'
#' @return Transformed data (a tibble)
#' @rdname transform_data_tb_for_cmprsn
#' @export 
#' @importFrom stats predict simulate rnorm sigma
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @keywords internal
transform_data_tb_for_cmprsn <- function (data_tb, model_mdl, depnt_var_nm_1L_chr = "utl_total_w", 
    source_data_nm_1L_chr = "Original", tf_type_1L_chr = "Predicted", 
    predn_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_) 
{
    if (tf_type_1L_chr == "Predicted") 
        new_data_dbl <- stats::predict(model_mdl, type = predn_type_1L_chr)
    if (tf_type_1L_chr == "Simulated" & !tfmn_for_bnml_1L_lgl) 
        new_data_dbl <- stats::simulate(model_mdl)$sim_1
    if (tf_type_1L_chr == "Simulated" & tfmn_for_bnml_1L_lgl) 
        new_data_dbl <- (stats::predict(model_mdl) + stats::rnorm(nrow(data_tb), 
            0, stats::sigma(model_mdl))) %>% calculate_dpnt_var_tfmn(tfmn_1L_chr = ifelse(family_1L_chr == 
            "quasibinomial(log)", "LOG", ifelse(family_1L_chr == 
            "quasibinomial(logit)", "LOGIT", ifelse(family_1L_chr == 
            "quasibinomial(cloglog)", "CLL", "NTF"))), tfmn_is_outp_1L_lgl = T)
    tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(tf_type_1L_chr), 
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
#' Transform dependent variable name for complementary log log
#' @description transform_depnt_var_nm_for_cll() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dependent variable name for complementary log log. Function argument depnt_var_nm_1L_chr specifies the object to be updated. The function returns Transformed dependent variable name (a character vector of length one).
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one)
#' @return Transformed dependent variable name (a character vector of length one)
#' @rdname transform_depnt_var_nm_for_cll
#' @export 

#' @keywords internal
transform_depnt_var_nm_for_cll <- function (depnt_var_nm_1L_chr) 
{
    tfd_depnt_var_nm_1L_chr <- paste0(depnt_var_nm_1L_chr, "_cloglog")
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
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param class_fn_1L_chr Class function (a character vector of length one), Default: 'youthvars::youthvars_aqol6d_adol'
#' @return Transformed dataset (a tibble)
#' @rdname transform_mdl_vars_with_clss
#' @export 
#' @importFrom tibble add_case
#' @importFrom purrr reduce
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom dplyr mutate
#' @importFrom rlang sym exec
#' @keywords internal
transform_mdl_vars_with_clss <- function (ds_tb, predictors_lup = NULL, prototype_lup = NULL, 
    depnt_var_nm_1L_chr = "aqol6d_total_w", class_fn_1L_chr = "youthvars::youthvars_aqol6d_adol") 
{
    if (is.null(predictors_lup)) 
        data("predictors_lup", package = "youthvars", envir = environment())
    if (is.null(prototype_lup)) 
        data("prototype_lup", package = "TTU", envir = environment())
    predictors_lup <- tibble::add_case(predictors_lup, short_name_chr = depnt_var_nm_1L_chr, 
        class_chr = "numeric", class_fn_chr = class_fn_1L_chr)
    tfd_ds_tb <- purrr::reduce(predictors_lup$short_name_chr, 
        .init = ds_tb, ~if (.y %in% names(.x)) {
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
            .x %>% dplyr::mutate(`:=`(!!rlang::sym(.y), rlang::exec(ready4fun::get_from_lup_obj(predictors_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .y, 
                target_var_nm_1L_chr = "class_fn_chr", evaluate_lgl = T), 
                !!rlang::sym(.y) %>% fn)))
        }
        else {
            .x
        })
    return(tfd_ds_tb)
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
#' @param ungroup_1L_lgl Ungroup (a logical vector of length one), Default: F
#' @param add_cll_tfmn_1L_lgl Add complementary log log transformation (a logical vector of length one), Default: T
#' @return Transformed for model input (a tibble)
#' @rdname transform_tb_to_mdl_inp
#' @export 
#' @importFrom dplyr select all_of group_by arrange mutate across first lag ungroup
#' @importFrom rlang sym
#' @importFrom purrr reduce
#' @importFrom stats na.omit
#' @keywords internal
transform_tb_to_mdl_inp <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr, 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline", drop_all_msng_1L_lgl = T, 
    scaling_fctr_dbl = 0.01, ungroup_1L_lgl = F, add_cll_tfmn_1L_lgl = T) 
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
                    "Baseline", 0, (. - dplyr::lag(.)) * scaling_fctr_dbl[idx_1L_int]))))
        })
    if (add_cll_tfmn_1L_lgl) {
        tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_depnt_var_nm_for_cll(depnt_var_nm_1L_chr)), 
            log(-log(1 - ifelse(!!rlang::sym(depnt_var_nm_1L_chr) == 
                1, 0.999, !!rlang::sym(depnt_var_nm_1L_chr))))))
    }
    if (drop_all_msng_1L_lgl) {
        tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% stats::na.omit()
    }
    if (ungroup_1L_lgl) {
        tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::ungroup()
    }
    return(tfd_for_mdl_inp_tb)
}
#' Transform time series model data
#' @description transform_ts_mdl_data() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform time series model data. Function argument mdl_ls specifies the object to be updated. Argument data_tb provides the object to be updated. The function returns Cnfdl (a list of models).
#' @param mdl_ls Model list (a list of models)
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @return Cnfdl (a list of models)
#' @rdname transform_ts_mdl_data
#' @export 
#' @importFrom dplyr select all_of summarise across everything
#' @importFrom purrr map flatten_chr
#' @keywords internal
transform_ts_mdl_data <- function (mdl_ls, data_tb, depnt_var_nm_1L_chr = "aqol6d_total_w", 
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
