#' Make fake ts data
#' @description make_fake_ts_data() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make fake ts data. The function returns Fk data (a tibble).
#' @param outp_smry_ls Output smry (a list)
#' @return Fk data (a tibble)
#' @rdname make_fake_ts_data
#' @export 
#' @importFrom synthpop syn
#' @importFrom purrr map_lgl
#' @importFrom dplyr mutate across all_of
make_fake_ts_data <- function (outp_smry_ls) 
{
    data_tb <- outp_smry_ls$scored_data_tb %>% transform_tb_to_mdl_inp(dep_var_nm_1L_chr = outp_smry_ls$dep_var_nm_1L_chr, 
        predr_vars_nms_chr = outp_smry_ls$predr_cmprsns_tb$predr_chr, 
        id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
        round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr)
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = names(data_tb)[names(data_tb) != 
        outp_smry_ls$id_var_nm_1L_chr], seed = outp_smry_ls$seed_1L_int)
    dep_vars_chr <- names(fk_data_ls$syn)[names(fk_data_ls$syn) %>% 
        purrr::map_lgl(~startsWith(.x, outp_smry_ls$dep_var_nm_1L_chr))]
    fk_data_tb <- fk_data_ls$syn %>% dplyr::mutate(dplyr::across(dplyr::all_of(dep_vars_chr), 
        ~NA_real_))
    return(fk_data_tb)
}
#' Make model
#' @description make_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make model. The function is called for its side effects and does not return a value.
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predr var name (a character vector of length one)
#' @param covar_var_nms_chr Covar var names (a character vector), Default: 'NA'
#' @param mdl_type_1L_chr Model type (a character vector of length one), Default: 'OLS_NTF'
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param control_1L_chr Control (a character vector of length one), Default: 'NA'
#' @param start_1L_chr Start (a character vector of length one), Default: NULL
#' @return Model (a model)
#' @rdname make_mdl
#' @export 
#' @importFrom utils data
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom stringi stri_locate_last_fixed
#' @importFrom stringr str_sub
make_mdl <- function (data_tb, dep_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF", 
    predr_var_nm_1L_chr, covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF", 
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- transform_ds_for_mdlng(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    if (is.null(start_1L_chr)) {
        start_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup, 
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
            target_var_nm_1L_chr = "start_chr", evaluate_lgl = F)
    }
    if (!is.na(control_1L_chr)) {
        idx_1L_int <- 1 + stringi::stri_locate_last_fixed(mdl_type_1L_chr, 
            "_")[1, 1] %>% as.vector()
        link_1L_chr <- stringr::str_sub(mdl_type_1L_chr, start = idx_1L_int)
        link_1L_chr <- ifelse(link_1L_chr == "LOG", "log", ifelse(link_1L_chr == 
            "LGT", "logit", ifelse(link_1L_chr == "CLL", "cloglog", 
            "ERROR")))
    }
    mdl_1L_chr <- paste0(ready4fun::get_from_lup_obj(mdl_types_lup, 
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
        target_var_nm_1L_chr = "fn_chr", evaluate_lgl = F), "(", 
        transform_dep_var_nm(dep_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr), 
        " ~ ", predr_var_nm_1L_chr, ifelse(is.na(covar_var_nms_chr[1]), 
            "", paste0(" + ", paste0(covar_var_nms_chr, collapse = " + "))), 
        ", data = data_tb", ifelse(!is.na(ready4fun::get_from_lup_obj(mdl_types_lup, 
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
            target_var_nm_1L_chr = "family_chr", evaluate_lgl = F)), 
            paste0(", family = ", ready4fun::get_from_lup_obj(mdl_types_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
                target_var_nm_1L_chr = "family_chr", evaluate_lgl = F)), 
            ""), ifelse(!is.na(start_1L_chr), ", ", ""), ifelse(!is.na(control_1L_chr), 
            paste0("link=\"", link_1L_chr, "\",control=", control_1L_chr, 
                "("), ""), ifelse(!is.na(start_1L_chr), paste0("start=c(", 
            start_1L_chr, ")"), ""), ifelse(!is.na(control_1L_chr), 
            ")", ""), ")")
    model_mdl <- eval(parse(text = mdl_1L_chr))
    return(model_mdl)
}
#' Make model names
#' @description make_mdl_nms_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make model names list. The function returns Model names (a list).
#' @param predr_vars_nms_ls Predr vars names (a list)
#' @param mdl_types_chr Model types (a character vector)
#' @return Model names (a list)
#' @rdname make_mdl_nms_ls
#' @export 
#' @importFrom purrr map2
make_mdl_nms_ls <- function (predr_vars_nms_ls, mdl_types_chr) 
{
    mdl_nms_ls <- purrr::map2(predr_vars_nms_ls, make_unique_ls_elmt_idx_int(predr_vars_nms_ls), 
        ~paste0(.x[1], "_", ifelse(is.na(.x[2]), "", paste0(.x[2], 
            "_")), .y, "_", mdl_types_chr))
    return(mdl_nms_ls)
}
#' Make predn dataset with one predr
#' @description make_predn_ds_with_one_predr() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make predn dataset with one predr. The function returns Predn dataset (a tibble).
#' @param model_mdl PARAM_DESCRIPTION
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predr var name (a character vector of length one)
#' @param predr_vals_dbl Predr values (a double vector)
#' @param pred_type_1L_chr Pred type (a character vector of length one), Default: NULL
#' @return Predn dataset (a tibble)
#' @rdname make_predn_ds_with_one_predr
#' @export 
#' @importFrom tibble tibble
#' @importFrom rlang sym
#' @importFrom dplyr mutate
#' @importFrom stats predict
make_predn_ds_with_one_predr <- function (model_mdl, dep_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF", 
    predr_var_nm_1L_chr, predr_vals_dbl, pred_type_1L_chr = NULL) 
{
    predn_ds_tb <- tibble::tibble(`:=`(!!rlang::sym(predr_var_nm_1L_chr), 
        predr_vals_dbl))
    predn_ds_tb <- predn_ds_tb %>% dplyr::mutate(`:=`(!!rlang::sym(dep_var_nm_1L_chr), 
        stats::predict(model_mdl, newdata = predn_ds_tb, type = pred_type_1L_chr) %>% 
            calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, 
                tfmn_is_outp_1L_lgl = T)))
    return(predn_ds_tb)
}
#' Make shareable model
#' @description make_shareable_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make shareable model. The function is called for its side effects and does not return a value.
#' @param data_tb Data (a tibble)
#' @param mdl_smry_tb Model smry (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'CLL'
#' @param mdl_type_1L_chr Model type (a character vector of length one), Default: 'OLS_CLL'
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param control_1L_chr Control (a character vector of length one), Default: 'NA'
#' @param start_1L_chr Start (a character vector of length one), Default: 'NA'
#' @param seed_1L_int Seed (an integer vector of length one), Default: 12345
#' @return Model (a model)
#' @rdname make_shareable_mdl
#' @export 
#' @importFrom utils data
#' @importFrom synthpop syn
#' @importFrom dplyr mutate case_when filter slice
#' @importFrom purrr map_chr
#' @importFrom stringr str_replace_all
#' @importFrom assertthat assert_that
make_shareable_mdl <- function (data_tb, mdl_smry_tb, dep_var_nm_1L_chr = "utl_total_w", 
    id_var_nm_1L_chr = "fkClientID", tfmn_1L_chr = "CLL", mdl_type_1L_chr = "OLS_CLL", 
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NA_character_, 
    seed_1L_int = 12345L) 
{
    if (is.null(mdl_types_lup)) 
        utils::data(mdl_types_lup, envir = environment())
    all_var_nms_chr <- names(data_tb)
    tfd_dep_var_nm_1L_chr <- all_var_nms_chr[all_var_nms_chr %>% 
        startsWith(dep_var_nm_1L_chr)]
    predr_var_nms_chr <- setdiff(all_var_nms_chr, c(tfd_dep_var_nm_1L_chr, 
        id_var_nm_1L_chr))
    if (length(predr_var_nms_chr) > 1) {
        covar_var_nms_chr <- predr_var_nms_chr[2:length(predr_var_nms_chr)]
    }
    else {
        covar_var_nms_chr <- NA_character_
    }
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = all_var_nms_chr[all_var_nms_chr != 
        id_var_nm_1L_chr], seed = seed_1L_int)
    model_mdl <- make_mdl(fk_data_ls$syn, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nms_chr[1], covar_var_nms_chr = covar_var_nms_chr, 
        tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr, 
        mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr, 
        start_1L_chr = start_1L_chr)
    par_nms_chr <- model_mdl$coefficients %>% names()
    mdl_smry_tb <- mdl_smry_tb %>% dplyr::mutate(Parameter = dplyr::case_when(Parameter == 
        "Intercept" ~ "(Intercept)", TRUE ~ purrr::map_chr(Parameter, 
        ~stringr::str_replace_all(.x, " ", "_")))) %>% dplyr::filter(Parameter %in% 
        par_nms_chr) %>% dplyr::slice(match(par_nms_chr, Parameter))
    assertthat::assert_that(all(par_nms_chr == mdl_smry_tb$Parameter), 
        msg = "Parameter names mismatch between data and model summary table")
    model_mdl$coefficients <- mdl_smry_tb$Estimate
    names(model_mdl$coefficients) <- par_nms_chr
    return(model_mdl)
}
#' Make smry of brm model
#' @description make_smry_of_brm_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make smry of brm model. The function returns Smry of brm model (a tibble).
#' @param mdl_ls Model (a list)
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param fn Function (a function), Default: calculate_rmse
#' @param mdl_nm_1L_chr Model name (a character vector of length one), Default: 'NA'
#' @param seed_1L_dbl Seed (a double vector of length one), Default: 23456
#' @return Smry of brm model (a tibble)
#' @rdname make_smry_of_brm_mdl
#' @export 
#' @importFrom stats predict
#' @importFrom brms bayes_R2
#' @importFrom psych describe
#' @importFrom dplyr pull mutate rename select
#' @importFrom rlang sym
#' @importFrom purrr map flatten_chr
make_smry_of_brm_mdl <- function (mdl_ls, data_tb, dep_var_nm_1L_chr = "utl_total_w", 
    predr_vars_nms_chr, fn = calculate_rmse, mdl_nm_1L_chr = NA_character_, 
    seed_1L_dbl = 23456) 
{
    if (is.na(mdl_nm_1L_chr)) 
        mdl_nm_1L_chr <- predr_vars_nms_chr[1]
    set.seed(seed_1L_dbl)
    predictions <- stats::predict(mdl_ls, summary = F)
    coef <- summary(mdl_ls, digits = 4)$fixed
    coef <- coef[1:nrow(coef), 1:4]
    R2 <- brms::bayes_R2(mdl_ls)
    RMSE <- psych::describe(apply(predictions, 1, fn, y_dbl = data_tb %>% 
        dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))), quant = c(0.25, 
        0.75), skew = F, ranges = F)
    RMSE <- cbind(RMSE$mean, RMSE$sd, RMSE$Q0.25, RMSE$Q0.75) %>% 
        as.vector()
    Sigma <- summary(mdl_ls, digits = 4)$spec_par[1:4]
    smry_of_brm_mdl_tb <- data.frame(round(rbind(coef, R2, RMSE, 
        Sigma), 3)) %>% dplyr::mutate(Parameter = c("Intercept", 
        purrr::map(predr_vars_nms_chr, ~paste0(.x, c(" baseline", 
            " change"))) %>% purrr::flatten_chr(), "R2", "RMSE", 
        "Sigma"), Model = mdl_nm_1L_chr) %>% dplyr::mutate(`95% CI` = paste(l.95..CI, 
        ",", u.95..CI)) %>% dplyr::rename(SE = Est.Error) %>% 
        dplyr::select(Model, Parameter, Estimate, SE, `95% CI`)
    return(smry_of_brm_mdl_tb)
}
#' Make smry of model
#' @description make_smry_of_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make smry of model. The function returns Smry of one predr model (a tibble).
#' @param data_tb Data (a tibble)
#' @param model_mdl PARAM_DESCRIPTION
#' @param n_folds_1L_int N folds (an integer vector of length one), Default: 10
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param start_1L_chr Start (a character vector of length one), Default: NULL
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predr var name (a character vector of length one)
#' @param covar_var_nms_chr Covar var names (a character vector), Default: 'NA'
#' @param mdl_type_1L_chr Model type (a character vector of length one), Default: 'OLS_NTF'
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param pred_type_1L_chr Pred type (a character vector of length one), Default: NULL
#' @return Smry of one predr model (a tibble)
#' @rdname make_smry_of_mdl
#' @export 
#' @importFrom utils data
#' @importFrom dplyr filter pull summarise_all mutate select everything
#' @importFrom rlang sym
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom purrr map_dfr
#' @importFrom stats predict
#' @importFrom tibble tibble
#' @importFrom caret R2 RMSE MAE
make_smry_of_mdl <- function (data_tb, model_mdl, n_folds_1L_int = 10, dep_var_nm_1L_chr = "utl_total_w", 
    start_1L_chr = NULL, tfmn_1L_chr = "NTF", predr_var_nm_1L_chr, 
    covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF", 
    mdl_types_lup = NULL, pred_type_1L_chr = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- data_tb %>% dplyr::filter(!is.na(!!rlang::sym(predr_var_nm_1L_chr)))
    data_tb <- transform_ds_for_mdlng(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    mdl_desc_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup, 
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
        target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F)
    folds_ls <- make_folds_ls(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        n_folds_1L_int = n_folds_1L_int)
    smry_of_one_predr_mdl_tb <- purrr::map_dfr(folds_ls, ~{
        model_mdl <- make_mdl(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
            start_1L_chr = start_1L_chr, tfmn_1L_chr = tfmn_1L_chr, 
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr, 
            mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup)
        pred_old_dbl <- stats::predict(model_mdl, type = pred_type_1L_chr)
        pred_new_dbl <- stats::predict(model_mdl, newdata = data_tb[.x, 
            ], type = pred_type_1L_chr) %>% calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, 
            tfmn_is_outp_1L_lgl = T)
        tibble::tibble(Rsquared = caret::R2(pred_old_dbl, data_tb[-.x, 
            ] %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr)), 
            form = "traditional"), RMSE = caret::RMSE(pred_old_dbl, 
            data_tb[-.x, ] %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))), 
            MAE = caret::MAE(pred_old_dbl, data_tb[-.x, ] %>% 
                dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))), 
            RsquaredP = caret::R2(pred_new_dbl, data_tb[.x, ] %>% 
                dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr)), 
                form = "traditional"), RMSEP = caret::RMSE(pred_new_dbl, 
                data_tb[.x, ] %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))), 
            MAEP = caret::MAE(pred_new_dbl, data_tb[.x, ] %>% 
                dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))))
    }) %>% dplyr::summarise_all(mean) %>% dplyr::mutate(Model = mdl_desc_1L_chr) %>% 
        dplyr::select(Model, dplyr::everything())
    return(smry_of_one_predr_mdl_tb)
}
#' Make smry of ts model
#' @description make_smry_of_ts_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make smry of ts model. The function returns Smry of ts model (a list).
#' @param data_tb Data (a tibble)
#' @param fn Function (a function)
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one), Default: 'NA'
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iters (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Smry of ts model (a list)
#' @rdname make_smry_of_ts_mdl
#' @export 
#' @importFrom rlang exec
make_smry_of_ts_mdl <- function (data_tb, fn, predr_vars_nms_chr, mdl_nm_1L_chr, path_to_write_to_1L_chr = NA_character_, 
    dep_var_nm_1L_chr = "utl_total_w", id_var_nm_1L_chr = "fkClientID", 
    round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", 
    backend_1L_chr = getOption("brms.backend", "rstan"), iters_1L_int = 4000L, 
    seed_1L_int = 1000L) 
{
    tfd_data_tb <- transform_tb_to_mdl_inp(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr)
    tfd_dep_var_nm_1L_chr <- ifelse(identical(fn, fit_clg_log_tfmn), 
        transform_dep_var_nm_for_cll(dep_var_nm_1L_chr), dep_var_nm_1L_chr)
    args_ls <- list(data_tb = tfd_data_tb, dep_var_nm_1L_chr = tfd_dep_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, iters_1L_int = iters_1L_int, 
        backend_1L_chr = backend_1L_chr, seed_1L_int = seed_1L_int)
    mdl_ls <- rlang::exec(fn, !!!args_ls)
    smry_of_ts_mdl_ls <- list(smry_of_ts_mdl_tb = make_smry_of_brm_mdl(mdl_ls, 
        data_tb = tfd_data_tb, dep_var_nm_1L_chr = tfd_dep_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, fn = ifelse(identical(fn, 
            fit_gsn_log_lnk), calculate_rmse, calculate_rmse_tfmn), 
        mdl_nm_1L_chr = mdl_nm_1L_chr))
    if (!is.na(path_to_write_to_1L_chr)) {
        smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr <- paste0(path_to_write_to_1L_chr, 
            "/", mdl_nm_1L_chr, ".RDS")
        if (file.exists(smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)) 
            file.remove(smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)
        saveRDS(mdl_ls, smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)
        smry_of_ts_mdl_ls$paths_to_mdl_plts_chr <- write_brm_model_plts(mdl_ls, 
            tfd_data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
            mdl_nm_1L_chr = mdl_nm_1L_chr, path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
            round_var_nm_1L_chr = round_var_nm_1L_chr, tfmn_fn = ifelse(identical(fn, 
                fit_gsn_log_lnk), function(x) {
                x
            }, function(x) {
                1 - exp(-exp(x))
            }))
    }
    return(smry_of_ts_mdl_ls)
}
