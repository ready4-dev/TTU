#' Write all algorithm outputs
#' @description write_all_alg_outps() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write all algorithm outputs. The function returns Output summary (a list).
#' @param scored_data_tb Scored data (a tibble)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param candidate_predrs_chr Candidate predictors (a character vector)
#' @param candidate_covar_nms_chr Candidate covariate names (a character vector)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one), Default: 'Baseline'
#' @param mdl_types_chr Model types (a character vector), Default: 'NA'
#' @param prefd_mdl_types_chr Preferred model types (a character vector), Default: 'NA'
#' @param choose_from_pfx_chr Choose from prefix (a character vector), Default: c("GLM", "OLS", "BET")
#' @param prefd_covars_chr Preferred covariates (a character vector), Default: 'NA'
#' @param seed_1L_int Seed (an integer vector of length one), Default: 12345
#' @param folds_1L_int Folds (an integer vector of length one), Default: 10
#' @param max_nbr_of_boruta_mdl_runs_int Maximum number of boruta model runs (an integer vector), Default: 300
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param ds_smry_ls Dataset summary (a list)
#' @return Output summary (a list)
#' @rdname write_all_alg_outps
#' @export 
#' @importFrom utils data
#' @importFrom youthvars transform_ds_for_tstng
#' @importFrom ready4fun get_from_lup_obj
write_all_alg_outps <- function (scored_data_tb, path_to_write_to_1L_chr, depnt_var_nm_1L_chr = "utl_total_w", 
    candidate_predrs_chr, candidate_covar_nms_chr, id_var_nm_1L_chr = "fkClientID", 
    round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", 
    mdl_types_chr = NA_character_, prefd_mdl_types_chr = NA_character_, 
    choose_from_pfx_chr = c("GLM", "OLS", "BET"), prefd_covars_chr = NA_character_, 
    seed_1L_int = 12345, folds_1L_int = 10L, max_nbr_of_boruta_mdl_runs_int = 300L, 
    mdl_types_lup = NULL, ds_smry_ls) 
{
    set.seed(seed_1L_int)
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    if (is.na(mdl_types_chr[1])) {
        mdl_types_chr <- mdl_types_lup$short_name_chr
    }
    bl_tb <- youthvars::transform_ds_for_tstng(scored_data_tb, 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, candidate_predrs_chr = candidate_predrs_chr, 
        dep_var_max_val_1L_dbl = 0.999, round_var_nm_1L_chr = round_var_nm_1L_chr, 
        round_val_1L_chr = round_bl_val_1L_chr)
    candidate_predrs_chr <- reorder_cndt_predrs_chr(candidate_predrs_chr, 
        data_tb = bl_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr)
    predr_var_nm_1L_chr <- candidate_predrs_chr[1]
    predr_var_desc_1L_chr <- candidate_predrs_lup %>% ready4fun::get_from_lup_obj(match_value_xx = predr_var_nm_1L_chr, 
        match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = "long_name_chr", 
        evaluate_lgl = F)
    predr_vals_dbl <- make_predr_vals(predr_var_nm_1L_chr, candidate_predrs_lup = candidate_predrs_lup)
    smry_of_sngl_predr_mdls_tb <- write_sngl_predr_multi_mdls_outps(data_tb = bl_tb, 
        folds_1L_int = folds_1L_int, mdl_types_chr = mdl_types_chr, 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, predr_var_nm_1L_chr = predr_var_nm_1L_chr, 
        predr_var_desc_1L_chr = predr_var_nm_1L_chr, predr_vals_dbl = predr_vals_dbl, 
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, new_dir_nm_1L_chr = "A_Candidate_Mdls_Cmprsn", 
        mdl_types_lup = mdl_types_lup, dictionary_tb = ds_smry_ls$dictionary_tb)
    if (is.na(prefd_mdl_types_chr[1])) {
        prefd_mdl_types_chr <- make_prefd_mdls_vec(smry_of_sngl_predr_mdls_tb, 
            choose_from_pfx_chr = choose_from_pfx_chr)
    }
    predr_cmprsn_tb <- write_predr_cmprsn_outps(data_tb = bl_tb, 
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, new_dir_nm_1L_chr = "B_Candidate_Predrs_Cmprsn", 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, candidate_predrs_chr = candidate_predrs_chr, 
        max_nbr_of_boruta_mdl_runs_int = max_nbr_of_boruta_mdl_runs_int)
    smry_of_mdl_sngl_predrs_tb <- write_mdl_type_multi_outps(data_tb = bl_tb, 
        folds_1L_int = folds_1L_int, predrs_var_nms_chr = predr_cmprsn_tb$predr_chr, 
        mdl_type_1L_chr = prefd_mdl_types_chr[1], depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, new_dir_nm_1L_chr = "C_Predrs_Sngl_Mdl_Cmprsn", 
        fl_nm_pfx_1L_chr = "C_PREDR", mdl_types_lup = mdl_types_lup)
    bl_tb <- scored_data_tb %>% youthvars::transform_ds_for_tstng(candidate_predrs_chr = candidate_predrs_chr, 
        covar_var_nms_chr = candidate_covar_nms_chr, remove_all_msng_1L_lgl = T, 
        round_val_1L_chr = round_bl_val_1L_chr)
    mdls_with_covars_smry_tb <- write_mdl_type_covars_mdls(bl_tb, 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, predrs_var_nms_chr = candidate_predrs_chr, 
        covar_var_nms_chr = candidate_covar_nms_chr, mdl_type_1L_chr = prefd_mdl_types_chr[1], 
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, new_dir_nm_1L_chr = "D_Predr_Covars_Cmprsn", 
        fl_nm_pfx_1L_chr = "D_CT", mdl_types_lup = mdl_types_lup)
    signt_covars_chr <- get_signft_covars(mdls_with_covars_smry_tb = mdls_with_covars_smry_tb, 
        covar_var_nms_chr = candidate_covar_nms_chr)
    if (is.na(prefd_covars_chr)) 
        prefd_covars_chr <- signt_covars_chr
    dud_tb <- write_mdl_type_multi_outps(data_tb = bl_tb, folds_1L_int = NULL, 
        start_1L_chr = NA_character_, predrs_var_nms_chr = predr_cmprsn_tb$predr_chr, 
        covar_var_nms_chr = candidate_covar_nms_chr, mdl_type_1L_chr = prefd_mdl_types_chr[2], 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = "E_Predrs_W_Covars_Sngl_Mdl_Cmprsn", 
        mdl_types_lup = mdl_types_lup, fl_nm_pfx_1L_chr = "E_CK_CV")
    predr_vars_nms_ls <- make_predr_vars_nms_ls(main_predrs_chr = predr_cmprsn_tb$predr_chr, 
        covars_ls = list(prefd_covars_chr))
    mdl_nms_ls <- make_mdl_nms_ls(predr_vars_nms_ls, mdl_types_chr = prefd_mdl_types_chr)
    outp_smry_ls <- list(scored_data_tb = scored_data_tb, smry_of_sngl_predr_mdls_tb = smry_of_sngl_predr_mdls_tb, 
        prefd_mdl_types_chr = prefd_mdl_types_chr, predr_cmprsn_tb = predr_cmprsn_tb, 
        smry_of_mdl_sngl_predrs_tb = smry_of_mdl_sngl_predrs_tb, 
        mdls_with_covars_smry_tb = mdls_with_covars_smry_tb, 
        signt_covars_chr = signt_covars_chr, prefd_covars_chr = prefd_covars_chr, 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, predr_vars_nms_ls = predr_vars_nms_ls, 
        mdl_nms_ls = mdl_nms_ls, id_var_nm_1L_chr = id_var_nm_1L_chr, 
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr, 
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, seed_1L_int = seed_1L_int, 
        folds_1L_int = folds_1L_int, max_nbr_of_boruta_mdl_runs_int = max_nbr_of_boruta_mdl_runs_int, 
        mdl_types_lup = mdl_types_lup, file_paths_chr = list.files(path_to_write_to_1L_chr, 
            recursive = T))
    return(outp_smry_ls)
}
#' Write box cox transformation
#' @description write_box_cox_tfmn() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write box cox transformation. The function returns Path to plot (a character vector of length one).
#' @param data_tb Data (a tibble)
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param covar_var_nms_chr Covariate variable names (a character vector), Default: 'NA'
#' @param fl_nm_pfx_1L_chr File name prefix (a character vector of length one), Default: 'A_RT'
#' @param height_1L_dbl Height (a double vector of length one), Default: 6
#' @param width_1L_dbl Width (a double vector of length one), Default: 6
#' @param start_1L_chr Start (a character vector of length one), Default: NULL
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @return Path to plot (a character vector of length one)
#' @rdname write_box_cox_tfmn
#' @export 
#' @importFrom utils data
#' @importFrom ready4show write_mdl_plt_fl
#' @importFrom MASS boxcox
#' @keywords internal
write_box_cox_tfmn <- function (data_tb, predr_var_nm_1L_chr, path_to_write_to_1L_chr, 
    depnt_var_nm_1L_chr = "aqol6d_total_w", covar_var_nms_chr = NA_character_, 
    fl_nm_pfx_1L_chr = "A_RT", height_1L_dbl = 6, width_1L_dbl = 6, 
    start_1L_chr = NULL, mdl_types_lup = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    mdl <- make_mdl(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr, 
        mdl_type_1L_chr = "OLS_NTF", mdl_types_lup = mdl_types_lup, 
        start_1L_chr = start_1L_chr)
    path_to_plot_1L_chr <- ready4show::write_mdl_plt_fl(plt_fn = MASS::boxcox, 
        fn_args_ls = list(mdl, plotit = T), path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
        plt_nm_1L_chr = paste0(fl_nm_pfx_1L_chr, "_", predr_var_nm_1L_chr, 
            "_", "BOXCOX"), height_1L_dbl = height_1L_dbl, width_1L_dbl = width_1L_dbl)
    return(path_to_plot_1L_chr)
}
#' Write bayesian regression model model plots
#' @description write_brm_model_plts() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write bayesian regression model model plots. The function returns Model plots paths (a list).
#' @param mdl_ls Model list (a list of models)
#' @param tfd_data_tb Transformed data (a tibble)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param depnt_var_desc_1L_chr Dependent variable description (a character vector of length one), Default: 'Utility score'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param tfmn_fn Transformation (a function), Default: function(x) {
#'    x
#'}
#' @param units_1L_chr Units (a character vector of length one), Default: 'in'
#' @param height_dbl Height (a double vector), Default: c(rep(6, 2), rep(5, 2))
#' @param width_dbl Width (a double vector), Default: c(rep(6, 2), rep(6, 2))
#' @param rsl_dbl Resolution (a double vector), Default: rep(300, 4)
#' @param args_ls Arguments (a list), Default: NULL
#' @param seed_1L_dbl Seed (a double vector of length one), Default: 23456
#' @return Model plots paths (a list)
#' @rdname write_brm_model_plts
#' @export 
#' @importFrom stats predict setNames
#' @importFrom purrr map discard
#' @importFrom ready4show write_mdl_plt_fl
#' @keywords internal
write_brm_model_plts <- function (mdl_ls, tfd_data_tb, mdl_nm_1L_chr, path_to_write_to_1L_chr, 
    depnt_var_nm_1L_chr = "utl_total_w", depnt_var_desc_1L_chr = "Utility score", 
    round_var_nm_1L_chr = "round", tfmn_fn = function(x) {
        x
    }, units_1L_chr = "in", height_dbl = c(rep(6, 2), rep(5, 
        2)), width_dbl = c(rep(6, 2), rep(6, 2)), rsl_dbl = rep(300, 
        4), args_ls = NULL, seed_1L_dbl = 23456) 
{
    set.seed(seed_1L_dbl)
    tfd_data_tb$Predicted <- stats::predict(mdl_ls)[, 1] %>% 
        tfmn_fn()
    plt_nms_chr <- paste0(mdl_nm_1L_chr, "_", c("coefs", "hetg", 
        "dnst", "sctr_plt"))
    mdl_plts_paths_ls <- purrr::map(1:4, ~{
        plt_fn <- fn_args_ls <- NULL
        if (.x %in% c(1, 2)) {
            plt <- plot(mdl_ls, ask = F, plot = F)
            if (length(plt) >= .x) {
                fn_args_ls <- list(mdl_ls = mdl_ls, idx_1L_int = as.integer(.x))
                plt_fn <- function(mdl_ls, idx_1L_int) {
                  plot(mdl_ls, ask = F, plot = F)[idx_1L_int]
                }
            }
        }
        else {
            if (.x == 3) {
                plt_fn <- plot_obsd_predd_dnst
                fn_args_ls <- list(tfd_data_tb = tfd_data_tb)
            }
            else {
                plt_fn <- plot_obsd_predd_sctr_cmprsn
                fn_args_ls <- list(tfd_data_tb = tfd_data_tb, 
                  depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                  depnt_var_desc_1L_chr = depnt_var_desc_1L_chr, 
                  round_var_nm_1L_chr = round_var_nm_1L_chr, 
                  args_ls = args_ls)
            }
        }
        ready4show::write_mdl_plt_fl(plt_fn, fn_args_ls = fn_args_ls, 
            path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
            plt_nm_1L_chr = plt_nms_chr[.x], units_1L_chr = units_1L_chr, 
            width_1L_dbl = width_dbl[.x], height_1L_dbl = height_dbl[.x], 
            rsl_1L_dbl = rsl_dbl[.x])
    }) %>% stats::setNames(plt_nms_chr) %>% purrr::discard(is.na)
    return(mdl_plts_paths_ls)
}
#' Write model comparison
#' @description write_mdl_cmprsn() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write model comparison. The function returns Model comparison (a list).
#' @param scored_data_tb Scored data (a tibble)
#' @param ds_smry_ls Dataset summary (a list)
#' @param mdl_smry_ls Model summary (a list)
#' @param output_data_dir_1L_chr Output data directory (a character vector of length one)
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1234
#' @return Model comparison (a list)
#' @rdname write_mdl_cmprsn
#' @export 
#' @importFrom youthvars transform_ds_for_tstng
#' @importFrom ready4fun get_from_lup_obj
write_mdl_cmprsn <- function (scored_data_tb, ds_smry_ls, mdl_smry_ls, output_data_dir_1L_chr, 
    seed_1L_int = 1234) 
{
    bl_tb <- youthvars::transform_ds_for_tstng(scored_data_tb, 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr, 
        dep_var_max_val_1L_dbl = 0.999, round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr, 
        round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr)
    ds_smry_ls$candidate_predrs_chr <- reorder_cndt_predrs_chr(ds_smry_ls$candidate_predrs_chr, 
        data_tb = bl_tb, depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr)
    mdl_smry_ls$predr_var_nm_1L_chr <- ds_smry_ls$candidate_predrs_chr[1]
    mdl_smry_ls$predr_var_desc_1L_chr <- ds_smry_ls$predictors_lup %>% 
        ready4fun::get_from_lup_obj(match_value_xx = mdl_smry_ls$predr_var_nm_1L_chr, 
            match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = "long_name_chr", 
            evaluate_lgl = F)
    mdl_smry_ls$predr_vals_dbl <- make_predr_vals(mdl_smry_ls$predr_var_nm_1L_chr, 
        candidate_predrs_lup = ds_smry_ls$predictors_lup)
    mdl_smry_ls$smry_of_sngl_predr_mdls_tb <- write_sngl_predr_multi_mdls_outps(data_tb = bl_tb, 
        folds_1L_int = mdl_smry_ls$folds_1L_int, mdl_types_chr = mdl_smry_ls$mdl_types_chr, 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        predr_var_nm_1L_chr = mdl_smry_ls$predr_var_nm_1L_chr, 
        predr_var_desc_1L_chr = mdl_smry_ls$predr_var_desc_1L_chr, 
        predr_vals_dbl = mdl_smry_ls$predr_vals_dbl, path_to_write_to_1L_chr = output_data_dir_1L_chr, 
        new_dir_nm_1L_chr = "A_Candidate_Mdls_Cmprsn", mdl_types_lup = mdl_smry_ls$mdl_types_lup, 
        dictionary_tb = ds_smry_ls$dictionary_tb)
    mdl_cmprsn_ls <- list(bl_tb = bl_tb, ds_smry_ls = ds_smry_ls, 
        mdl_smry_ls = mdl_smry_ls)
    return(mdl_cmprsn_ls)
}
#' Write model plots
#' @description write_mdl_plts() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write model plots. The function is called for its side effects and does not return a value. WARNING: This function writes R scripts to your local environment. Make sure to only use if you want this behaviour
#' @param data_tb Data (a tibble)
#' @param model_mdl Model (a model)
#' @param mdl_fl_nm_1L_chr Model file name (a character vector of length one), Default: 'OLS_NTF'
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param predr_var_desc_1L_chr Predictor variable description (a character vector of length one)
#' @param predr_vals_dbl Predictor values (a double vector)
#' @param covar_var_nms_chr Covariate variable names (a character vector), Default: 'NA'
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @param tfmn_for_bnml_1L_lgl Transformation for binomial (a logical vector of length one), Default: F
#' @param family_1L_chr Family (a character vector of length one), Default: 'NA'
#' @param plt_idxs_int Plot indices (an integer vector), Default: 1:5
#' @return NULL
#' @rdname write_mdl_plts
#' @export 
#' @importFrom purrr pwalk
#' @importFrom ready4show write_mdl_plt_fl
#' @keywords internal
write_mdl_plts <- function (data_tb, model_mdl, mdl_fl_nm_1L_chr = "OLS_NTF", depnt_var_nm_1L_chr = "utl_total_w", 
    tfmn_1L_chr = "NTF", predr_var_nm_1L_chr, predr_var_desc_1L_chr, 
    predr_vals_dbl, covar_var_nms_chr = NA_character_, path_to_write_to_1L_chr, 
    predn_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_, 
    plt_idxs_int = 1:5) 
{
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    tfd_data_tb <- transform_data_tb_for_cmprsn(data_tb, model_mdl = model_mdl, 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, predn_type_1L_chr = predn_type_1L_chr, 
        tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, family_1L_chr = family_1L_chr)
    if (1 %in% plt_idxs_int) {
        predn_ds_tb <- make_predn_ds_with_one_predr(model_mdl, 
            depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr, 
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, predr_vals_dbl = predr_vals_dbl, 
            predn_type_1L_chr = predn_type_1L_chr)
    }
    else {
        predn_ds_tb <- NULL
    }
    purrr::pwalk(list(plt_fn_ls = list(plot_lnr_cmprsn, plot_auto_lm, 
        plot_obsd_predd_dnst, plot_obsd_predd_dnst, plot_sctr_plt_cmprsn)[plt_idxs_int], 
        fn_args_ls_ls = list(list(data_tb = data_tb, predn_ds_tb = predn_ds_tb, 
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, predr_var_desc_1L_chr = predr_var_desc_1L_chr), 
            list(model_mdl, which_dbl = 1:6, ncol_1L_int = 3L, 
                label_size_1L_int = 3), list(tfd_data_tb = tfd_data_tb), 
            list(tfd_data_tb = transform_data_tb_for_cmprsn(data_tb, 
                model_mdl = model_mdl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                tf_type_1L_chr = ifelse(!4 %in% plt_idxs_int, 
                  "Predicted", "Simulated"), predn_type_1L_chr = NULL, 
                tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, 
                family_1L_chr = family_1L_chr), predd_val_var_nm_1L_chr = "Simulated"), 
            list(tfd_data_tb = tfd_data_tb))[plt_idxs_int], plt_nm_sfx_chr = c("_LNR_CMPRSN", 
            "_AUTOPLT", "_PRED_DNSTY", "_SIM_DNSTY", "_PRED_SCTR")[plt_idxs_int], 
        size_ls = list(c(6, 6), c(4, 7), c(6, 6), c(6, 6), c(6, 
            6))[plt_idxs_int]), ~ready4show::write_mdl_plt_fl(plt_fn = ..1, 
        fn_args_ls = ..2, path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
        plt_nm_1L_chr = paste0(mdl_fl_nm_1L_chr, ifelse(!is.na(covar_var_nms_chr[1]), 
            paste("_", paste0(covar_var_nms_chr[1:min(length(covar_var_nms_chr), 
                3)], collapse = "")), ""), ..3), height_1L_dbl = ..4[1], 
        width_1L_dbl = ..4[2]))
}
#' Write model type covariates models
#' @description write_mdl_type_covars_mdls() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write model type covariates models. The function returns Summary of models with covariates (a tibble).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predrs_var_nms_chr Predictors variable names (a character vector)
#' @param covar_var_nms_chr Covariate variable names (a character vector)
#' @param mdl_type_1L_chr Model type (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'D_Covars_Selection'
#' @param fl_nm_pfx_1L_chr File name prefix (a character vector of length one), Default: 'D_CT'
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param start_1L_chr Start (a character vector of length one), Default: 'NA'
#' @return Summary of models with covariates (a tibble)
#' @rdname write_mdl_type_covars_mdls
#' @export 
#' @importFrom utils data
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom caret R2
#' @importFrom dplyr pull arrange desc
#' @importFrom rlang sym
#' @importFrom stats predict AIC BIC
#' @keywords internal
write_mdl_type_covars_mdls <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predrs_var_nms_chr, 
    covar_var_nms_chr, mdl_type_1L_chr, path_to_write_to_1L_chr, 
    new_dir_nm_1L_chr = "D_Covars_Selection", fl_nm_pfx_1L_chr = "D_CT", 
    mdl_types_lup = NULL, start_1L_chr = NA_character_) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    output_dir_1L_chr <- output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    smry_of_mdls_with_covars_tb <- purrr::map_dfr(predrs_var_nms_chr, 
        ~{
            model_mdl <- make_mdl(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                predr_var_nm_1L_chr = .x, covar_var_nms_chr = covar_var_nms_chr, 
                tfmn_1L_chr = "NTF", mdl_type_1L_chr = mdl_type_1L_chr, 
                control_1L_chr = NA_character_, mdl_types_lup = mdl_types_lup, 
                start_1L_chr = start_1L_chr)
            mdl_fl_nm_1L_chr <- paste0(fl_nm_pfx_1L_chr, "_", 
                .x, "_", mdl_type_1L_chr)
            saveRDS(model_mdl, paste0(output_dir_1L_chr, "/", 
                mdl_fl_nm_1L_chr, ".RDS"))
            tibble::tibble(variable = .x, Rsquare = caret::R2(data_tb %>% 
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)), 
                stats::predict(model_mdl), form = "traditional"), 
                AIC = stats::AIC(model_mdl), BIC = stats::BIC(model_mdl), 
                Significant = paste(names(which(summary(model_mdl)$coefficients[, 
                  4] < 0.01)), collapse = " "))
        })
    smry_of_mdls_with_covars_tb <- smry_of_mdls_with_covars_tb %>% 
        dplyr::arrange(dplyr::desc(AIC))
    saveRDS(smry_of_mdls_with_covars_tb, paste0(output_dir_1L_chr, 
        "/", paste0(fl_nm_pfx_1L_chr, "_", "SMRY", "_", mdl_type_1L_chr), 
        ".RDS"))
    return(smry_of_mdls_with_covars_tb)
}
#' Write model type multi outputs
#' @description write_mdl_type_multi_outps() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write model type multi outputs. The function returns Summary of model single predictors (a tibble).
#' @param data_tb Data (a tibble)
#' @param folds_1L_int Folds (an integer vector of length one), Default: 10
#' @param predrs_var_nms_chr Predictors variable names (a character vector)
#' @param covar_var_nms_chr Covariate variable names (a character vector), Default: 'NA'
#' @param start_1L_chr Start (a character vector of length one), Default: NULL
#' @param mdl_type_1L_chr Model type (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one)
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param fl_nm_pfx_1L_chr File name prefix (a character vector of length one), Default: 'C_PREDR'
#' @param plt_idxs_int Plot indices (an integer vector), Default: c(3, 5)
#' @return Summary of model single predictors (a tibble)
#' @rdname write_mdl_type_multi_outps
#' @export 
#' @importFrom utils data
#' @importFrom purrr map_dfr
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom dplyr select mutate everything arrange desc
#' @keywords internal
write_mdl_type_multi_outps <- function (data_tb, folds_1L_int = 10, predrs_var_nms_chr, covar_var_nms_chr = NA_character_, 
    start_1L_chr = NULL, mdl_type_1L_chr, depnt_var_nm_1L_chr = "utl_total_w", 
    path_to_write_to_1L_chr, new_dir_nm_1L_chr, mdl_types_lup = NULL, 
    fl_nm_pfx_1L_chr = "C_PREDR", plt_idxs_int = c(3, 5)) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    output_dir_1L_chr <- output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    smry_of_mdl_sngl_predrs_tb <- purrr::map_dfr(predrs_var_nms_chr, 
        ~{
            tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
                target_var_nm_1L_chr = "tfmn_chr", evaluate_lgl = F)
            mdl_smry_tb <- write_mdl_type_sngl_outps(data_tb, 
                folds_1L_int = folds_1L_int, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                tfmn_1L_chr = tfmn_1L_chr, start_1L_chr = start_1L_chr, 
                predr_var_nm_1L_chr = .x, predr_var_desc_1L_chr = NA_character_, 
                predr_vals_dbl = NA_real_, covar_var_nms_chr = covar_var_nms_chr, 
                mdl_type_1L_chr = mdl_type_1L_chr, path_to_write_to_1L_chr = output_dir_1L_chr, 
                mdl_types_lup = mdl_types_lup, mdl_fl_nm_1L_chr = paste0(fl_nm_pfx_1L_chr, 
                  "_", .x, "_", mdl_type_1L_chr), plt_idxs_int = plt_idxs_int)
            if (!is.null(folds_1L_int)) {
                mdl_smry_tb <- mdl_smry_tb %>% dplyr::select((-Model)) %>% 
                  dplyr::mutate(Predictor = .x) %>% dplyr::select(Predictor, 
                  dplyr::everything())
            }
            mdl_smry_tb
        })
    if (!is.null(folds_1L_int)) {
        smry_of_mdl_sngl_predrs_tb <- smry_of_mdl_sngl_predrs_tb %>% 
            dplyr::arrange(dplyr::desc(RsquaredP))
    }
    return(smry_of_mdl_sngl_predrs_tb)
}
#' Write model type single outputs
#' @description write_mdl_type_sngl_outps() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write model type single outputs. The function returns Summary of one predictor model (a tibble).
#' @param data_tb Data (a tibble)
#' @param folds_1L_int Folds (an integer vector of length one), Default: 10
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param start_1L_chr Start (a character vector of length one), Default: NULL
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param predr_var_desc_1L_chr Predictor variable description (a character vector of length one)
#' @param predr_vals_dbl Predictor values (a double vector)
#' @param covar_var_nms_chr Covariate variable names (a character vector), Default: 'NA'
#' @param mdl_type_1L_chr Model type (a character vector of length one), Default: 'OLS_NTF'
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param mdl_fl_nm_1L_chr Model file name (a character vector of length one)
#' @param plt_idxs_int Plot indices (an integer vector), Default: NA
#' @return Summary of one predictor model (a tibble)
#' @rdname write_mdl_type_sngl_outps
#' @export 
#' @importFrom utils data
#' @importFrom purrr map_chr
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @importFrom tibble tibble
#' @keywords internal
write_mdl_type_sngl_outps <- function (data_tb, folds_1L_int = 10, depnt_var_nm_1L_chr = "utl_total_w", 
    start_1L_chr = NULL, tfmn_1L_chr = "NTF", predr_var_nm_1L_chr, 
    predr_var_desc_1L_chr, predr_vals_dbl, covar_var_nms_chr = NA_character_, 
    mdl_type_1L_chr = "OLS_NTF", mdl_types_lup = NULL, path_to_write_to_1L_chr, 
    mdl_fl_nm_1L_chr, plt_idxs_int = NA_integer_) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    arg_vals_chr <- c("control_chr", "family_chr", "predn_type_chr") %>% 
        purrr::map_chr(~ready4fun::get_from_lup_obj(mdl_types_lup, 
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
            target_var_nm_1L_chr = .x, evaluate_lgl = F))
    control_1L_chr <- arg_vals_chr[1]
    family_1L_chr <- arg_vals_chr[2]
    predn_type_1L_chr <- arg_vals_chr[3]
    if (is.na(predn_type_1L_chr)) 
        predn_type_1L_chr <- NULL
    if (is.na(plt_idxs_int[1])) {
        plt_idxs_int <- 1:5
        if (!is.na(control_1L_chr)) {
            if (control_1L_chr %>% startsWith("betareg")) 
                plt_idxs_int <- c(1, 3, 5)
        }
    }
    tfmn_for_bnml_1L_lgl <- ready4fun::get_from_lup_obj(mdl_types_lup, 
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
        target_var_nm_1L_chr = "tfmn_for_bnml_lgl", evaluate_lgl = F)
    data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_depnt_var_nm(depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)), !!rlang::sym(depnt_var_nm_1L_chr) %>% 
        calculate_dpnt_var_tfmn()))
    model_mdl <- make_mdl(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr, predr_var_nm_1L_chr = predr_var_nm_1L_chr, 
        covar_var_nms_chr = covar_var_nms_chr, mdl_type_1L_chr = mdl_type_1L_chr, 
        mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr, 
        start_1L_chr = start_1L_chr)
    write_mdl_plts(data_tb, model_mdl = model_mdl, mdl_fl_nm_1L_chr = mdl_fl_nm_1L_chr, 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, predr_var_desc_1L_chr = predr_var_desc_1L_chr, 
        predr_vals_dbl = predr_vals_dbl, covar_var_nms_chr = covar_var_nms_chr, 
        path_to_write_to_1L_chr = path_to_write_to_1L_chr, predn_type_1L_chr = predn_type_1L_chr, 
        tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, family_1L_chr = family_1L_chr, 
        plt_idxs_int = plt_idxs_int)
    if (!is.null(folds_1L_int)) {
        smry_of_one_predr_mdl_tb <- make_smry_of_mdl_outp(data_tb, 
            model_mdl = model_mdl, folds_1L_int = folds_1L_int, 
            depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr, 
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr, 
            mdl_types_lup = mdl_types_lup, predn_type_1L_chr = predn_type_1L_chr)
    }
    else {
        smry_of_one_predr_mdl_tb <- tibble::tibble()
    }
    saveRDS(model_mdl, paste0(path_to_write_to_1L_chr, "/", mdl_fl_nm_1L_chr, 
        ".RDS"))
    return(smry_of_one_predr_mdl_tb)
}
#' Write models with covariates comparison
#' @description write_mdls_with_covars_cmprsn() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write models with covariates comparison. The function returns Output summary (a list).
#' @param scored_data_tb Scored data (a tibble)
#' @param bl_tb Baseline (a tibble)
#' @param ds_smry_ls Dataset summary (a list)
#' @param mdl_smry_ls Model summary (a list)
#' @param output_data_dir_1L_chr Output data directory (a character vector of length one)
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1234
#' @param session_data_ls Session data (a list)
#' @return Output summary (a list)
#' @rdname write_mdls_with_covars_cmprsn
#' @export 

write_mdls_with_covars_cmprsn <- function (scored_data_tb, bl_tb, ds_smry_ls, mdl_smry_ls, output_data_dir_1L_chr, 
    seed_1L_int = 1234, session_data_ls) 
{
    empty_tb <- write_mdl_type_multi_outps(data_tb = bl_tb, folds_1L_int = NULL, 
        start_1L_chr = NA_character_, predrs_var_nms_chr = mdl_smry_ls$predr_cmprsn_tb$predr_chr, 
        covar_var_nms_chr = mdl_smry_ls$prefd_covars_chr, mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1], 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        path_to_write_to_1L_chr = output_data_dir_1L_chr, new_dir_nm_1L_chr = "E_Predrs_W_Covars_Sngl_Mdl_Cmprsn", 
        mdl_types_lup = mdl_smry_ls$mdl_types_lup, fl_nm_pfx_1L_chr = "E_CK_CV")
    mdl_smry_ls$predr_vars_nms_ls <- make_predr_vars_nms_ls(main_predrs_chr = mdl_smry_ls$predr_cmprsn_tb$predr_chr, 
        covars_ls = list(mdl_smry_ls$prefd_covars_chr))
    mdl_smry_ls$mdl_nms_ls <- make_mdl_nms_ls(mdl_smry_ls$predr_vars_nms_ls, 
        mdl_types_chr = mdl_smry_ls$prefd_mdl_types_chr)
    outp_smry_ls <- list(scored_data_tb = scored_data_tb, smry_of_sngl_predr_mdls_tb = mdl_smry_ls$smry_of_sngl_predr_mdls_tb, 
        prefd_mdl_types_chr = mdl_smry_ls$prefd_mdl_types_chr, 
        predr_cmprsn_tb = mdl_smry_ls$predr_cmprsn_tb, smry_of_mdl_sngl_predrs_tb = mdl_smry_ls$smry_of_mdl_sngl_predrs_tb, 
        mdls_with_covars_smry_tb = mdl_smry_ls$mdls_with_covars_smry_tb, 
        signt_covars_chr = mdl_smry_ls$signt_covars_chr, prefd_covars_chr = mdl_smry_ls$prefd_covars_chr, 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        predr_vars_nms_ls = mdl_smry_ls$predr_vars_nms_ls, mdl_nms_ls = mdl_smry_ls$mdl_nms_ls, 
        id_var_nm_1L_chr = ds_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr, 
        round_bl_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr, 
        path_to_write_to_1L_chr = output_data_dir_1L_chr, seed_1L_int = seed_1L_int, 
        folds_1L_int = mdl_smry_ls$folds_1L_int, max_nbr_of_boruta_mdl_runs_int = mdl_smry_ls$max_nbr_of_boruta_mdl_runs_int, 
        mdl_types_lup = mdl_smry_ls$mdl_types_lup, file_paths_chr = list.files(output_data_dir_1L_chr, 
            recursive = T), session_data_ls = session_data_ls)
    saveRDS(outp_smry_ls, paste0(outp_smry_ls$path_to_write_to_1L_chr, 
        "/I_ALL_OUTPUT_.RDS"))
    return(outp_smry_ls)
}
#' Write new output directory
#' @description write_new_outp_dir() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write new output directory. The function returns Output directory (a character vector of length one).
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one)
#' @return Output directory (a character vector of length one)
#' @rdname write_new_outp_dir
#' @export 

#' @keywords internal
write_new_outp_dir <- function (path_to_write_to_1L_chr, new_dir_nm_1L_chr) 
{
    output_dir_1L_chr <- paste0(path_to_write_to_1L_chr, "/", 
        new_dir_nm_1L_chr)
    if (!dir.exists(output_dir_1L_chr)) 
        dir.create(output_dir_1L_chr)
    return(output_dir_1L_chr)
}
#' Write predictor and covariates comparison
#' @description write_predr_and_covars_cmprsn() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write predictor and covariates comparison. The function returns Predictor and covariates comparison (a list).
#' @param scored_data_tb Scored data (a tibble)
#' @param bl_tb Baseline (a tibble)
#' @param ds_smry_ls Dataset summary (a list)
#' @param mdl_smry_ls Model summary (a list)
#' @param output_data_dir_1L_chr Output data directory (a character vector of length one)
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1234
#' @return Predictor and covariates comparison (a list)
#' @rdname write_predr_and_covars_cmprsn
#' @export 
#' @importFrom youthvars transform_ds_for_tstng
write_predr_and_covars_cmprsn <- function (scored_data_tb, bl_tb, ds_smry_ls, mdl_smry_ls, output_data_dir_1L_chr, 
    seed_1L_int = 1234) 
{
    mdl_smry_ls$predr_cmprsn_tb <- write_predr_cmprsn_outps(data_tb = bl_tb, 
        path_to_write_to_1L_chr = output_data_dir_1L_chr, new_dir_nm_1L_chr = "B_Candidate_Predrs_Cmprsn", 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr, 
        max_nbr_of_boruta_mdl_runs_int = mdl_smry_ls$max_nbr_of_boruta_mdl_runs_int)
    mdl_smry_ls$smry_of_mdl_sngl_predrs_tb <- write_mdl_type_multi_outps(data_tb = bl_tb, 
        folds_1L_int = mdl_smry_ls$folds_1L_int, predrs_var_nms_chr = mdl_smry_ls$predr_cmprsn_tb$predr_chr, 
        mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1], 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        path_to_write_to_1L_chr = output_data_dir_1L_chr, new_dir_nm_1L_chr = "C_Predrs_Sngl_Mdl_Cmprsn", 
        fl_nm_pfx_1L_chr = "C_PREDR", mdl_types_lup = mdl_smry_ls$mdl_types_lup)
    bl_tb <- scored_data_tb %>% youthvars::transform_ds_for_tstng(candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr, 
        covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr, 
        remove_all_msng_1L_lgl = T, round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr)
    mdl_smry_ls$mdls_with_covars_smry_tb <- write_mdl_type_covars_mdls(bl_tb, 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        predrs_var_nms_chr = ds_smry_ls$candidate_predrs_chr, 
        covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr, 
        mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1], 
        path_to_write_to_1L_chr = output_data_dir_1L_chr, new_dir_nm_1L_chr = "D_Predr_Covars_Cmprsn", 
        fl_nm_pfx_1L_chr = "D_CT", mdl_types_lup = mdl_smry_ls$mdl_types_lup)
    mdl_smry_ls$signt_covars_chr <- get_signft_covars(mdls_with_covars_smry_tb = mdl_smry_ls$mdls_with_covars_smry_tb, 
        covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr)
    predr_and_covars_cmprsn_ls <- list(bl_tb = bl_tb, ds_smry_ls = ds_smry_ls, 
        mdl_smry_ls = mdl_smry_ls)
    return(predr_and_covars_cmprsn_ls)
}
#' Write predictor and model testing results
#' @description write_predr_and_mdl_tstng_results() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write predictor and model testing results. The function returns Output summary (a list).
#' @param scored_data_tb Scored data (a tibble)
#' @param ds_smry_ls Dataset summary (a list)
#' @param mdl_smry_ls Model summary (a list)
#' @param session_data_ls Session data (a list)
#' @param output_data_dir_1L_chr Output data directory (a character vector of length one)
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1234
#' @return Output summary (a list)
#' @rdname write_predr_and_mdl_tstng_results
#' @export 

#' @keywords internal
write_predr_and_mdl_tstng_results <- function (scored_data_tb, ds_smry_ls, mdl_smry_ls, session_data_ls, 
    output_data_dir_1L_chr, seed_1L_int = 1234) 
{
    cmprsn_ls <- write_mdl_cmprsn(scored_data_tb = scored_data_tb, 
        ds_smry_ls = ds_smry_ls, mdl_smry_ls = mdl_smry_ls, output_data_dir_1L_chr = output_data_dir_1L_chr, 
        seed_1L_int = seed_1L_int)
    if (ifelse(is.null(cmprsn_ls$mdl_smry_ls$prefd_mdl_types_chr), 
        T, is.na(cmprsn_ls$mdl_smry_ls$prefd_mdl_types_chr[1]))) {
        cmprsn_ls$mdl_smry_ls$prefd_mdl_types_chr <- make_prefd_mdls_vec(cmprsn_ls$mdl_smry_ls$smry_of_sngl_predr_mdls_tb, 
            choose_from_pfx_chr = cmprsn_ls$mdl_smry_ls$choose_from_pfx_chr)
    }
    cmprsn_ls <- write_predr_and_covars_cmprsn(scored_data_tb = scored_data_tb, 
        bl_tb = cmprsn_ls$bl_tb, ds_smry_ls = cmprsn_ls$ds_smry_ls, 
        mdl_smry_ls = cmprsn_ls$mdl_smry_ls, output_data_dir_1L_chr = output_data_dir_1L_chr, 
        seed_1L_int = seed_1L_int)
    if (ifelse(is.null(cmprsn_ls$mdl_smry_ls$prefd_covars_chr), 
        T, is.na(cmprsn_ls$mdl_smry_ls$prefd_covars_chr[1]))) {
        cmprsn_ls$mdl_smry_ls$prefd_covars_chr <- cmprsn_ls$mdl_smry_ls$signt_covars_chr
    }
    outp_smry_ls <- write_mdls_with_covars_cmprsn(scored_data_tb = scored_data_tb, 
        bl_tb = cmprsn_ls$bl_tb, ds_smry_ls = cmprsn_ls$ds_smry_ls, 
        mdl_smry_ls = cmprsn_ls$mdl_smry_ls, output_data_dir_1L_chr = output_data_dir_1L_chr, 
        seed_1L_int = seed_1L_int, session_data_ls = session_data_ls)
    return(outp_smry_ls)
}
#' Write predictor comparison outputs
#' @description write_predr_cmprsn_outps() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write predictor comparison outputs. The function returns Confirmed predictors (a tibble).
#' @param data_tb Data (a tibble)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'B_Candidate_Predrs_Cmprsn'
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param candidate_predrs_chr Candidate predictors (a character vector)
#' @param max_nbr_of_boruta_mdl_runs_int Maximum number of boruta model runs (an integer vector), Default: 300
#' @return Confirmed predictors (a tibble)
#' @rdname write_predr_cmprsn_outps
#' @export 
#' @importFrom randomForest randomForest varImpPlot
#' @importFrom stats as.formula
#' @importFrom Boruta Boruta
#' @importFrom purrr pwalk
#' @importFrom ready4show write_mdl_plt_fl
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange desc filter
#' @keywords internal
write_predr_cmprsn_outps <- function (data_tb, path_to_write_to_1L_chr, new_dir_nm_1L_chr = "B_Candidate_Predrs_Cmprsn", 
    depnt_var_nm_1L_chr = "utl_total_w", candidate_predrs_chr, 
    max_nbr_of_boruta_mdl_runs_int = 300L) 
{
    if (length(candidate_predrs_chr) > 1) {
        covar_var_nms_chr <- candidate_predrs_chr[2:length(candidate_predrs_chr)]
    }
    else {
        covar_var_nms_chr <- NA_character_
    }
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_var_nm_1L_chr = candidate_predrs_chr[1], covar_var_nms_chr = covar_var_nms_chr)
    rf_mdl <- randomForest::randomForest(stats::as.formula(paste0(depnt_var_nm_1L_chr, 
        " ~ .")), data = data_tb, importance = TRUE)
    boruta_mdl <- Boruta::Boruta(stats::as.formula(paste0(depnt_var_nm_1L_chr, 
        " ~ .")), data = data_tb, maxRuns = max_nbr_of_boruta_mdl_runs_int)
    output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    purrr::pwalk(list(fn_ls = list(randomForest::varImpPlot, 
        plot), fn_args_ls_ls = list(list(rf_mdl, main = ""), 
        list(boruta_mdl, cex = 1.5, cex.axis = 0.8, las = 2, 
            xlab = "", main = "")), plt_nm_sfx_chr = c("_RF_VAR_IMP", 
        "_BORUTA_VAR_IMP"), size_ls = list(c(6, 6), c(4, 6))), 
        ~ready4show::write_mdl_plt_fl(plt_fn = ..1, fn_args_ls = ..2, 
            path_to_write_to_1L_chr = output_dir_1L_chr, plt_nm_1L_chr = paste0("B_PRED_CMPRSN", 
                ..3), height_1L_dbl = ..4[1], width_1L_dbl = ..4[2]))
    confirmed_predrs_chr <- names(boruta_mdl$finalDecision)[boruta_mdl$finalDecision == 
        "Confirmed"]
    confirmed_predrs_tb <- rf_mdl$importance %>% tibble::as_tibble(rownames = "predr_chr") %>% 
        dplyr::arrange(dplyr::desc(`%IncMSE`)) %>% dplyr::filter(predr_chr %in% 
        confirmed_predrs_chr)
    return(confirmed_predrs_tb)
}
#' Write results to comma separated variables file
#' @description write_results_to_csv() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write results to comma separated variables file. The function returns Datasets (a tibble).
#' @param synth_data_spine_ls Synthetic data spine (a list)
#' @param output_dir_1L_chr Output directory (a character vector of length one), Default: '.'
#' @return Datasets (a tibble)
#' @rdname write_results_to_csv
#' @export 
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map_dfr map_dfc map map_dbl walk2
#' @importFrom stats setNames
#' @importFrom dplyr mutate select everything
#' @keywords internal
write_results_to_csv <- function (synth_data_spine_ls, output_dir_1L_chr = ".") 
{
    measurements_tb <- tibble::tibble(timepoint_nms_chr = synth_data_spine_ls$timepoint_nms_chr, 
        nbr_obs_dbl = synth_data_spine_ls$nbr_obs_dbl)
    var_summ_res_tb <- suppressMessages(purrr::map_dfr(1:length(synth_data_spine_ls$timepoint_nms_chr), 
        ~{
            idx_dbl <- .x
            suppressWarnings({
                synth_data_spine_ls[c(4:6)] %>% purrr::map_dfc(~.x[idx_dbl])
            }) %>% stats::setNames(c("Mean", "SD", "N_Missing")) %>% 
                dplyr::mutate(var_names_chr = synth_data_spine_ls$var_names_chr, 
                  timepoint_nms_chr = synth_data_spine_ls$timepoint_nms_chr[idx_dbl]) %>% 
                dplyr::select(timepoint_nms_chr, var_names_chr, 
                  dplyr::everything())
        }))
    cor_tb_ls <- synth_data_spine_ls$cor_mat_ls %>% purrr::map(~tibble::as_tibble(.x) %>% 
        stats::setNames(synth_data_spine_ls$var_names_chr) %>% 
        dplyr::mutate(var_names_chr = synth_data_spine_ls$var_names_chr) %>% 
        dplyr::select(var_names_chr, dplyr::everything())) %>% 
        stats::setNames(paste0(synth_data_spine_ls$timepoint_nms_chr, 
            "_correlations_tb"))
    var_class_pars_tb <- synth_data_spine_ls[7:9] %>% tibble::as_tibble() %>% 
        dplyr::mutate(min_dbl = purrr::map_dbl(min_max_ls, ~.x[1]), 
            max_dbl = purrr::map_dbl(min_max_ls, ~.x[2])) %>% 
        dplyr::select(var_names_chr, dplyr::everything(), -min_max_ls)
    output_ls <- list(measurements_tb = measurements_tb, var_summ_res_tb = var_summ_res_tb, 
        var_class_pars_tb = var_class_pars_tb) %>% append(cor_tb_ls)
    dss_tb <- tibble::tibble(ds_obj_nm_chr = names(output_ls), 
        title_chr = c("Brief summary table of the number of observations for which data was collected at each study timepoint.", 
            "Summary statistics (Mean, SD and Number Missing) for AQoL6D health utility and six mental health outcome measures for each study timepoint.", 
            "Brief information about the data structure (whether discrete and allowable range) of AQoL6D health utility and six mental health outcome variables.", 
            paste0("Correlation matrix for AQoL6D health utility and six mental health outcome measures at the ", 
                synth_data_spine_ls$timepoint_nms_chr, " study timepoint.")))
    purrr::walk2(output_ls, names(output_ls), ~write.csv(.x, 
        file = paste0(output_dir_1L_chr, "/", .y, ".csv"), row.names = F))
    return(dss_tb)
}
#' Write shareable models
#' @description write_shareable_mdls() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write shareable models. The function returns Output summary (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'G_Shareable'
#' @param shareable_title_detail_1L_chr Shareable title detail (a character vector of length one), Default: ''
#' @return Output summary (a list)
#' @rdname write_shareable_mdls
#' @export 
#' @importFrom purrr map_chr flatten_chr map map_lgl map_int map2
#' @importFrom stringr str_locate
#' @importFrom dplyr filter rename
#' @importFrom rlang sym
#' @importFrom stats setNames
write_shareable_mdls <- function (outp_smry_ls, new_dir_nm_1L_chr = "G_Shareable", shareable_title_detail_1L_chr = "") 
{
    output_dir_1L_chr <- write_new_outp_dir(outp_smry_ls$path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    incld_mdl_paths_chr <- outp_smry_ls$file_paths_chr %>% purrr::map_chr(~{
        file_path_1L_chr <- .x
        mdl_fl_nms_chr <- paste0(outp_smry_ls$mdl_nms_ls %>% 
            purrr::flatten_chr(), ".RDS")
        mdl_fl_nms_locn_ls <- mdl_fl_nms_chr %>% purrr::map(~stringr::str_locate(file_path_1L_chr, 
            .x))
        match_lgl <- mdl_fl_nms_locn_ls %>% purrr::map_lgl(~!(is.na(.x[[1, 
            1]]) | is.na(.x[[1, 2]])))
        if (any(match_lgl)) {
            file_path_1L_chr
        }
        else {
            NA_character_
        }
    })
    incld_mdl_paths_chr <- incld_mdl_paths_chr[!is.na(incld_mdl_paths_chr)]
    ranked_mdl_nms_chr <- outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr()
    sorted_mdl_nms_chr <- sort(ranked_mdl_nms_chr)
    rank_idcs_int <- purrr::map_int(sorted_mdl_nms_chr, ~which(ranked_mdl_nms_chr == 
        .x))
    incld_mdl_paths_chr <- incld_mdl_paths_chr[order(rank_idcs_int)]
    shareable_mdls_ls <- outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr() %>% 
        purrr::map2(incld_mdl_paths_chr, ~{
            model_mdl <- readRDS(paste0(outp_smry_ls$path_to_write_to_1L_chr, 
                "/", .y))
            mdl_smry_tb <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Model == 
                .x)
            data_tb <- model_mdl$data
            if (endsWith(.x, "OLS_CLL")) 
                data_tb <- data_tb %>% dplyr::rename(`:=`(!!rlang::sym(paste0(outp_smry_ls$depnt_var_nm_1L_chr, 
                  "_CLL")), !!rlang::sym(paste0(outp_smry_ls$depnt_var_nm_1L_chr, 
                  "_cloglog"))))
            shareable_mdl <- make_shareable_mdl(data_tb = data_tb, 
                mdl_smry_tb = mdl_smry_tb, depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
                tfmn_1L_chr = ifelse(endsWith(.x, "OLS_CLL"), 
                  "CLL", "NTF"), mdl_type_1L_chr = ifelse(endsWith(.x, 
                  "OLS_CLL"), "OLS_CLL", "GLM_GSN_LOG"), mdl_types_lup = outp_smry_ls$mdl_types_lup, 
                control_1L_chr = NA_character_, start_1L_chr = NA_character_, 
                seed_1L_int = outp_smry_ls$seed_1L_int)
            saveRDS(shareable_mdl, paste0(output_dir_1L_chr, 
                "/", .x, ".RDS"))
            shareable_mdl
        }) %>% stats::setNames(outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr())
    outp_smry_ls$shareable_mdls_ls <- shareable_mdls_ls
    outp_smry_ls$shareable_mdls_tb <- NULL
    if (!is.null(outp_smry_ls$dv_ls)) {
        outp_smry_ls$shareable_mdls_tb <- write_shareable_mdls_to_dv(outp_smry_ls, 
            new_dir_nm_1L_chr = new_dir_nm_1L_chr, shareable_title_detail_1L_chr = shareable_title_detail_1L_chr)
    }
    return(outp_smry_ls)
}
#' Write shareable models to dataverse
#' @description write_shareable_mdls_to_dv() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write shareable models to dataverse. The function returns Shareable models (a tibble).
#' @param outp_smry_ls Output summary (a list)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'G_Shareable'
#' @param shareable_title_detail_1L_chr Shareable title detail (a character vector of length one), Default: ''
#' @return Shareable models (a tibble)
#' @rdname write_shareable_mdls_to_dv
#' @export 
#' @importFrom tibble tibble
#' @importFrom ready4use write_fls_to_dv_ds get_fl_id_from_dv_ls
#' @importFrom dataverse get_dataset
#' @importFrom dplyr mutate
#' @importFrom purrr map_int
#' @keywords internal
write_shareable_mdls_to_dv <- function (outp_smry_ls, new_dir_nm_1L_chr = "G_Shareable", shareable_title_detail_1L_chr = "") 
{
    output_dir_1L_chr <- write_new_outp_dir(outp_smry_ls$path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    shareable_mdls_tb <- tibble::tibble(ds_obj_nm_chr = names(outp_smry_ls$shareable_mdls_ls), 
        title_chr = paste0("A shareable (contains no confidential data) statistical model, ", 
            names(outp_smry_ls$shareable_mdls_ls), ".", shareable_title_detail_1L_chr))
    ready4use::write_fls_to_dv_ds(shareable_mdls_tb, dv_nm_1L_chr = outp_smry_ls$dv_ls$dv_nm_1L_chr, 
        ds_url_1L_chr = outp_smry_ls$dv_ls$ds_url_1L_chr, parent_dv_dir_1L_chr = outp_smry_ls$dv_ls$parent_dv_dir_1L_chr, 
        paths_to_dirs_chr = output_dir_1L_chr, paths_are_rltv_1L_lgl = F, 
        inc_fl_types_chr = ".RDS")
    ds_ls <- dataverse::get_dataset(outp_smry_ls$dv_ls$ds_url_1L_chr)
    shareable_mdls_tb <- shareable_mdls_tb %>% dplyr::mutate(dv_nm_chr = outp_smry_ls$dv_ls$dv_nm_1L_chr, 
        fl_ids_int = ds_obj_nm_chr %>% purrr::map_int(~ready4use::get_fl_id_from_dv_ls(ds_ls, 
            fl_nm_1L_chr = paste0(.x, ".RDS")) %>% as.integer()))
    return(shareable_mdls_tb)
}
#' Write single predictor multi models outputs
#' @description write_sngl_predr_multi_mdls_outps() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write single predictor multi models outputs. The function returns Summary of single predictor models (a tibble).
#' @param data_tb Data (a tibble)
#' @param mdl_types_chr Model types (a character vector)
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param predr_var_desc_1L_chr Predictor variable description (a character vector of length one)
#' @param predr_vals_dbl Predictor values (a double vector)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'A_Candidate_Mdls_Cmprsn'
#' @param start_1L_chr Start (a character vector of length one), Default: NULL
#' @param covar_var_nms_chr Covariate variable names (a character vector), Default: 'NA'
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param folds_1L_int Folds (an integer vector of length one), Default: 10
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param fl_nm_pfx_1L_chr File name prefix (a character vector of length one), Default: 'A_RT_'
#' @param plt_idxs_int Plot indices (an integer vector), Default: NA
#' @param dictionary_tb Dictionary (a tibble)
#' @return Summary of single predictor models (a tibble)
#' @rdname write_sngl_predr_multi_mdls_outps
#' @export 
#' @importFrom utils data
#' @importFrom ready4show write_mdl_plt_fl
#' @importFrom purrr map_dfr
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom dplyr arrange desc
#' @keywords internal
write_sngl_predr_multi_mdls_outps <- function (data_tb, mdl_types_chr, predr_var_nm_1L_chr, predr_var_desc_1L_chr, 
    predr_vals_dbl, path_to_write_to_1L_chr, new_dir_nm_1L_chr = "A_Candidate_Mdls_Cmprsn", 
    start_1L_chr = NULL, covar_var_nms_chr = NA_character_, depnt_var_nm_1L_chr = "utl_total_w", 
    folds_1L_int = 10, mdl_types_lup = NULL, fl_nm_pfx_1L_chr = "A_RT_", 
    plt_idxs_int = NA_integer_, dictionary_tb) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    ready4show::write_mdl_plt_fl(plt_fn = make_tfmn_cmprsn_plt, 
        fn_args_ls = list(data_tb = data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
            dictionary_tb = dictionary_tb), path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
        plt_nm_1L_chr = "A_TFMN_CMPRSN_DNSTY", height_1L_dbl = 6, 
        width_1L_dbl = 10)
    smry_of_sngl_predr_mdls_tb <- purrr::map_dfr(mdl_types_chr, 
        ~{
            tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .x, 
                target_var_nm_1L_chr = "tfmn_chr", evaluate_lgl = F)
            write_mdl_type_sngl_outps(data_tb, folds_1L_int = folds_1L_int, 
                depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, start_1L_chr = start_1L_chr, 
                tfmn_1L_chr = tfmn_1L_chr, predr_var_nm_1L_chr = predr_var_nm_1L_chr, 
                predr_var_desc_1L_chr = predr_var_desc_1L_chr, 
                predr_vals_dbl = predr_vals_dbl, covar_var_nms_chr = covar_var_nms_chr, 
                mdl_type_1L_chr = .x, path_to_write_to_1L_chr = output_dir_1L_chr, 
                mdl_types_lup = mdl_types_lup, mdl_fl_nm_1L_chr = paste0(fl_nm_pfx_1L_chr, 
                  predr_var_nm_1L_chr, "_", .x), plt_idxs_int = plt_idxs_int)
        })
    if (!is.null(folds_1L_int)) 
        smry_of_sngl_predr_mdls_tb <- smry_of_sngl_predr_mdls_tb %>% 
            dplyr::arrange(dplyr::desc(RsquaredP))
    return(smry_of_sngl_predr_mdls_tb)
}
#' Write time series models
#' @description write_ts_mdls() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write time series models. The function returns Models summary (a tibble).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_vars_nms_ls Predictor variables names (a list)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one), Default: 'Baseline'
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param fn_ls Function list (a list of functions)
#' @param mdl_nms_ls Model names (a list)
#' @param mdl_smry_dir_1L_chr Model summary directory (a character vector of length one)
#' @param predictors_lup Predictors (a lookup table)
#' @param iters_1L_int Iterations (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Models summary (a tibble)
#' @rdname write_ts_mdls
#' @export 
#' @importFrom purrr map_dfr map2_dfr
#' @keywords internal
write_ts_mdls <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_ls, 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline", backend_1L_chr = getOption("brms.backend", 
        "rstan"), fn_ls, mdl_nms_ls, mdl_smry_dir_1L_chr, predictors_lup, 
    iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    if (!dir.exists(mdl_smry_dir_1L_chr)) 
        dir.create(mdl_smry_dir_1L_chr)
    mdls_smry_tb <- purrr::map_dfr(1:length(mdl_nms_ls), ~{
        idx_1L_int <- .x
        purrr::map2_dfr(fn_ls, mdl_nms_ls[[idx_1L_int]], ~{
            smry_ls <- make_smry_of_ts_mdl_outp(data_tb = data_tb, 
                fn = .x, predr_vars_nms_chr = predr_vars_nms_ls[[idx_1L_int]], 
                mdl_nm_1L_chr = .y, path_to_write_to_1L_chr = mdl_smry_dir_1L_chr, 
                depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
                round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr, 
                predictors_lup = predictors_lup, backend_1L_chr = backend_1L_chr, 
                iters_1L_int = iters_1L_int, seed_1L_int = seed_1L_int)
            Sys.sleep(5)
            smry_ls$smry_of_ts_mdl_tb
        })
    })
    saveRDS(mdls_smry_tb, paste0(mdl_smry_dir_1L_chr, "/mdls_smry_tb.RDS"))
    return(mdls_smry_tb)
}
#' Write time series models from algorithm output
#' @description write_ts_mdls_from_alg_outp() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write time series models from algorithm output. The function returns Output summary (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param fn_ls Function list (a list of functions)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'F_TS_Mdls'
#' @param predictors_lup Predictors (a lookup table)
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iterations (an integer vector of length one), Default: 4000
#' @return Output summary (a list)
#' @rdname write_ts_mdls_from_alg_outp
#' @export 

write_ts_mdls_from_alg_outp <- function (outp_smry_ls, fn_ls, new_dir_nm_1L_chr = "F_TS_Mdls", 
    predictors_lup, backend_1L_chr = getOption("brms.backend", 
        "rstan"), iters_1L_int = 4000L) 
{
    output_dir_1L_chr <- output_dir_1L_chr <- write_new_outp_dir(outp_smry_ls$path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    mdls_smry_tb <- write_ts_mdls(data_tb = outp_smry_ls$scored_data_tb, 
        depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
        predr_vars_nms_ls = outp_smry_ls$predr_vars_nms_ls, id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
        round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
        round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr, 
        fn_ls = fn_ls, mdl_nms_ls = outp_smry_ls$mdl_nms_ls, 
        mdl_smry_dir_1L_chr = output_dir_1L_chr, predictors_lup = predictors_lup, 
        backend_1L_chr = backend_1L_chr, iters_1L_int = iters_1L_int, 
        seed_1L_int = outp_smry_ls$seed_1L_int)
    outp_smry_ls$mdls_smry_tb <- mdls_smry_tb
    outp_smry_ls$file_paths_chr <- list.files(outp_smry_ls$path_to_write_to_1L_chr, 
        recursive = T)
    return(outp_smry_ls)
}
