#' Write analyses
#' @description write_analyses() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write analyses. The function is called for its side effects and does not return a value. WARNING: This function writes R scripts to your local environment. Make sure to only use if you want this behaviour
#' @param input_params_ls Input parameters (a list)
#' @param abstract_args_ls Abstract arguments (a list), Default: NULL
#' @param start_at_int Start at (an integer vector), Default: c(2, 1)
#' @return NULL
#' @rdname write_analyses
#' @export 
#' @importFrom purrr walk pluck
write_analyses <- function (input_params_ls, abstract_args_ls = NULL, start_at_int = c(2, 
    1)) 
{
    write_report(params_ls = input_params_ls$params_ls, paths_ls = input_params_ls$path_params_ls$paths_ls, 
        rprt_nm_1L_chr = "AAA_PMRY_ANLYS_MTH", abstract_args_ls = abstract_args_ls, 
        header_yaml_args_ls = input_params_ls$header_yaml_args_ls)
    if (!is.null(input_params_ls$scndry_anlys_params_ls)) {
        references_int <- 1:length(input_params_ls$scndry_anlys_params_ls)
        references_int %>% purrr::walk(~{
            changes_ls <- input_params_ls$scndry_anlys_params_ls %>% 
                purrr::pluck(.x)
            if (is.null(changes_ls$candidate_covar_nms_chr)) 
                changes_ls$candidate_covar_nms_chr <- input_params_ls$params_ls$candidate_covar_nms_chr %>% 
                  transform_names(input_params_ls$rename_lup, 
                    invert_1L_lgl = T)
            if (is.null(changes_ls$candidate_predrs_chr)) {
                changes_ls$candidate_covar_nms_chr <- changes_ls$candidate_covar_nms_chr[!changes_ls$candidate_covar_nms_chr %in% 
                  changes_ls$candidate_predrs_chr]
            }
            write_scndry_analysis(valid_params_ls_ls = input_params_ls, 
                candidate_covar_nms_chr = changes_ls$candidate_covar_nms_chr, 
                candidate_predrs_chr = changes_ls$candidate_predrs_chr, 
                header_yaml_args_ls = input_params_ls$header_yaml_args_ls, 
                path_params_ls = input_params_ls$path_params_ls, 
                prefd_covars_chr = changes_ls$prefd_covars_chr, 
                reference_1L_int = .x, start_at_int = start_at_int, 
                rprt_nm_1L_chr = "AAA_SUPLRY_ANLYS_MTH", abstract_args_ls = abstract_args_ls)
        })
    }
}
#' Write box cox transformation
#' @description write_box_cox_tfmn() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write box cox transformation. The function returns Path to plot (a character vector of length one).
#' @param data_tb Data (a tibble)
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
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
    depnt_var_nm_1L_chr = "utl_total_w", covar_var_nms_chr = NA_character_, 
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
#' Write complete study program output
#' @description write_csp_output() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write complete study program output. The function is called for its side effects and does not return a value. WARNING: This function writes R scripts to your local environment. Make sure to only use if you want this behaviour
#' @param path_to_csp_1L_chr Path to complete study program (a character vector of length one)
#' @param dv_ds_doi_1L_chr Dataverse dataset digital object identifier (a character vector of length one), Default: NULL
#' @param execute_1L_lgl Execute (a logical vector of length one), Default: T
#' @return NULL
#' @rdname write_csp_output
#' @export 
#' @importFrom purrr map_chr
#' @importFrom stringr str_sub
#' @importFrom knitr purl
#' @importFrom DescTools SplitPath
#' @importFrom rmarkdown render
#' @importFrom dataverse add_dataset_file
write_csp_output <- function (path_to_csp_1L_chr, dv_ds_doi_1L_chr = NULL, execute_1L_lgl = T) 
{
    readLines(path_to_csp_1L_chr) %>% purrr::map_chr(~ifelse(.x == 
        "knitr::opts_chunk$set(eval = F)", "knitr::opts_chunk$set(eval = T)", 
        .x)) %>% writeLines(con = path_to_csp_1L_chr)
    path_to_r_script_1L_chr <- stringr::str_sub(path_to_csp_1L_chr, 
        end = -3)
    knitr::purl(path_to_csp_1L_chr, path_to_r_script_1L_chr)
    readLines(path_to_csp_1L_chr) %>% purrr::map_chr(~ifelse(.x == 
        "knitr::opts_chunk$set(eval = T)", "knitr::opts_chunk$set(eval = F)", 
        .x)) %>% writeLines(con = path_to_csp_1L_chr)
    if (execute_1L_lgl) {
        old_wd_1L_chr <- getwd()
        path_info_ls <- DescTools::SplitPath(path_to_r_script_1L_chr)
        setwd(path_info_ls$dirname)
        source(path_info_ls$fullfilename)
        setwd(old_wd_1L_chr)
    }
    rmarkdown::render(path_to_csp_1L_chr)
    if (!is.null(dv_ds_doi_1L_chr)) {
        dataverse::add_dataset_file(paste0(stringr::str_sub(path_to_csp_1L_chr, 
            end = -4), "pdf"), dataset = dv_ds_doi_1L_chr, description = "Methods Report 1: Complete Study Program")
    }
}
#' Write main output directory
#' @description write_main_outp_dir() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write main output directory. The function returns Paths (a list).
#' @param params_ls Parameters (a list), Default: NULL
#' @param use_fake_data_1L_lgl Use fake data (a logical vector of length one), Default: F
#' @param R_fl_nm_1L_chr R file name (a character vector of length one), Default: 'aaaaaaaaaa.txt'
#' @return Paths (a list)
#' @rdname write_main_outp_dir
#' @export 
#' @importFrom purrr pluck
#' @importFrom ready4show make_paths_ls
#' @importFrom here i_am here
#' @importFrom youthvars write_all_outp_dirs
write_main_outp_dir <- function (params_ls = NULL, use_fake_data_1L_lgl = F, R_fl_nm_1L_chr = "aaaaaaaaaa.txt") 
{
    file.create(R_fl_nm_1L_chr)
    R_fl_nm_1L_chr <- list.files() %>% purrr::pluck(1)
    paths_ls <- ready4show::make_paths_ls(append(params_ls, list(use_fake_data_1L_lgl = use_fake_data_1L_lgl)), 
        depth_1L_int = 0)
    paths_ls$path_to_current_1L_chr <- ifelse(!is.null(paths_ls$path_to_current_1L_chr), 
        paths_ls$path_to_current_1L_chr, params_ls$path_to_current_1L_chr)
    here::i_am(paste0(paths_ls$path_from_top_level_1L_chr, "/", 
        paths_ls$path_to_current_1L_chr, "/", R_fl_nm_1L_chr))
    dir.create(paste0(here::here(paths_ls$path_from_top_level_1L_chr), 
        "/", paths_ls$write_to_dir_nm_1L_chr))
    paths_ls$R_fl_nm_1L_chr <- R_fl_nm_1L_chr
    paths_ls <- youthvars::write_all_outp_dirs(paths_ls)
    return(paths_ls)
}
#' Write manuscript
#' @description write_manuscript() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write manuscript. The function returns Results (a list).
#' @param abstract_args_ls Abstract arguments (a list), Default: NULL
#' @param input_params_ls Input parameters (a list), Default: NULL
#' @param results_ls Results (a list), Default: NULL
#' @param figures_in_body_lgl Figures in body (a logical vector), Default: NULL
#' @param output_type_1L_chr Output type (a character vector of length one), Default: NULL
#' @param tables_in_body_lgl Tables in body (a logical vector), Default: NULL
#' @param title_1L_chr Title (a character vector of length one), Default: 'Scientific manuscript'
#' @param version_1L_chr Version (a character vector of length one), Default: '0.5'
#' @param write_to_dv_1L_lgl Write to dataverse (a logical vector of length one), Default: F
#' @return Results (a list)
#' @rdname write_manuscript
#' @export 
#' @importFrom utils unzip
#' @importFrom ready4show write_header_fls
#' @importFrom rmarkdown render
#' @importFrom ready4 write_to_dv_with_wait
#' @importFrom tibble tibble
write_manuscript <- function (abstract_args_ls = NULL, input_params_ls = NULL, results_ls = NULL, 
    figures_in_body_lgl = NULL, output_type_1L_chr = NULL, tables_in_body_lgl = NULL, 
    title_1L_chr = "Scientific manuscript", version_1L_chr = "0.5", 
    write_to_dv_1L_lgl = F) 
{
    mkdn_data_dir_1L_chr <- ifelse(!is.null(input_params_ls), 
        input_params_ls$path_params_ls$paths_ls$mkdn_data_dir_1L_chr, 
        results_ls$path_params_ls$paths_ls$mkdn_data_dir_1L_chr)
    outp_dir_1L_chr <- ifelse(!is.null(input_params_ls), input_params_ls$path_params_ls$paths_ls$output_data_dir_1L_chr, 
        results_ls$path_params_ls$paths_ls$output_data_dir_1L_chr)
    output_type_1L_chr <- ifelse(!is.null(output_type_1L_chr), 
        output_type_1L_chr, ifelse(!is.null(input_params_ls), 
            input_params_ls$output_format_ls$manuscript_outp_1L_chr, 
            results_ls$output_format_ls$manuscript_outp_1L_chr))
    path_to_ms_mkdn_1L_dir <- paste0(mkdn_data_dir_1L_chr, "/ttu_lng_ss-", 
        version_1L_chr)
    path_to_results_dir_1L_chr <- ifelse(!is.null(input_params_ls), 
        input_params_ls$path_params_ls$paths_ls$reports_dir_1L_chr, 
        results_ls$path_params_ls$paths_ls$reports_dir_1L_chr)
    if (!dir.exists(path_to_ms_mkdn_1L_dir)) {
        tmp_fl <- tempfile()
        download.file(paste0("https://github.com/ready4-dev/ttu_lng_ss/archive/refs/tags/v", 
            version_1L_chr, ".zip"), tmp_fl)
        utils::unzip(tmp_fl, exdir = mkdn_data_dir_1L_chr)
        unlink(tmp_fl)
    }
    if (!is.null(input_params_ls)) {
        header_yaml_args_ls <- input_params_ls$header_yaml_args_ls
    }
    else {
        header_yaml_args_ls <- results_ls$header_yaml_args_ls
    }
    if (is.null(results_ls)) {
        results_ls <- make_results_ls(dv_ds_nm_and_url_chr = input_params_ls$path_params_ls$dv_ds_nm_and_url_chr, 
            output_format_ls = input_params_ls$output_format_ls, 
            params_ls_ls = input_params_ls, path_params_ls = input_params_ls$path_params_ls, 
            study_descs_ls = input_params_ls$study_descs_ls, 
            var_nm_change_lup = input_params_ls$study_descs_ls$var_nm_change_lup, 
            version_1L_chr = version_1L_chr)
    }
    if (is.null(abstract_args_ls)) {
        abstract_args_ls <- make_abstract_args_ls(results_ls)
    }
    ready4show::write_header_fls(path_to_header_dir_1L_chr = paste0(path_to_ms_mkdn_1L_dir, 
        "/Header"), header_yaml_args_ls = header_yaml_args_ls, 
        abstract_args_ls = abstract_args_ls)
    params_ls <- list(output_type_1L_chr = output_type_1L_chr, 
        results_ls = results_ls)
    if (!is.null(figures_in_body_lgl)) 
        params_ls$figures_in_body_lgl <- figures_in_body_lgl
    if (!is.null(tables_in_body_lgl)) 
        params_ls$tables_in_body_lgl <- tables_in_body_lgl
    rmarkdown::render(paste0(path_to_ms_mkdn_1L_dir, "/", output_type_1L_chr, 
        "/", output_type_1L_chr, ".Rmd"), output_format = NULL, 
        params = params_ls, output_file = paste0("TTU_Study_Manuscript", 
            ifelse(output_type_1L_chr == "Word", ".docx", ".pdf")), 
        output_dir = path_to_results_dir_1L_chr)
    if (write_to_dv_1L_lgl) {
        if (!is.null(input_params_ls)) {
            paths_ls <- input_params_ls$path_params_ls$paths_ls
        }
        else {
            paths_ls <- results_ls$path_params_ls$paths_ls
        }
        ready4::write_to_dv_with_wait(dss_tb = tibble::tibble(ds_obj_nm_chr = "TTU_Study_Manuscript", 
            title_chr = title_1L_chr), dv_nm_1L_chr = ifelse(!is.null(input_params_ls), 
            input_params_ls$path_params_ls$dv_ds_nm_and_url_chr[1], 
            results_ls$path_params_ls$dv_ds_nm_and_url_chr[1]), 
            ds_url_1L_chr = ifelse(!is.null(input_params_ls), 
                input_params_ls$path_params_ls$dv_ds_nm_and_url_chr[2], 
                results_ls$path_params_ls$dv_ds_nm_and_url_chr[2]), 
            parent_dv_dir_1L_chr = paths_ls$dv_dir_1L_chr, paths_to_dirs_chr = paths_ls$reports_dir_1L_chr, 
            inc_fl_types_chr = ifelse(output_type_1L_chr == "Word", 
                ".docx", ".pdf"), paths_are_rltv_1L_lgl = F)
    }
    results_ls$path_params_ls$paths_ls$path_to_ms_mkdn_1L_dir <- path_to_ms_mkdn_1L_dir
    saveRDS(results_ls, paste0(outp_dir_1L_chr, "/results_ls.RDS"))
    return(results_ls)
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
write_mdl_cmprsn <- function (scored_data_tb, ds_smry_ls, mdl_smry_ls, output_data_dir_1L_chr, 
    seed_1L_int = 1234) 
{
    bl_tb <- youthvars::transform_ds_for_tstng(scored_data_tb, 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr, 
        depnt_var_max_val_1L_dbl = 0.999, round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr, 
        round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr)
    ds_smry_ls$candidate_predrs_chr <- reorder_cndt_predrs_chr(ds_smry_ls$candidate_predrs_chr, 
        data_tb = bl_tb, depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr)
    mdl_smry_ls <- add_prefd_predr_var_to_mdl_smry_ls(mdl_smry_ls, 
        ds_smry_ls = ds_smry_ls)
    mdl_smry_ls$smry_of_sngl_predr_mdls_tb <- write_sngl_predr_multi_mdls_outps(data_tb = bl_tb, 
        folds_1L_int = mdl_smry_ls$folds_1L_int, mdl_types_chr = mdl_smry_ls$mdl_types_chr, 
        depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        predr_var_nm_1L_chr = mdl_smry_ls$predr_var_nm_1L_chr, 
        predr_var_desc_1L_chr = mdl_smry_ls$predr_var_desc_1L_chr, 
        predr_vals_dbl = mdl_smry_ls$predr_vals_dbl, path_to_write_to_1L_chr = output_data_dir_1L_chr, 
        new_dir_nm_1L_chr = "A_Candidate_Mdls_Cmprsn", mdl_types_lup = mdl_smry_ls$mdl_types_lup, 
        dictionary_tb = ds_smry_ls$dictionary_tb)
    mdl_smry_ls$prefd_mdl_types_chr <- make_prefd_mdls_vec(mdl_smry_ls$smry_of_sngl_predr_mdls_tb, 
        choose_from_pfx_chr = mdl_smry_ls$choose_from_pfx_chr)
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
#' @param depnt_var_desc_1L_chr Dependent variable description (a character vector of length one), Default: 'Utility score'
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
    depnt_var_desc_1L_chr = "Utility score", tfmn_1L_chr = "NTF", 
    predr_var_nm_1L_chr, predr_var_desc_1L_chr, predr_vals_dbl, 
    covar_var_nms_chr = NA_character_, path_to_write_to_1L_chr, 
    predn_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_, 
    plt_idxs_int = 1:5) 
{
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    tfd_data_tb <- transform_data_tb_for_cmprsn(data_tb, model_mdl = model_mdl, 
        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, predn_type_1L_chr = predn_type_1L_chr, 
        tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, family_1L_chr = family_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)
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
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, predr_var_desc_1L_chr = predr_var_desc_1L_chr, 
            depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, depnt_var_desc_1L_chr = depnt_var_desc_1L_chr), 
            list(model_mdl, which_dbl = 1:6, ncol_1L_int = 3L, 
                label_size_1L_int = 3), list(tfd_data_tb = tfd_data_tb, 
                depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, depnt_var_desc_1L_chr = depnt_var_desc_1L_chr), 
            list(tfd_data_tb = transform_data_tb_for_cmprsn(data_tb, 
                model_mdl = model_mdl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                new_data_is_1L_chr = ifelse(!4 %in% plt_idxs_int, 
                  "Predicted", "Simulated"), predn_type_1L_chr = NULL, 
                tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, 
                family_1L_chr = family_1L_chr, tfmn_1L_chr = tfmn_1L_chr), 
                depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, depnt_var_desc_1L_chr = depnt_var_desc_1L_chr, 
                predd_val_var_nm_1L_chr = "Simulated"), list(tfd_data_tb = tfd_data_tb, 
                depnt_var_nm_1L_chr = depnt_var_nm_1L_chr))[plt_idxs_int], 
        plt_nm_sfx_chr = c("_LNR_CMPRSN", "_AUTOPLT", "_PRED_DNSTY", 
            "_SIM_DNSTY", "_PRED_SCTR")[plt_idxs_int], size_ls = list(c(6, 
            6), c(4, 7), c(6, 6), c(6, 6), c(6, 6))[plt_idxs_int]), 
        ~ready4show::write_mdl_plt_fl(plt_fn = ..1, fn_args_ls = ..2, 
            path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
            plt_nm_1L_chr = paste0(mdl_fl_nm_1L_chr, ifelse(!is.na(covar_var_nms_chr[1]), 
                paste("_", paste0(covar_var_nms_chr[1:min(length(covar_var_nms_chr), 
                  3)], collapse = "")), ""), ..3), height_1L_dbl = ..4[1], 
            width_1L_dbl = ..4[2]))
}
#' Write model summary report
#' @description write_mdl_smry_rprt() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write model summary report. The function returns Input parameters (a list).
#' @param input_params_ls Input parameters (a list), Default: NULL
#' @param header_yaml_args_ls Header yaml arguments (a list), Default: NULL
#' @param path_params_ls Path parameters (a list), Default: NULL
#' @param use_fake_data_1L_lgl Use fake data (a logical vector of length one), Default: FALSE
#' @param output_format_ls Output format (a list), Default: NULL
#' @param abstract_args_ls Abstract arguments (a list), Default: NULL
#' @param dv_ds_nm_and_url_chr Dataverse dataset name and url (a character vector), Default: NULL
#' @param reference_int Reference (an integer vector), Default: 0
#' @param rprt_lup Report (a lookup table), Default: NULL
#' @param rcrd_nm_1L_chr Record name (a character vector of length one), Default: 'AAA_RPRT_WRTNG_MTH'
#' @param rprt_nm_1L_chr Report name (a character vector of length one), Default: 'AAA_TTU_MDL_CTG'
#' @param start_at_int Start at (an integer vector), Default: c(2, 1)
#' @param use_shareable_mdls_1L_lgl Use shareable models (a logical vector of length one), Default: F
#' @return Input parameters (a list)
#' @rdname write_mdl_smry_rprt
#' @export 
#' @importFrom purrr map pluck map_lgl map_chr reduce
#' @importFrom here here
#' @importFrom ready4 write_to_dv_with_wait
#' @importFrom tibble tibble
#' @importFrom dplyr filter pull bind_rows distinct mutate
#' @importFrom stats setNames
#' @keywords internal
write_mdl_smry_rprt <- function (input_params_ls = NULL, header_yaml_args_ls = NULL, 
    path_params_ls = NULL, use_fake_data_1L_lgl = FALSE, output_format_ls = NULL, 
    abstract_args_ls = NULL, dv_ds_nm_and_url_chr = NULL, reference_int = 0, 
    rprt_lup = NULL, rcrd_nm_1L_chr = "AAA_RPRT_WRTNG_MTH", rprt_nm_1L_chr = "AAA_TTU_MDL_CTG", 
    start_at_int = c(2, 1), use_shareable_mdls_1L_lgl = F) 
{
    if (missing(header_yaml_args_ls)) {
        header_yaml_args_ls <- input_params_ls$header_yaml_args_ls
    }
    else {
        warning("The argument header_yaml_args_ls is soft deprecated. We recommend passing the header information as part of the list passed to the input_params_ls argument.")
    }
    if (missing(path_params_ls)) {
        path_params_ls <- input_params_ls$path_params_ls
    }
    else {
        warning("The argument path_params_ls is soft deprecated. We recommend passing the paths information as part of the list passed to the input_params_ls argument.")
    }
    if (missing(output_format_ls)) {
        output_format_ls <- input_params_ls$output_format_ls
    }
    else {
        warning("The argument output_format_ls is soft deprecated. We recommend passing the output format information as part of the list passed to the input_params_ls argument.")
    }
    if (missing(use_fake_data_1L_lgl)) {
        use_fake_data_1L_lgl <- input_params_ls$params_ls$use_fake_data_1L_lgl
    }
    else {
        warning("The argument use_fake_data_1L_lgl is soft deprecated. We recommend declaring whether dataset is fake as part of the list passed to the input_params_ls argument.")
    }
    if (missing(reference_int)) {
        reference_int <- 0:(ifelse(is.null(input_params_ls$scndry_anlys_params_ls), 
            0, length(input_params_ls$scndry_anlys_params_ls)))
    }
    else {
        warning("The argument reference_int is soft deprecated. It is unnecessary if supplying a valid value to the input_params_ls argument.")
    }
    paths_ls <- path_params_ls$paths_ls
    if (is.null(rprt_lup)) 
        data("rprt_lup", package = "TTU", envir = environment())
    rprt_lups_ls <- purrr::map(reference_int, ~{
        if (.x == 0) {
            reference_1L_int <- NULL
        }
        else {
            reference_1L_int <- .x
        }
        rprt_lup <- rprt_lup %>% transform_rprt_lup(add_suplry_rprt_1L_lgl = !is.null(reference_1L_int), 
            add_sharing_rprt_1L_lgl = T, start_at_int = start_at_int, 
            reference_1L_int = reference_1L_int)
        if (is.null(reference_1L_int)) {
            path_to_outp_fl_1L_chr <- paste0(paths_ls$output_data_dir_1L_chr, 
                "/I_ALL_OUTPUT_.RDS")
            if (use_shareable_mdls_1L_lgl) {
                main_rprt_append_ls <- list(rltv_path_to_data_dir_1L_chr = "../Output/G_Shareable/Models")
            }
            else {
                main_rprt_append_ls <- NULL
            }
            rcrd_rprt_append_ls <- path_params_ls[1:2]
        }
        else {
            path_to_outp_fl_1L_chr <- here::here(paths_ls$path_from_top_level_1L_chr, 
                paths_ls$write_to_dir_nm_1L_chr, paste0("secondary_", 
                  reference_1L_int), "Output", "I_ALL_OUTPUT_.RDS")
            main_rprt_append_ls <- list(existing_predrs_ls = readRDS(paste0(paths_ls$output_data_dir_1L_chr, 
                "/I_ALL_OUTPUT_.RDS")) %>% purrr::pluck("predr_vars_nms_ls"))
            if (use_shareable_mdls_1L_lgl) {
                main_rprt_append_ls$rltv_path_to_data_dir_1L_chr <- "../Output/G_Shareable/Models"
            }
            else {
                main_rprt_append_ls$rltv_path_to_data_dir_1L_chr <- NULL
            }
            paths_ls <- transform_paths_ls_for_scndry(paths_ls, 
                reference_1L_int = reference_1L_int, remove_prmry_1L_lgl = T)
            rcrd_rprt_append_ls <- list(transform_paths_ls = list(fn = transform_paths_ls_for_scndry, 
                args_ls = list(reference_1L_int = reference_1L_int, 
                  remove_prmry_1L_lgl = T))) %>% append(path_params_ls[1:2])
        }
        write_rprt_with_rcrd(path_to_outp_fl_1L_chr = path_to_outp_fl_1L_chr, 
            paths_ls = paths_ls, header_yaml_args_ls = header_yaml_args_ls, 
            use_fake_data_1L_lgl = use_fake_data_1L_lgl, reference_1L_int = reference_1L_int, 
            start_at_int = start_at_int, rprt_nm_1L_chr = rprt_lup$rprt_nms_chr[purrr::map_lgl(rprt_lup$rprt_nms_chr, 
                ~startsWith(.x, rprt_nm_1L_chr))], rcrd_nm_1L_chr = rcrd_nm_1L_chr, 
            output_type_1L_chr = output_format_ls$supplementary_outp_1L_chr, 
            rprt_output_type_1L_chr = output_format_ls$supplementary_outp_1L_chr, 
            nbr_of_digits_1L_int = output_format_ls$supplementary_digits_1L_int, 
            abstract_args_ls = abstract_args_ls, rcrd_rprt_append_ls = rcrd_rprt_append_ls, 
            rprt_lup = rprt_lup, main_rprt_append_ls = main_rprt_append_ls)
        if (!is.null(dv_ds_nm_and_url_chr)) {
            ready4::write_to_dv_with_wait(dss_tb = tibble::tibble(ds_obj_nm_chr = c(rprt_nm_1L_chr, 
                rcrd_nm_1L_chr), title_chr = rprt_lup %>% dplyr::filter(rprt_nms_chr %in% 
                c(rprt_nm_1L_chr, rcrd_nm_1L_chr)) %>% dplyr::pull(title_chr)), 
                dv_nm_1L_chr = dv_ds_nm_and_url_chr[1], ds_url_1L_chr = dv_ds_nm_and_url_chr[2], 
                parent_dv_dir_1L_chr = paths_ls$dv_dir_1L_chr, 
                paths_to_dirs_chr = paths_ls$reports_dir_1L_chr, 
                inc_fl_types_chr = ".pdf", paths_are_rltv_1L_lgl = F)
        }
        rprt_lup
    }) %>% stats::setNames(reference_int %>% purrr::map_chr(~ifelse(.x == 
        0, "Primary", paste0("secondary_", .x))))
    consolidated_mdl_ings_ls <- reference_int %>% purrr::reduce(.init = paste0(paths_ls$output_data_dir_1L_chr, 
        "/G_Shareable/Ingredients/mdl_ingredients.RDS") %>% readRDS(), 
        ~if (.y > 0) {
            ingredients_ls <- here::here(paths_ls$path_from_top_level_1L_chr, 
                paths_ls$write_to_dir_nm_1L_chr, paste0("secondary_", 
                  .y), "Output", "G_Shareable", "Ingredients", 
                "mdl_ingredients.RDS") %>% readRDS()
            .x <- append(.x, list(ingredients_ls) %>% setNames(paste0("secondary_", 
                .y)))
            .x$dictionary_tb <- dplyr::bind_rows(.x$dictionary_tb, 
                ingredients_ls$dictionary_tb) %>% dplyr::distinct()
            .x$mdls_lup <- dplyr::bind_rows(.x$mdls_lup, ingredients_ls$mdls_lup %>% 
                dplyr::mutate(source_chr = paste0("Secondary Analysis ", 
                  LETTERS[.y]))) %>% dplyr::distinct()
            .x$mdls_smry_tb <- dplyr::bind_rows(.x$mdls_smry_tb, 
                ingredients_ls$mdls_smry_tb) %>% dplyr::distinct()
            .x$predictors_lup <- dplyr::bind_rows(.x$predictors_lup, 
                ingredients_ls$predictors_lup) %>% dplyr::distinct()
            .x
        }
        else {
            .x$Primary <- .x
            .x$mdls_lup <- .x$mdls_lup %>% dplyr::mutate(source_chr = "Primary Analysis")
            .x
        })
    saveRDS(consolidated_mdl_ings_ls, paste0(paths_ls$output_data_dir_1L_chr, 
        "/G_Shareable/Ingredients/mdl_ingredients.RDS"))
    if (!is.null(input_params_ls)) {
        input_params_ls$rprt_lups_ls <- rprt_lups_ls
    }
    else {
        input_params_ls <- rprt_lups_ls
    }
    return(input_params_ls)
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
#' @importFrom purrr map_chr map_dfr
#' @importFrom ready4 get_from_lup_obj
#' @importFrom tibble tibble
#' @importFrom caret R2
#' @importFrom stats predict AIC BIC
#' @importFrom dplyr pull arrange desc
#' @importFrom rlang sym
write_mdl_type_covars_mdls <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predrs_var_nms_chr, 
    covar_var_nms_chr, mdl_type_1L_chr, path_to_write_to_1L_chr, 
    new_dir_nm_1L_chr = "D_Covars_Selection", fl_nm_pfx_1L_chr = "D_CT", 
    mdl_types_lup = NULL, start_1L_chr = NA_character_) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    arg_vals_chr <- c("control_chr", "predn_type_chr", "tfmn_chr") %>% 
        purrr::map_chr(~ready4::get_from_lup_obj(mdl_types_lup, 
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
            target_var_nm_1L_chr = .x, evaluate_1L_lgl = F))
    control_1L_chr <- arg_vals_chr[1]
    predn_type_1L_chr <- arg_vals_chr[2]
    tfmn_1L_chr <- arg_vals_chr[3]
    if (is.na(predn_type_1L_chr)) 
        predn_type_1L_chr <- NULL
    data_tb <- data_tb %>% add_tfd_var_to_ds(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)
    output_dir_1L_chr <- output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    smry_of_mdls_with_covars_tb <- purrr::map_dfr(predrs_var_nms_chr, 
        ~{
            model_mdl <- make_mdl(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                predr_var_nm_1L_chr = .x, covar_var_nms_chr = covar_var_nms_chr, 
                tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr, 
                control_1L_chr = control_1L_chr, mdl_types_lup = mdl_types_lup, 
                start_1L_chr = start_1L_chr)
            mdl_fl_nm_1L_chr <- paste0(fl_nm_pfx_1L_chr, "_", 
                .x, "_", mdl_type_1L_chr)
            saveRDS(model_mdl, paste0(output_dir_1L_chr, "/", 
                mdl_fl_nm_1L_chr, ".RDS"))
            if ("summary.betareg" %in% class(summary(model_mdl))) {
                coefficients_mat <- summary(model_mdl)$coefficients$mean
            }
            else {
                coefficients_mat <- summary(model_mdl)$coefficients
            }
            tibble::tibble(variable = .x, Rsquare = caret::R2(stats::predict(model_mdl, 
                type = predn_type_1L_chr) %>% calculate_depnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, 
                tfmn_is_outp_1L_lgl = T), data_tb %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)), 
                form = "traditional"), AIC = stats::AIC(model_mdl), 
                BIC = stats::BIC(model_mdl), Significant = paste(names(which(coefficients_mat[, 
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
#' @importFrom ready4 get_from_lup_obj
#' @importFrom dplyr select mutate everything arrange desc
write_mdl_type_multi_outps <- function (data_tb, folds_1L_int = 10, predrs_var_nms_chr, covar_var_nms_chr = NA_character_, 
    start_1L_chr = NULL, mdl_type_1L_chr, depnt_var_nm_1L_chr = "utl_total_w", 
    path_to_write_to_1L_chr, new_dir_nm_1L_chr, mdl_types_lup = NULL, 
    fl_nm_pfx_1L_chr = "C_PREDR", plt_idxs_int = c(3, 5)) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    smry_of_mdl_sngl_predrs_tb <- purrr::map_dfr(predrs_var_nms_chr, 
        ~{
            tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
                target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F)
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
#' @importFrom ready4 get_from_lup_obj
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
        purrr::map_chr(~ready4::get_from_lup_obj(mdl_types_lup, 
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
            target_var_nm_1L_chr = .x, evaluate_1L_lgl = F))
    control_1L_chr <- arg_vals_chr[1]
    family_1L_chr <- arg_vals_chr[2]
    predn_type_1L_chr <- arg_vals_chr[3]
    if (is.na(predn_type_1L_chr)) 
        predn_type_1L_chr <- NULL
    if (is.na(plt_idxs_int[1])) {
        plt_idxs_int <- 1:5
        if (!is.na(control_1L_chr)) {
            if (control_1L_chr %>% startsWith("betareg")) 
                plt_idxs_int <- c(1, 3, 4, 5)
        }
    }
    tfmn_for_bnml_1L_lgl <- ready4::get_from_lup_obj(mdl_types_lup, 
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
        target_var_nm_1L_chr = "tfmn_for_bnml_lgl", evaluate_1L_lgl = F)
    data_tb <- data_tb %>% add_tfd_var_to_ds(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)
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
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr, 
            mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup, 
            start_1L_chr = start_1L_chr, predn_type_1L_chr = predn_type_1L_chr)
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
#' @param session_data_ls Session data (a list), Default: NULL
#' @return Output summary (a list)
#' @rdname write_mdls_with_covars_cmprsn
#' @export 
write_mdls_with_covars_cmprsn <- function (scored_data_tb, bl_tb, ds_smry_ls, mdl_smry_ls, output_data_dir_1L_chr, 
    seed_1L_int = 1234, session_data_ls = NULL) 
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
    outp_smry_ls <- list(scored_data_tb = scored_data_tb, dictionary_tb = ds_smry_ls$dictionary_tb, 
        predictors_lup = ds_smry_ls$predictors_lup, smry_of_sngl_predr_mdls_tb = mdl_smry_ls$smry_of_sngl_predr_mdls_tb, 
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
        fl_nm_pfx_1L_chr = "C_PREDR", start_1L_chr = NA_character_, 
        mdl_types_lup = mdl_smry_ls$mdl_types_lup)
    bl_tb <- scored_data_tb %>% youthvars::transform_ds_for_tstng(depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr, 
        candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr, 
        covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr, 
        remove_all_msng_1L_lgl = T, round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr, 
        round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr)
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
#' Write report
#' @description write_report() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write report. The function is called for its side effects and does not return a value. WARNING: This function writes R scripts to your local environment. Make sure to only use if you want this behaviour
#' @param params_ls Parameters (a list)
#' @param paths_ls Paths (a list)
#' @param rprt_nm_1L_chr Report name (a character vector of length one)
#' @param abstract_args_ls Abstract arguments (a list), Default: NULL
#' @param header_yaml_args_ls Header yaml arguments (a list), Default: NULL
#' @param rprt_lup Report (a lookup table), Default: NULL
#' @return NULL
#' @rdname write_report
#' @export 
#' @importFrom ready4show make_rprt_type_ls write_rprt_from_tmpl
#' @importFrom here i_am here
#' @importFrom rlang exec
write_report <- function (params_ls, paths_ls, rprt_nm_1L_chr, abstract_args_ls = NULL, 
    header_yaml_args_ls = NULL, rprt_lup = NULL) 
{
    if (is.null(rprt_lup)) 
        data("rprt_lup", package = "TTU", envir = environment())
    rprt_type_ls <- rprt_lup %>% ready4show::make_rprt_type_ls(rprt_nm_1L_chr = rprt_nm_1L_chr)
    here::i_am(paste0(paths_ls$path_from_top_level_1L_chr, "/", 
        paths_ls$path_to_current_1L_chr, "/", paths_ls$R_fl_nm_1L_chr))
    args_ls <- list(rprt_type_ls = rprt_type_ls, params_ls = params_ls, 
        output_type_1L_chr = params_ls$output_type_1L_chr, path_to_prjs_dir_1L_chr = here::here(paths_ls$path_from_top_level_1L_chr), 
        prj_dir_1L_chr = paths_ls$write_to_dir_nm_1L_chr, header_yaml_args_ls = header_yaml_args_ls, 
        abstract_args_ls = abstract_args_ls, reports_dir_1L_chr = "Reports", 
        rltv_path_to_data_dir_1L_chr = "../Output", nm_of_mkdn_dir_1L_chr = "Markdown")
    rlang::exec(ready4show::write_rprt_from_tmpl, !!!args_ls)
}
#' Write reporting directory
#' @description write_reporting_dir() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write reporting directory. The function returns Path to complete study program (a character vector of length one).
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one), Default: getwd()
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'TTU_Project'
#' @param overwrite_1L_lgl Overwrite (a logical vector of length one), Default: FALSE
#' @return Path to complete study program (a character vector of length one)
#' @rdname write_reporting_dir
#' @export 
write_reporting_dir <- function (path_to_write_to_1L_chr = getwd(), new_dir_nm_1L_chr = "TTU_Project", 
    overwrite_1L_lgl = FALSE) 
{
    path_to_prjt_dir_1L_chr <- paste0(path_to_write_to_1L_chr, 
        "/", new_dir_nm_1L_chr)
    if (!dir.exists(path_to_prjt_dir_1L_chr)) 
        dir.create(path_to_prjt_dir_1L_chr)
    path_to_RMD_dir_1L_chr <- system.file("Project/CSP", package = "TTU")
    file.copy(path_to_RMD_dir_1L_chr, path_to_prjt_dir_1L_chr, 
        recursive = T, overwrite = overwrite_1L_lgl)
    path_to_csp_1L_chr <- paste0(path_to_prjt_dir_1L_chr, "/CSP/CSP.Rmd")
    return(path_to_csp_1L_chr)
}
#' Write report with record
#' @description write_rprt_with_rcrd() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write report with record. The function is called for its side effects and does not return a value. WARNING: This function writes R scripts to your local environment. Make sure to only use if you want this behaviour
#' @param path_to_outp_fl_1L_chr Path to output file (a character vector of length one)
#' @param paths_ls Paths (a list)
#' @param header_yaml_args_ls Header yaml arguments (a list), Default: NULL
#' @param rprt_lup Report (a lookup table), Default: NULL
#' @param use_fake_data_1L_lgl Use fake data (a logical vector of length one), Default: F
#' @param rprt_nm_1L_chr Report name (a character vector of length one), Default: 'AAA_TTU_MDL_CTG'
#' @param rcrd_nm_1L_chr Record name (a character vector of length one), Default: 'AAA_RPRT_WRTNG_MTH'
#' @param reference_1L_int Reference (an integer vector of length one), Default: NULL
#' @param start_at_int Start at (an integer vector), Default: c(2, 1)
#' @param output_type_1L_chr Output type (a character vector of length one), Default: 'PDF'
#' @param rprt_output_type_1L_chr Report output type (a character vector of length one), Default: 'PDF'
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @param abstract_args_ls Abstract arguments (a list), Default: NULL
#' @param main_rprt_append_ls Main report append (a list), Default: NULL
#' @param rcrd_rprt_append_ls Record report append (a list), Default: NULL
#' @return NULL
#' @rdname write_rprt_with_rcrd
#' @export 
#' @importFrom ready4 get_from_lup_obj
#' @keywords internal
write_rprt_with_rcrd <- function (path_to_outp_fl_1L_chr, paths_ls, header_yaml_args_ls = NULL, 
    rprt_lup = NULL, use_fake_data_1L_lgl = F, rprt_nm_1L_chr = "AAA_TTU_MDL_CTG", 
    rcrd_nm_1L_chr = "AAA_RPRT_WRTNG_MTH", reference_1L_int = NULL, 
    start_at_int = c(2, 1), output_type_1L_chr = "PDF", rprt_output_type_1L_chr = "PDF", 
    nbr_of_digits_1L_int = 2L, abstract_args_ls = NULL, main_rprt_append_ls = NULL, 
    rcrd_rprt_append_ls = NULL) 
{
    if (is.null(rprt_lup)) {
        data("rprt_lup", package = "TTU", envir = environment())
        rprt_lup <- rprt_lup %>% transform_rprt_lup(add_suplry_rprt_1L_lgl = !is.null(reference_1L_int), 
            add_sharing_rprt_1L_lgl = F, start_at_int = start_at_int, 
            reference_1L_int = reference_1L_int)
    }
    params_ls <- list(abstract_args_ls = NULL, eval_1L_lgl = F, 
        header_yaml_args_ls = header_yaml_args_ls, output_type_1L_chr = rprt_output_type_1L_chr, 
        nbr_of_digits_1L_int = nbr_of_digits_1L_int, rprt_lup = rprt_lup, 
        rprt_nm_1L_chr = rprt_nm_1L_chr, rprt_output_type_1L_chr = output_type_1L_chr, 
        rprt_subtitle_1L_chr = ready4::get_from_lup_obj(rprt_lup, 
            match_value_xx = rprt_nm_1L_chr, match_var_nm_1L_chr = "rprt_nms_chr", 
            target_var_nm_1L_chr = "title_chr", evaluate_1L_lgl = F), 
        subtitle_1L_chr = ready4::get_from_lup_obj(rprt_lup, 
            match_value_xx = "AAA_RPRT_WRTNG_MTH", match_var_nm_1L_chr = "rprt_nms_chr", 
            target_var_nm_1L_chr = "title_chr", evaluate_1L_lgl = F), 
        use_fake_data_1L_lgl = use_fake_data_1L_lgl) %>% append(rcrd_rprt_append_ls)
    params_ls %>% write_report(paths_ls = paths_ls, rprt_nm_1L_chr = rcrd_nm_1L_chr, 
        abstract_args_ls = NULL, header_yaml_args_ls = header_yaml_args_ls, 
        rprt_lup = rprt_lup)
    list(outp_smry_ls = append(readRDS(path_to_outp_fl_1L_chr), 
        list(rprt_lup = rprt_lup)), output_type_1L_chr = output_type_1L_chr, 
        subtitle_1L_chr = ready4::get_from_lup_obj(rprt_lup, 
            match_value_xx = rprt_nm_1L_chr, match_var_nm_1L_chr = "rprt_nms_chr", 
            target_var_nm_1L_chr = "title_chr", evaluate_1L_lgl = F)) %>% 
        append(main_rprt_append_ls) %>% write_report(paths_ls = paths_ls, 
        rprt_nm_1L_chr = rprt_nm_1L_chr, abstract_args_ls = abstract_args_ls, 
        header_yaml_args_ls = header_yaml_args_ls, rprt_lup = rprt_lup)
}
#' Write secondary analysis
#' @description write_scndry_analysis() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write secondary analysis. The function is called for its side effects and does not return a value. WARNING: This function writes R scripts to your local environment. Make sure to only use if you want this behaviour
#' @param predictors_lup Predictors (a lookup table), Default: NULL
#' @param valid_params_ls_ls Valid parameters (a list of lists)
#' @param candidate_covar_nms_chr Candidate covariate names (a character vector)
#' @param candidate_predrs_chr Candidate predictors (a character vector), Default: NULL
#' @param header_yaml_args_ls Header yaml arguments (a list)
#' @param path_params_ls Path parameters (a list)
#' @param prefd_covars_chr Preferred covariates (a character vector), Default: 'NA'
#' @param reference_1L_int Reference (an integer vector of length one)
#' @param start_at_int Start at (an integer vector), Default: c(2, 1)
#' @param rprt_nm_1L_chr Report name (a character vector of length one), Default: 'AAA_SUPLRY_ANLYS_MTH'
#' @param abstract_args_ls Abstract arguments (a list), Default: NULL
#' @return NULL
#' @rdname write_scndry_analysis
#' @export 
#' @importFrom purrr map_chr pluck
#' @importFrom ready4 get_from_lup_obj
#' @importFrom dplyr filter
#' @keywords internal
write_scndry_analysis <- function (predictors_lup = NULL, valid_params_ls_ls, candidate_covar_nms_chr, 
    candidate_predrs_chr = NULL, header_yaml_args_ls, path_params_ls, 
    prefd_covars_chr = NA_character_, reference_1L_int, start_at_int = c(2, 
        1), rprt_nm_1L_chr = "AAA_SUPLRY_ANLYS_MTH", abstract_args_ls = NULL) 
{
    analysis_params_ls <- valid_params_ls_ls$params_ls %>% append(path_params_ls[1:2])
    rename_lup <- valid_params_ls_ls$rename_lup
    if (!is.null(predictors_lup)) {
        predictors_lup$short_name_chr <- predictors_lup$short_name_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_1L_lgl = F)))
        analysis_params_ls$predictors_lup <- predictors_lup
    }
    if (!is.null(candidate_predrs_chr)) {
        candidate_predrs_chr <- candidate_predrs_chr %>% purrr::map_chr(~ifelse(!.x %in% 
            rename_lup$old_nms_chr, .x, ready4::get_from_lup_obj(rename_lup, 
            match_value_xx = .x, match_var_nm_1L_chr = "old_nms_chr", 
            target_var_nm_1L_chr = "new_nms_chr", evaluate_1L_lgl = F)))
        analysis_params_ls$ds_descvs_ls$candidate_predrs_chr <- candidate_predrs_chr
    }
    if (!is.null(candidate_covar_nms_chr)) {
        candidate_covar_nms_chr <- candidate_covar_nms_chr %>% 
            purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr, 
                .x, ready4::get_from_lup_obj(rename_lup, match_value_xx = .x, 
                  match_var_nm_1L_chr = "old_nms_chr", target_var_nm_1L_chr = "new_nms_chr", 
                  evaluate_1L_lgl = F)))
    }
    if (ifelse(is.null(prefd_covars_chr), F, !is.na(prefd_covars_chr))) {
        prefd_covars_chr <- prefd_covars_chr %>% purrr::map_chr(~ifelse(!.x %in% 
            rename_lup$old_nms_chr, .x, ready4::get_from_lup_obj(rename_lup, 
            match_value_xx = .x, match_var_nm_1L_chr = "old_nms_chr", 
            target_var_nm_1L_chr = "new_nms_chr", evaluate_1L_lgl = F)))
    }
    data("rprt_lup", package = "TTU", envir = environment())
    rprt_lup <- rprt_lup %>% transform_rprt_lup(start_at_int = start_at_int, 
        reference_1L_int = reference_1L_int) %>% dplyr::filter(rprt_nms_chr == 
        "AAA_SUPLRY_ANLYS_MTH")
    analysis_params_ls$subtitle_1L_chr <- ready4::get_from_lup_obj(rprt_lup, 
        match_value_xx = "AAA_SUPLRY_ANLYS_MTH", match_var_nm_1L_chr = "rprt_nms_chr", 
        target_var_nm_1L_chr = "title_chr", evaluate_1L_lgl = F)
    analysis_params_ls$prefd_covars_chr <- prefd_covars_chr
    analysis_params_ls$candidate_covar_nms_chr <- candidate_covar_nms_chr
    path_params_ls$paths_ls <- write_scndry_analysis_dir(path_params_ls$paths_ls, 
        reference_1L_int = reference_1L_int)
    params_ls <- list(candidate_predrs_chr = candidate_predrs_chr, 
        transform_paths_ls = list(fn = transform_paths_ls_for_scndry, 
            args_ls = list(reference_1L_int = reference_1L_int))) %>% 
        append(analysis_params_ls)
    params_ls$utl_class_fn_1L_chr <- params_ls$raw_ds_tfmn_fn <- NULL
    params_ls_ls <- transform_params_ls_to_valid(params_ls)
    params_ls <- params_ls_ls %>% purrr::pluck("params_ls") %>% 
        append(list(rename_lup = params_ls_ls$rename_lup))
    params_ls %>% write_report(paths_ls = path_params_ls$paths_ls, 
        rprt_nm_1L_chr = rprt_nm_1L_chr, abstract_args_ls = abstract_args_ls, 
        header_yaml_args_ls = header_yaml_args_ls, rprt_lup = rprt_lup)
}
#' Write secondary analysis directory
#' @description write_scndry_analysis_dir() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write secondary analysis directory. The function returns Paths (a list).
#' @param paths_ls Paths (a list)
#' @param reference_1L_int Reference (an integer vector of length one), Default: 1
#' @return Paths (a list)
#' @rdname write_scndry_analysis_dir
#' @export 
#' @importFrom here here
#' @keywords internal
write_scndry_analysis_dir <- function (paths_ls, reference_1L_int = 1) 
{
    paths_ls <- transform_paths_ls_for_scndry(paths_ls, reference_1L_int = reference_1L_int)
    paste0(here::here(paths_ls$path_from_top_level_1L_chr), "/", 
        paths_ls$write_to_dir_nm_1L_chr) %>% dir.create()
    return(paths_ls)
}
#' Write shareable directory
#' @description write_shareable_dir() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write shareable directory. The function returns Output directory (a character vector).
#' @param outp_smry_ls Output summary (a list)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'G_Shareable'
#' @param sub_dirs_chr Sub directories (a character vector), Default: c("Ingredients", "Models", "Table_Predn_Tools")
#' @return Output directory (a character vector)
#' @rdname write_shareable_dir
#' @export 
#' @importFrom purrr map_chr
#' @keywords internal
write_shareable_dir <- function (outp_smry_ls, new_dir_nm_1L_chr = "G_Shareable", sub_dirs_chr = c("Ingredients", 
    "Models", "Table_Predn_Tools")) 
{
    output_dir_chr <- write_new_outp_dir(outp_smry_ls$path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    output_dir_chr <- c(output_dir_chr, sub_dirs_chr %>% purrr::map_chr(~write_new_outp_dir(output_dir_chr, 
        new_dir_nm_1L_chr = .x)))
    return(output_dir_chr)
}
#' Write shareable models
#' @description write_shareable_mdls() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write shareable models. The function returns Output summary (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'G_Shareable'
#' @param shareable_title_detail_1L_chr Shareable title detail (a character vector of length one), Default: ''
#' @param write_mdls_to_dv_1L_lgl Write models to dataverse (a logical vector of length one), Default: F
#' @return Output summary (a list)
#' @rdname write_shareable_mdls
#' @export 
#' @importFrom purrr map_chr flatten_chr map map_lgl map_int map2 map_dbl map2_dfr discard
#' @importFrom stringr str_locate str_remove_all
#' @importFrom dplyr filter select mutate
#' @importFrom ready4 get_from_lup_obj
#' @importFrom stats setNames
#' @importFrom tibble tibble
write_shareable_mdls <- function (outp_smry_ls, new_dir_nm_1L_chr = "G_Shareable", shareable_title_detail_1L_chr = "", 
    write_mdls_to_dv_1L_lgl = F) 
{
    output_dir_chr <- write_shareable_dir(outp_smry_ls = outp_smry_ls, 
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
    fake_ds_tb <- make_fake_ts_data(outp_smry_ls, depnt_vars_are_NA_1L_lgl = F)
    mdl_types_lup <- outp_smry_ls$mdl_types_lup
    shareable_mdls_ls <- outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr() %>% 
        purrr::map2(incld_mdl_paths_chr, ~{
            model_mdl <- readRDS(paste0(outp_smry_ls$path_to_write_to_1L_chr, 
                "/", .y))
            mdl_smry_tb <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Model == 
                .x)
            mdl_nm_1L_chr <- .x
            mdl_type_1L_chr <- get_mdl_type_from_nm(mdl_nm_1L_chr, 
                mdl_types_lup = mdl_types_lup)
            tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F)
            predn_type_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                target_var_nm_1L_chr = "predn_type_chr", evaluate_1L_lgl = F)
            if (is.na(predn_type_1L_chr)) 
                predn_type_1L_chr <- NULL
            control_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                target_var_nm_1L_chr = "control_chr", evaluate_1L_lgl = F)
            sd_dbl <- mdl_smry_tb %>% dplyr::filter(Parameter == 
                "SD (Intercept)") %>% dplyr::select(Estimate, 
                SE) %>% t() %>% as.vector()
            mdl_fake_ds_tb <- fake_ds_tb %>% add_tfd_var_to_ds(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                tfmn_1L_chr = tfmn_1L_chr, depnt_var_max_val_1L_dbl = 0.999) %>% 
                dplyr::select(names(model_mdl$data))
            model_mdl$data <- mdl_fake_ds_tb
            table_predn_mdl <- make_shareable_mdl(fake_ds_tb = mdl_fake_ds_tb, 
                mdl_smry_tb = mdl_smry_tb, depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
                tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr, 
                mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr, 
                start_1L_chr = NA_character_, seed_1L_int = outp_smry_ls$seed_1L_int)
            saveRDS(table_predn_mdl, paste0(output_dir_chr[4], 
                "/", .x, ".RDS"))
            saveRDS(model_mdl, paste0(output_dir_chr[3], "/", 
                .x, ".RDS"))
            scaling_fctr_dbl <- outp_smry_ls$predr_vars_nms_ls %>% 
                purrr::flatten_chr() %>% unique() %>% purrr::map_dbl(~ifelse(.x %in% 
                outp_smry_ls$predictors_lup$short_name_chr, ready4::get_from_lup_obj(outp_smry_ls$predictors_lup, 
                target_var_nm_1L_chr = "mdl_scaling_dbl", match_value_xx = .x, 
                match_var_nm_1L_chr = "short_name_chr", evaluate_1L_lgl = F), 
                1))
            write_ts_mdl_plts(brms_mdl = model_mdl, table_predn_mdl = table_predn_mdl, 
                tfd_data_tb = outp_smry_ls$scored_data_tb %>% 
                  transform_tb_to_mdl_inp(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                    predr_vars_nms_chr = outp_smry_ls$predr_vars_nms_ls %>% 
                      purrr::flatten_chr() %>% unique(), id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
                    round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
                    round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr, 
                    scaling_fctr_dbl = scaling_fctr_dbl), depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                mdl_nm_1L_chr = mdl_nm_1L_chr, path_to_write_to_1L_chr = output_dir_chr[3], 
                predn_type_1L_chr = predn_type_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
                sd_dbl = sd_dbl, sfx_1L_chr = " from table", 
                tfmn_1L_chr = tfmn_1L_chr, utl_min_val_1L_dbl = ifelse(!is.null(outp_smry_ls$utl_min_val_1L_dbl), 
                  outp_smry_ls$utl_min_val_1L_dbl, -1))
            table_predn_mdl
        }) %>% stats::setNames(outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr())
    outp_smry_ls$shareable_mdls_ls <- shareable_mdls_ls
    outp_smry_ls$shareable_mdls_tb <- NULL
    ingredients_ls <- list(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
        dictionary_tb = outp_smry_ls$dictionary_tb %>% dplyr::filter(var_nm_chr %in% 
            names(fake_ds_tb)), id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
        fake_ds_tb = fake_ds_tb, mdls_lup = outp_smry_ls$shareable_mdls_ls %>% 
            purrr::map2_dfr(names(outp_smry_ls$shareable_mdls_ls), 
                ~{
                  if (inherits(.x, "betareg")) {
                    coeffs_dbl <- .x$coefficients$mean
                  } else {
                    coeffs_dbl <- .x$coefficients
                  }
                  mdl_type_1L_chr = get_mdl_type_from_nm(.y, 
                    mdl_types_lup = outp_smry_ls$mdl_types_lup)
                  tibble::tibble(mdl_nms_chr = .y) %>% dplyr::mutate(predrs_ls = list(coeffs_dbl %>% 
                    names() %>% stringr::str_remove_all("_change") %>% 
                    stringr::str_remove_all("_baseline") %>% 
                    unique() %>% purrr::discard(~.x == "(Intercept)")), 
                    mdl_type_chr = mdl_type_1L_chr, tfmn_chr = ready4::get_from_lup_obj(outp_smry_ls$mdl_types_lup, 
                      match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                      target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F))
                }), mdls_smry_tb = outp_smry_ls$mdls_smry_tb, 
        mdl_types_lup = mdl_types_lup, predictors_lup = outp_smry_ls$predictors_lup, 
        round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
        seed_1L_int = outp_smry_ls$seed_1L_int, utl_min_val_1L_dbl = ifelse(!is.null(outp_smry_ls$utl_min_val_1L_dbl), 
            outp_smry_ls$utl_min_val_1L_dbl, -1))
    saveRDS(ingredients_ls, paste0(output_dir_chr[2], "/mdl_ingredients", 
        ".RDS"))
    if (!is.null(outp_smry_ls$dv_ls)) {
        write_shareable_mdls_to_dv(outp_smry_ls, new_dir_nm_1L_chr = new_dir_nm_1L_chr, 
            share_ingredients_1L_lgl = T, output_dir_chr = output_dir_chr)
        if (write_mdls_to_dv_1L_lgl) {
            outp_smry_ls$shareable_mdls_tb <- write_shareable_mdls_to_dv(outp_smry_ls, 
                new_dir_nm_1L_chr = new_dir_nm_1L_chr, shareable_title_detail_1L_chr = shareable_title_detail_1L_chr, 
                share_ingredients_1L_lgl = F, output_dir_chr = output_dir_chr)
        }
    }
    return(outp_smry_ls)
}
#' Write shareable models to dataverse
#' @description write_shareable_mdls_to_dv() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write shareable models to dataverse. The function returns Shareable models (a tibble).
#' @param outp_smry_ls Output summary (a list)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'G_Shareable'
#' @param shareable_title_detail_1L_chr Shareable title detail (a character vector of length one), Default: ''
#' @param share_ingredients_1L_lgl Share ingredients (a logical vector of length one), Default: T
#' @param output_dir_chr Output directory (a character vector), Default: 'NA'
#' @return Shareable models (a tibble)
#' @rdname write_shareable_mdls_to_dv
#' @export 
#' @importFrom tibble tibble
#' @importFrom ready4 write_to_dv_with_wait get_fl_id_from_dv_ls
#' @importFrom dataverse get_dataset
#' @importFrom dplyr mutate
#' @importFrom purrr map_int
#' @keywords internal
write_shareable_mdls_to_dv <- function (outp_smry_ls, new_dir_nm_1L_chr = "G_Shareable", shareable_title_detail_1L_chr = "", 
    share_ingredients_1L_lgl = T, output_dir_chr = NA_character_) 
{
    if (is.na(output_dir_chr[1])) 
        output_dir_chr <- write_shareable_dir(outp_smry_ls = outp_smry_ls, 
            new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    if (share_ingredients_1L_lgl) {
        shareable_mdls_tb <- tibble::tibble(ds_obj_nm_chr = "mdl_ingredients", 
            title_chr = "An R object that can be used to construct model objects from tables of coefficients. Contains a synthetic dataset.")
    }
    else {
        shareable_mdls_tb <- tibble::tibble(ds_obj_nm_chr = names(outp_smry_ls$shareable_mdls_ls), 
            title_chr = paste0("A shareable (contains no confidential data) statistical model, ", 
                names(outp_smry_ls$shareable_mdls_ls), ".", shareable_title_detail_1L_chr))
    }
    ready4::write_to_dv_with_wait(shareable_mdls_tb, dv_nm_1L_chr = outp_smry_ls$dv_ls$dv_nm_1L_chr, 
        ds_url_1L_chr = outp_smry_ls$dv_ls$ds_url_1L_chr, parent_dv_dir_1L_chr = outp_smry_ls$dv_ls$parent_dv_dir_1L_chr, 
        paths_to_dirs_chr = output_dir_chr[ifelse(share_ingredients_1L_lgl, 
            2, 3)], paths_are_rltv_1L_lgl = F, inc_fl_types_chr = ".RDS")
    if (!share_ingredients_1L_lgl) {
        ds_ls <- dataverse::get_dataset(outp_smry_ls$dv_ls$ds_url_1L_chr)
        shareable_mdls_tb <- shareable_mdls_tb %>% dplyr::mutate(dv_nm_chr = outp_smry_ls$dv_ls$dv_nm_1L_chr, 
            fl_ids_int = ds_obj_nm_chr %>% purrr::map_int(~ready4::get_fl_id_from_dv_ls(ds_ls, 
                fl_nm_1L_chr = paste0(.x, ".RDS")) %>% as.integer()))
    }
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
#' @importFrom ready4 get_from_lup_obj
#' @importFrom dplyr arrange desc
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
            dictionary_tb = dictionary_tb), path_to_write_to_1L_chr = output_dir_1L_chr, 
        plt_nm_1L_chr = "A_TFMN_CMPRSN_DNSTY", height_1L_dbl = 6, 
        width_1L_dbl = 10)
    smry_of_sngl_predr_mdls_tb <- purrr::map_dfr(mdl_types_chr, 
        ~{
            tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .x, 
                target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F)
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
#' Write study output dataset
#' @description write_study_outp_ds() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write study output dataset. The function returns Dataverse dataset name and url (a character vector).
#' @param input_params_ls Input parameters (a list)
#' @param dv_ds_nm_and_url_chr Dataverse dataset name and url (a character vector), Default: NULL
#' @param rprt_lups_ls Report lookup tables (a list), Default: NULL
#' @param output_format_ls Output format (a list), Default: NULL
#' @param path_params_ls Path parameters (a list), Default: NULL
#' @param abstract_args_ls Abstract arguments (a list), Default: NULL
#' @param dv_mdl_desc_1L_chr Dataverse model description (a character vector of length one), Default: 'This is a longitudinal transfer to utility model designed for use with the youthu R package.'
#' @param header_yaml_args_ls Header yaml arguments (a list), Default: NULL
#' @param inc_fl_types_chr Include file types (a character vector), Default: '.pdf'
#' @param purge_data_1L_lgl Purge data (a logical vector of length one), Default: FALSE
#' @param start_at_int Start at (an integer vector), Default: c(2, 1)
#' @param use_fake_data_1L_lgl Use fake data (a logical vector of length one), Default: NULL
#' @return Dataverse dataset name and url (a character vector)
#' @rdname write_study_outp_ds
#' @export 
#' @importFrom purrr walk2
#' @importFrom stringr str_remove
#' @importFrom rlang exec
#' @importFrom dplyr filter
#' @importFrom ready4 get_from_lup_obj write_to_dv_with_wait
#' @importFrom tibble tibble
#' @keywords internal
write_study_outp_ds <- function (input_params_ls, dv_ds_nm_and_url_chr = NULL, rprt_lups_ls = NULL, 
    output_format_ls = NULL, path_params_ls = NULL, abstract_args_ls = NULL, 
    dv_mdl_desc_1L_chr = "This is a longitudinal transfer to utility model designed for use with the youthu R package.", 
    header_yaml_args_ls = NULL, inc_fl_types_chr = ".pdf", purge_data_1L_lgl = FALSE, 
    start_at_int = c(2, 1), use_fake_data_1L_lgl = NULL) 
{
    if (missing(header_yaml_args_ls)) {
        header_yaml_args_ls <- input_params_ls$header_yaml_args_ls
    }
    else {
        warning("The argument header_yaml_args_ls is soft deprecated. We recommend passing the header information as part of the list passed to the input_params_ls argument.")
    }
    if (missing(path_params_ls)) {
        path_params_ls <- input_params_ls$path_params_ls
    }
    else {
        warning("The argument path_params_ls is soft deprecated. We recommend passing the paths information as part of the list passed to the input_params_ls argument.")
    }
    if (missing(output_format_ls)) {
        output_format_ls <- input_params_ls$output_format_ls
    }
    else {
        warning("The argument output_format_ls is soft deprecated. We recommend passing the output format information as part of the list passed to the input_params_ls argument.")
    }
    if (missing(use_fake_data_1L_lgl)) {
        use_fake_data_1L_lgl <- input_params_ls$params_ls$use_fake_data_1L_lgl
    }
    else {
        warning("The argument use_fake_data_1L_lgl is soft deprecated. We recommend declaring whether dataset is fake as part of the list passed to the input_params_ls argument.")
    }
    if (missing(dv_ds_nm_and_url_chr)) {
        dv_ds_nm_and_url_chr <- input_params_ls$path_params_ls$dv_ds_nm_and_url_chr
    }
    else {
        warning("The argument dv_ds_nm_and_url_chr is soft deprecated. We recommend declaring dataverse details as part of the list passed to the input_params_ls argument.")
    }
    if (missing(rprt_lups_ls)) {
        rprt_lups_ls <- input_params_ls$rprt_lups_ls
    }
    else {
        warning("The argument rprt_lups_ls is soft deprecated. We recommend declaring report detail as part of the list passed to the input_params_ls argument.")
    }
    paths_ls <- path_params_ls$paths_ls
    rprt_lups_ls %>% purrr::walk2(names(rprt_lups_ls), ~{
        rprt_lup <- .x
        reference_1L_int <- ifelse(.y == "Primary", 0, as.numeric(stringr::str_remove(.y, 
            "secondary_")))
        if (is.null(rprt_lup)) {
            data("rprt_lup", package = "TTU", envir = environment())
            rprt_lup <- transform_rprt_lup(rprt_lup, add_suplry_rprt_1L_lgl = reference_1L_int > 
                0, add_sharing_rprt_1L_lgl = T, start_at_int = start_at_int, 
                reference_1L_int = reference_1L_int)
        }
        if (reference_1L_int == 0) {
            included_rprts_chr <- rprt_lup$rprt_nms_chr[rprt_lup$rprt_nms_chr != 
                "AAA_SHARING_MTH"]
            transform_paths_ls <- NULL
        }
        else {
            included_rprts_chr <- c("AAA_SUPLRY_ANLYS_MTH", paste0("AAA_TTU_MDL_CTG-", 
                reference_1L_int))[min(2, reference_1L_int):2]
            transform_paths_ls = list(fn = transform_paths_ls_for_scndry, 
                args_ls = list(reference_1L_int = reference_1L_int, 
                  remove_prmry_1L_lgl = T, remove_mkdn_1L_lgl = T))
            paths_ls <- rlang::exec(transform_paths_ls$fn, paths_ls, 
                !!!transform_paths_ls$args_ls)
        }
        params_ls <- list(dv_ds_nm_and_url_chr = dv_ds_nm_and_url_chr, 
            dv_mdl_desc_1L_chr = dv_mdl_desc_1L_chr, inc_fl_types_chr = inc_fl_types_chr, 
            nbr_of_digits_1L_int = output_format_ls$supplementary_digits_1L_int, 
            output_type_1L_chr = output_format_ls$supplementary_outp_1L_chr, 
            rprt_lup = rprt_lup %>% dplyr::filter(rprt_nms_chr %in% 
                included_rprts_chr), share_mdls_1L_lgl = (reference_1L_int == 
                0), subtitle_1L_chr = ready4::get_from_lup_obj(rprt_lup, 
                match_value_xx = "AAA_SHARING_MTH", match_var_nm_1L_chr = "rprt_nms_chr", 
                target_var_nm_1L_chr = "title_chr", evaluate_1L_lgl = F), 
            transform_paths_ls = transform_paths_ls, use_fake_data_1L_lgl = use_fake_data_1L_lgl) %>% 
            append(path_params_ls[1:2])
        params_ls %>% write_report(paths_ls = paths_ls, rprt_nm_1L_chr = "AAA_SHARING_MTH", 
            abstract_args_ls = abstract_args_ls, header_yaml_args_ls = header_yaml_args_ls, 
            rprt_lup = rprt_lup)
    })
    ready4::write_to_dv_with_wait(dss_tb = tibble::tibble(ds_obj_nm_chr = "AAA_SHARING_MTH", 
        title_chr = rprt_lups_ls[[1]] %>% ready4::get_from_lup_obj(match_value_xx = "AAA_SHARING_MTH", 
            match_var_nm_1L_chr = "rprt_nms_chr", target_var_nm_1L_chr = "title_chr", 
            evaluate_1L_lgl = F)), dv_nm_1L_chr = dv_ds_nm_and_url_chr[1], 
        ds_url_1L_chr = dv_ds_nm_and_url_chr[2], parent_dv_dir_1L_chr = paths_ls$dv_dir_1L_chr, 
        paths_to_dirs_chr = paths_ls$reports_dir_1L_chr, inc_fl_types_chr = inc_fl_types_chr, 
        paths_are_rltv_1L_lgl = F)
    return(dv_ds_nm_and_url_chr)
}
#' Write to delete dataset copies
#' @description write_to_delete_ds_copies() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write to delete dataset copies. The function is called for its side effects and does not return a value. WARNING: This function writes R scripts to your local environment. Make sure to only use if you want this behaviour
#' @param input_params_ls Input parameters (a list), Default: NULL
#' @param paths_ls Paths (a list), Default: NULL
#' @return NULL
#' @rdname write_to_delete_ds_copies
#' @export 
#' @importFrom purrr map_chr walk
#' @importFrom here here
write_to_delete_ds_copies <- function (input_params_ls = NULL, paths_ls = NULL) 
{
    if (is.null(paths_ls)) 
        paths_ls <- input_params_ls$path_params_ls$paths_ls
    paths_to_outp_chr <- c(paste0(paths_ls$output_data_dir_1L_chr, 
        "/I_ALL_OUTPUT_.RDS"))
    secondary_refs_int <- NULL
    if (!is.null(input_params_ls$scndry_anlys_params_ls)) {
        1:length(input_params_ls$scndry_anlys_params_ls)
        paths_to_outp_chr <- c(paths_to_outp_chr, secondary_refs_int %>% 
            purrr::map_chr(~here::here(paths_ls$path_from_top_level_1L_chr, 
                paths_ls$write_to_dir_nm_1L_chr, paste0("secondary_", 
                  .x), "Output", "I_ALL_OUTPUT_.RDS")))
    }
    paths_to_outp_chr %>% purrr::walk(~{
        outp_smry_ls <- readRDS(.x)
        write_to_delete_mdl_fls(outp_smry_ls)
        outp_smry_ls$scored_data_tb <- NULL
        saveRDS(outp_smry_ls, .x)
    })
}
#' Write to delete model files
#' @description write_to_delete_mdl_fls() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write to delete model files. The function is called for its side effects and does not return a value. WARNING: This function writes R scripts to your local environment. Make sure to only use if you want this behaviour
#' @param outp_smry_ls Output summary (a list)
#' @return NULL
#' @rdname write_to_delete_mdl_fls
#' @export 
#' @importFrom purrr map_lgl walk
write_to_delete_mdl_fls <- function (outp_smry_ls) 
{
    paths_to_mdls_chr <- outp_smry_ls$file_paths_chr[outp_smry_ls$file_paths_chr %>% 
        purrr::map_lgl(~endsWith(.x, ".RDS") & (startsWith(.x, 
            "A_Candidate_Mdls_Cmprsn") | startsWith(.x, "C_Predrs_Sngl_Mdl_Cmprsn") | 
            startsWith(.x, "D_Predr_Covars_Cmprsn") | startsWith(.x, 
            "E_Predrs_W_Covars_Sngl_Mdl_Cmprsn") | startsWith(.x, 
            "F_TS_Mdls")) & !endsWith(.x, "mdls_smry_tb.RDS"))]
    paths_to_mdls_chr %>% purrr::walk(~unlink(paste0(outp_smry_ls$path_to_write_to_1L_chr, 
        "/", .x)))
}
#' Write time series model plots
#' @description write_ts_mdl_plts() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write time series model plots. The function returns Model plots paths (a list).
#' @param brms_mdl Bayesian regression models (a model)
#' @param table_predn_mdl Table prediction (a model), Default: NULL
#' @param tfd_data_tb Transformed data (a tibble)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param depnt_var_desc_1L_chr Dependent variable description (a character vector of length one), Default: 'Utility score'
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param sd_dbl Standard deviation (a double vector), Default: NA
#' @param sfx_1L_chr Suffix (a character vector of length one), Default: ' from table'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param units_1L_chr Units (a character vector of length one), Default: 'in'
#' @param height_dbl Height (a double vector), Default: c(rep(6, 2), rep(5, 8))
#' @param width_dbl Width (a double vector), Default: c(rep(6, 2), rep(6, 8))
#' @param rsl_dbl Resolution (a double vector), Default: rep(300, 10)
#' @param args_ls Arguments (a list), Default: NULL
#' @param seed_1L_dbl Seed (a double vector of length one), Default: 23456
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: -1
#' @return Model plots paths (a list)
#' @rdname write_ts_mdl_plts
#' @export 
#' @importFrom purrr map discard
#' @importFrom ready4show write_mdl_plt_fl
#' @importFrom stats setNames
#' @keywords internal
write_ts_mdl_plts <- function (brms_mdl, table_predn_mdl = NULL, tfd_data_tb, mdl_nm_1L_chr, 
    path_to_write_to_1L_chr, depnt_var_nm_1L_chr = "utl_total_w", 
    depnt_var_desc_1L_chr = "Utility score", predn_type_1L_chr = NULL, 
    round_var_nm_1L_chr = "round", sd_dbl = NA_real_, sfx_1L_chr = " from table", 
    tfmn_1L_chr = "NTF", units_1L_chr = "in", height_dbl = c(rep(6, 
        2), rep(5, 8)), width_dbl = c(rep(6, 2), rep(6, 8)), 
    rsl_dbl = rep(300, 10), args_ls = NULL, seed_1L_dbl = 23456, 
    utl_min_val_1L_dbl = -1) 
{
    set.seed(seed_1L_dbl)
    tfd_data_tb <- transform_ds_for_all_cmprsn_plts(tfd_data_tb = tfd_data_tb, 
        model_mdl = brms_mdl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        is_brms_mdl_1L_lgl = inherits(brms_mdl, "brmsfit"), predn_type_1L_chr = predn_type_1L_chr, 
        sd_dbl = NA_real_, sfx_1L_chr = ifelse(!is.null(table_predn_mdl), 
            " from brmsfit", sfx_1L_chr), tfmn_1L_chr = tfmn_1L_chr, 
        utl_min_val_1L_dbl = utl_min_val_1L_dbl)
    if (!is.null(table_predn_mdl)) {
        tfd_data_tb <- transform_ds_for_all_cmprsn_plts(tfd_data_tb = tfd_data_tb, 
            model_mdl = table_predn_mdl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
            is_brms_mdl_1L_lgl = F, predn_type_1L_chr = predn_type_1L_chr, 
            sd_dbl = sd_dbl, sfx_1L_chr = ifelse(!is.null(brms_mdl), 
                " from table", sfx_1L_chr), tfmn_1L_chr = tfmn_1L_chr, 
            utl_min_val_1L_dbl = utl_min_val_1L_dbl)
    }
    plt_nms_chr <- paste0(mdl_nm_1L_chr, "_", c("coefs", "hetg", 
        "dnst", "sctr_plt", "sim_dnst", "sim_sctr", "cnstrd_dnst", 
        "cnstrd_sctr_plt", "cnstrd_sim_dnst", "cnstrd_sim_sctr"))
    mdl_plts_paths_ls <- purrr::map(ifelse(inherits(brms_mdl, 
        "brmsfit"), 1, 3):10, ~{
        plt_fn <- fn_args_ls <- NULL
        if (.x %in% c(1, 2)) {
            plt <- plot(brms_mdl, ask = F, plot = F)
            if (length(plt) >= .x) {
                fn_args_ls <- list(brms_mdl = brms_mdl, idx_1L_int = as.integer(.x))
                plt_fn <- function(brms_mdl, idx_1L_int) {
                  plot(brms_mdl, ask = F, plot = F)[idx_1L_int]
                }
            }
        }
        else {
            if (.x %in% c(3, 5, 7, 9)) {
                plt_fn <- plot_obsd_predd_dnst
                fn_args_ls <- list(tfd_data_tb = tfd_data_tb, 
                  depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                  depnt_var_desc_1L_chr = depnt_var_desc_1L_chr, 
                  predd_val_var_nm_1L_chr = ifelse(.x %in% c(3, 
                    7), transform_predd_var_nm("Predicted", sfx_1L_chr = ifelse(!is.null(table_predn_mdl), 
                    " from brmsfit", sfx_1L_chr), utl_min_val_1L_dbl = ifelse(.x == 
                    3, NA_real_, utl_min_val_1L_dbl)), transform_predd_var_nm("Simulated", 
                    sfx_1L_chr = ifelse(!is.null(table_predn_mdl), 
                      " from brmsfit", sfx_1L_chr), utl_min_val_1L_dbl = ifelse(.x == 
                      5, NA_real_, utl_min_val_1L_dbl))), cmprsn_predd_var_nm_1L_chr = ifelse(is.null(table_predn_mdl), 
                    NA_character_, ifelse(.x %in% c(3, 7), transform_predd_var_nm("Predicted", 
                      sfx_1L_chr = " from table", utl_min_val_1L_dbl = ifelse(.x == 
                        3, NA_real_, utl_min_val_1L_dbl)), transform_predd_var_nm("Simulated", 
                      sfx_1L_chr = " from table", utl_min_val_1L_dbl = ifelse(.x == 
                        5, NA_real_, utl_min_val_1L_dbl)))))
            }
            else {
                plt_fn <- plot_obsd_predd_sctr_cmprsn
                fn_args_ls <- list(tfd_data_tb = tfd_data_tb, 
                  depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                  depnt_var_desc_1L_chr = depnt_var_desc_1L_chr, 
                  round_var_nm_1L_chr = round_var_nm_1L_chr, 
                  predd_val_var_nm_1L_chr = ifelse(.x %in% c(4, 
                    8), transform_predd_var_nm("Predicted", sfx_1L_chr = ifelse(!is.null(table_predn_mdl), 
                    " from brmsfit", sfx_1L_chr), utl_min_val_1L_dbl = ifelse(.x == 
                    4, NA_real_, utl_min_val_1L_dbl)), transform_predd_var_nm("Simulated", 
                    sfx_1L_chr = ifelse(!is.null(table_predn_mdl), 
                      " from brmsfit", sfx_1L_chr), utl_min_val_1L_dbl = ifelse(.x == 
                      6, NA_real_, utl_min_val_1L_dbl))), args_ls = args_ls)
            }
        }
        ready4show::write_mdl_plt_fl(plt_fn, fn_args_ls = fn_args_ls, 
            path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
            plt_nm_1L_chr = plt_nms_chr[.x], units_1L_chr = units_1L_chr, 
            width_1L_dbl = width_dbl[.x], height_1L_dbl = height_dbl[.x], 
            rsl_1L_dbl = rsl_dbl[.x])
    }) %>% stats::setNames(plt_nms_chr[ifelse(inherits(brms_mdl, 
        "brmsfit"), 1, 3):10]) %>% purrr::discard(is.na)
    return(mdl_plts_paths_ls)
}
#' Write time series models
#' @description write_ts_mdls() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write time series models. The function returns Models summary (a tibble).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_vars_nms_ls Predictor variables names (a list)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one), Default: 'Baseline'
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: -1
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param mdl_nms_ls Model names (a list)
#' @param mdl_smry_dir_1L_chr Model summary directory (a character vector of length one)
#' @param predictors_lup Predictors (a lookup table)
#' @param iters_1L_int Iterations (an integer vector of length one), Default: 4000
#' @param mdl_types_lup Model types (a lookup table)
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @param prior_ls Prior (a list), Default: NULL
#' @param control_ls Control (a list), Default: NULL
#' @return Models summary (a tibble)
#' @rdname write_ts_mdls
#' @export 
#' @importFrom purrr map_dfr
#' @keywords internal
write_ts_mdls <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_ls, 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline", utl_min_val_1L_dbl = -1, 
    backend_1L_chr = getOption("brms.backend", "rstan"), mdl_nms_ls, 
    mdl_smry_dir_1L_chr, predictors_lup, iters_1L_int = 4000L, 
    mdl_types_lup, seed_1L_int = 1000L, prior_ls = NULL, control_ls = NULL) 
{
    if (!dir.exists(mdl_smry_dir_1L_chr)) 
        dir.create(mdl_smry_dir_1L_chr)
    mdls_smry_tb <- purrr::map_dfr(1:length(mdl_nms_ls), ~{
        idx_1L_int <- .x
        purrr::map_dfr(mdl_nms_ls[[idx_1L_int]], ~{
            smry_ls <- make_smry_of_ts_mdl_outp(data_tb = data_tb, 
                predr_vars_nms_chr = predr_vars_nms_ls[[idx_1L_int]], 
                mdl_nm_1L_chr = .x, path_to_write_to_1L_chr = mdl_smry_dir_1L_chr, 
                depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
                round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr, 
                predictors_lup = predictors_lup, utl_min_val_1L_dbl = utl_min_val_1L_dbl, 
                backend_1L_chr = backend_1L_chr, iters_1L_int = iters_1L_int, 
                mdl_types_lup = mdl_types_lup, seed_1L_int = seed_1L_int, 
                prior_ls = prior_ls, control_ls = control_ls)
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
#' @param predictors_lup Predictors (a lookup table)
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: -1
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iterations (an integer vector of length one), Default: 4000
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'F_TS_Mdls'
#' @param prior_ls Prior (a list), Default: NULL
#' @param control_ls Control (a list), Default: NULL
#' @return Output summary (a list)
#' @rdname write_ts_mdls_from_alg_outp
#' @export 
write_ts_mdls_from_alg_outp <- function (outp_smry_ls, predictors_lup, utl_min_val_1L_dbl = -1, 
    backend_1L_chr = getOption("brms.backend", "rstan"), iters_1L_int = 4000L, 
    new_dir_nm_1L_chr = "F_TS_Mdls", prior_ls = NULL, control_ls = NULL) 
{
    output_dir_1L_chr <- write_new_outp_dir(outp_smry_ls$path_to_write_to_1L_chr, 
        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
    mdls_smry_tb <- write_ts_mdls(data_tb = outp_smry_ls$scored_data_tb, 
        depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
        predr_vars_nms_ls = outp_smry_ls$predr_vars_nms_ls, id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
        round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
        round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr, 
        mdl_nms_ls = outp_smry_ls$mdl_nms_ls, mdl_smry_dir_1L_chr = output_dir_1L_chr, 
        predictors_lup = predictors_lup, utl_min_val_1L_dbl = utl_min_val_1L_dbl, 
        backend_1L_chr = backend_1L_chr, iters_1L_int = iters_1L_int, 
        mdl_types_lup = outp_smry_ls$mdl_types_lup, seed_1L_int = outp_smry_ls$seed_1L_int, 
        prior_ls = prior_ls, control_ls = control_ls)
    outp_smry_ls$mdls_smry_tb <- mdls_smry_tb
    outp_smry_ls$utl_min_val_1L_dbl <- utl_min_val_1L_dbl
    outp_smry_ls$file_paths_chr <- list.files(outp_smry_ls$path_to_write_to_1L_chr, 
        recursive = T)
    return(outp_smry_ls)
}
