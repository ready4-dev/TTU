#' Make analysis dataset summary
#' @description make_analysis_ds_smry_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make analysis dataset summary list. The function returns Analysis dataset summary (a list).
#' @param ds_descvs_ls Dataset descriptives (a list)
#' @param candidate_covar_nms_chr Candidate covariate names (a character vector)
#' @param predictors_lup Predictors (a lookup table)
#' @return Analysis dataset summary (a list)
#' @rdname make_analysis_ds_smry_ls
#' @export 

#' @keywords internal
make_analysis_ds_smry_ls <- function (ds_descvs_ls, candidate_covar_nms_chr, predictors_lup) 
{
    analysis_ds_smry_ls <- list(candidate_predrs_chr = ds_descvs_ls$candidate_predrs_chr, 
        candidate_covar_nms_chr = candidate_covar_nms_chr, depnt_var_nm_1L_chr = ds_descvs_ls$utl_wtd_var_nm_1L_chr, 
        id_var_nm_1L_chr = ds_descvs_ls$id_var_nm_1L_chr, predictors_lup = predictors_lup, 
        round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr, 
        round_bl_val_1L_chr = ds_descvs_ls$round_vals_chr[1], 
        dictionary_tb = ds_descvs_ls$dictionary_tb)
    return(analysis_ds_smry_ls)
}
#' Make bayesian regression models model print list
#' @description make_brms_mdl_print_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make bayesian regression models model print list. The function returns Bayesian regression models model print (a list).
#' @param mdl_ls Model list (a list of models)
#' @param label_stub_1L_chr Label stub (a character vector of length one)
#' @param caption_1L_chr Caption (a character vector of length one)
#' @param output_type_1L_chr Output type (a character vector of length one), Default: 'PDF'
#' @param digits_1L_dbl Digits (a double vector of length one), Default: 2
#' @param big_mark_1L_chr Big mark (a character vector of length one), Default: ' '
#' @return Bayesian regression models model print (a list)
#' @rdname make_brms_mdl_print_ls
#' @export 
#' @importFrom utils capture.output
#' @importFrom purrr map_dbl map_chr map2_chr discard map
#' @importFrom dplyr mutate all_of across case_when
#' @importFrom Hmisc latexTranslate
#' @importFrom stringr str_replace
#' @keywords internal
make_brms_mdl_print_ls <- function (mdl_ls, label_stub_1L_chr, caption_1L_chr, output_type_1L_chr = "PDF", 
    digits_1L_dbl = 2, big_mark_1L_chr = " ") 
{
    smry_mdl_ls <- summary(mdl_ls, digits = 4)
    mdl_smry_chr <- smry_mdl_ls %>% utils::capture.output()
    idx_dbl <- c("Formula: ", "Samples: ", "Group-Level Effects: ", 
        "Population-Level Effects: ", "Family Specific Parameters: ", 
        "Samples were drawn using ") %>% purrr::map_dbl(~mdl_smry_chr %>% 
        startsWith(.x) %>% which())
    data_tb <- make_brms_mdl_smry_tbl(smry_mdl_ls, grp_1L_chr = mdl_smry_chr[idx_dbl[3]], 
        popl_1L_chr = mdl_smry_chr[idx_dbl[4]], fam_1L_chr = mdl_smry_chr[idx_dbl[5]])
    bold_lgl <- data_tb$Parameter %in% c(mdl_smry_chr[idx_dbl[3]], 
        mdl_smry_chr[idx_dbl[4]], mdl_smry_chr[idx_dbl[5]])
    if (output_type_1L_chr == "PDF") {
        data_tb <- data_tb %>% dplyr::mutate(Parameter = purrr::map_chr(Parameter, 
            ~.x %>% Hmisc::latexTranslate()))
    }
    data_tb <- data_tb %>% dplyr::mutate(Parameter = Parameter %>% 
        purrr::map2_chr(dplyr::all_of(bold_lgl), ~ifelse(.y & 
            output_type_1L_chr == "PDF", paste0("\\textbf{", 
            .x, "}"), .x)))
    if (output_type_1L_chr != "PDF") {
        data_tb <- data_tb %>% dplyr::mutate(dplyr::across(c(Bulk_ESS, 
            Tail_ESS), ~format(., big.mark = big_mark_1L_chr)))
    }
    if (output_type_1L_chr == "HTML") {
        data_tb <- data_tb %>% dplyr::mutate(dplyr::across(where(is.numeric), 
            ~format(round(., digits = digits_1L_dbl), digits = digits_1L_dbl, 
                nsmall = digits_1L_dbl)))
    }
    data_tb <- data_tb %>% dplyr::mutate(dplyr::across(where(is.character), 
        ~dplyr::case_when(is.na(.) ~ "", . == "NA" ~ "", endsWith(., 
            " NA") ~ "", TRUE ~ .)))
    end_matter_1L_chr <- trimws(mdl_smry_chr[idx_dbl[6]:length(mdl_smry_chr)]) %>% 
        paste0(collapse = " ")
    brms_mdl_print_ls <- list(part_1 = mdl_smry_chr[idx_dbl[1]], 
        part_2 = "\n\n", part_3 = c(trimws(mdl_smry_chr[1:(idx_dbl[2] - 
            1)][-idx_dbl[1]]), paste0(trimws(mdl_smry_chr[idx_dbl[2]]), 
            " ", trimws(mdl_smry_chr[idx_dbl[2] + 1]), collapse = " ")) %>% 
            paste0(collapse = ifelse(output_type_1L_chr == "PDF", 
                "\n\n", "\n")), part_4 = "\n\n", part_5 = list(data_tb = data_tb, 
            output_type_1L_chr = output_type_1L_chr, caption_1L_chr = caption_1L_chr, 
            mkdn_tbl_ref_1L_chr = paste0("tab:", label_stub_1L_chr), 
            merge_row_idx_int = as.integer(which(bold_lgl)), 
            digits_dbl = c(ifelse(output_type_1L_chr == "PDF", 
                0, NA_real_) %>% purrr::discard(is.na), names(data_tb) %>% 
                purrr::map_dbl(~ifelse(.x %in% c("Bulk_ESS", 
                  "Tail_ESS"), 0, digits_1L_dbl))), big_mark_1L_chr = big_mark_1L_chr, 
            hline_after_ls = c(-1, 0), sanitize_fn = force, footnotes_chr = NA_character_), 
        part_6 = end_matter_1L_chr)
    if (output_type_1L_chr != "PDF") {
        brms_mdl_print_ls$part_5$footnotes_chr <- c(paste0(brms_mdl_print_ls$part_1, 
            ifelse(output_type_1L_chr == "Word", "", "\n")), 
            brms_mdl_print_ls$part_3, brms_mdl_print_ls$part_6)
        brms_mdl_print_ls$part_6 <- NULL
    }
    else {
        footnotes_chr <- c(mdl_smry_chr[idx_dbl[1]], trimws(mdl_smry_chr[1:(idx_dbl[2] - 
            1)][-idx_dbl[1]]), trimws(mdl_smry_chr[idx_dbl[2]]), 
            trimws(mdl_smry_chr[idx_dbl[2] + 1]), trimws(mdl_smry_chr[idx_dbl[6]:length(mdl_smry_chr)])) %>% 
            Hmisc::latexTranslate()
        footnotes_chr[1] <- footnotes_chr[1] %>% stringr::str_replace("~", 
            "\\\\textasciitilde")
        brms_mdl_print_ls$part_5$add_to_row_ls <- list(pos = purrr::map(c(0, 
            rep(nrow(data_tb), length(footnotes_chr) + 1)), ~.x), 
            command = c(names(data_tb) %>% Hmisc::latexTranslate() %>% 
                paste0(collapse = " & ") %>% paste0("\\\\\n"), 
                c("\\toprule\n"), footnotes_chr %>% purrr::map_chr(~paste0("\\multicolumn{", 
                  ncol(data_tb), "}{l}{", paste0("{\\footnotesize ", 
                    .x, "}\n", collapse = ","), "}\\\\\n"))))
    }
    return(brms_mdl_print_ls)
}
#' Make bayesian regression models model summary table
#' @description make_brms_mdl_smry_tbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make bayesian regression models model summary table. The function returns Bayesian regression models model summary (a tibble).
#' @param smry_mdl_ls Summary (a list of models)
#' @param grp_1L_chr Group (a character vector of length one)
#' @param popl_1L_chr Population (a character vector of length one)
#' @param fam_1L_chr Fam (a character vector of length one)
#' @return Bayesian regression models model summary (a tibble)
#' @rdname make_brms_mdl_smry_tbl
#' @export 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @keywords internal
make_brms_mdl_smry_tbl <- function (smry_mdl_ls, grp_1L_chr, popl_1L_chr, fam_1L_chr) 
{
    brms_mdl_smry_tb <- purrr::map(1:length(smry_mdl_ls$random), 
        ~make_mdl_smry_elmt_tbl(ctg_chr = c(ifelse(.x == 1, grp_1L_chr, 
            character(0)), paste0(names(smry_mdl_ls$ngrps)[.x], 
            " (Number of levels: ", smry_mdl_ls$ngrps[.x][[1]], 
            ")")), mat = smry_mdl_ls$random[.x][[1]])) %>% dplyr::bind_rows(make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$fixed, 
        ctg_chr = popl_1L_chr), make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$spec_pars, 
        ctg_chr = fam_1L_chr))
    return(brms_mdl_smry_tb)
}
#' Make fake time series data
#' @description make_fake_ts_data() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make fake time series data. The function returns Fk data (a tibble).
#' @param outp_smry_ls Output summary (a list)
#' @return Fk data (a tibble)
#' @rdname make_fake_ts_data
#' @export 
#' @importFrom synthpop syn
#' @importFrom purrr map_lgl
#' @importFrom dplyr mutate across all_of
make_fake_ts_data <- function (outp_smry_ls) 
{
    data_tb <- outp_smry_ls$scored_data_tb %>% transform_tb_to_mdl_inp(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
        predr_vars_nms_chr = outp_smry_ls$predr_cmprsn_tb$predr_chr, 
        id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
        round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr)
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = names(data_tb)[names(data_tb) != 
        outp_smry_ls$id_var_nm_1L_chr], seed = outp_smry_ls$seed_1L_int)
    dep_vars_chr <- names(fk_data_ls$syn)[names(fk_data_ls$syn) %>% 
        purrr::map_lgl(~startsWith(.x, outp_smry_ls$depnt_var_nm_1L_chr))]
    fk_data_tb <- fk_data_ls$syn %>% dplyr::mutate(dplyr::across(dplyr::all_of(dep_vars_chr), 
        ~NA_real_))
    return(fk_data_tb)
}
#' Make folds
#' @description make_folds_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make folds list. The function returns Folds (a list).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param folds_1L_int Folds (an integer vector of length one), Default: 10
#' @return Folds (a list)
#' @rdname make_folds_ls
#' @export 
#' @importFrom caret createFolds
#' @importFrom dplyr pull
#' @importFrom rlang sym
#' @keywords internal
make_folds_ls <- function (data_tb, depnt_var_nm_1L_chr = "aqol6d_total_w", folds_1L_int = 10L) 
{
    folds_ls <- caret::createFolds(data_tb %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)), 
        k = folds_1L_int, list = TRUE, returnTrain = FALSE)
    return(folds_ls)
}
#' Make knit parameters
#' @description make_knit_pars_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make knit parameters list. The function returns Knit parameters (a list).
#' @param rltv_path_to_data_dir_1L_chr Relative path to data directory (a character vector of length one)
#' @param mdl_types_chr Model types (a character vector)
#' @param predr_vars_nms_ls Predictor variables names (a list)
#' @param output_type_1L_chr Output type (a character vector of length one), Default: 'HTML'
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param plt_types_lup Plot types (a lookup table), Default: NULL
#' @param plt_types_chr Plot types (a character vector), Default: 'NA'
#' @param section_type_1L_chr Section type (a character vector of length one), Default: '#'
#' @return Knit parameters (a list)
#' @rdname make_knit_pars_ls
#' @export 
#' @importFrom utils data
#' @importFrom purrr pmap map map_chr flatten_chr map2
#' @importFrom stringr str_detect
#' @importFrom stats setNames
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
make_knit_pars_ls <- function (rltv_path_to_data_dir_1L_chr, mdl_types_chr, predr_vars_nms_ls, 
    output_type_1L_chr = "HTML", mdl_types_lup = NULL, plt_types_lup = NULL, 
    plt_types_chr = NA_character_, section_type_1L_chr = "#") 
{
    if (is.null(mdl_types_lup)) 
        utils::data(mdl_types_lup, envir = environment())
    if (is.null(plt_types_lup)) 
        utils::data(plt_types_lup, envir = environment())
    if (is.na(plt_types_chr)) {
        plt_types_chr <- plt_types_lup$short_name_chr
    }
    paths_to_all_data_fls_chr <- list.files(rltv_path_to_data_dir_1L_chr, 
        full.names = T)
    lab_idx_dbl <- 1:(length(mdl_types_chr) * length(predr_vars_nms_ls))
    knit_pars_ls <- purrr::pmap(list(predr_vars_nms_ls, split(lab_idx_dbl, 
        ceiling(seq_along(lab_idx_dbl)/length(mdl_types_chr))), 
        make_mdl_nms_ls(predr_vars_nms_ls, mdl_types_chr = mdl_types_chr)), 
        ~{
            mdl_nms_chr <- ..3
            mdl_data_paths_ls <- mdl_nms_chr %>% purrr::map(~paths_to_all_data_fls_chr[stringr::str_detect(paths_to_all_data_fls_chr, 
                .x)]) %>% stats::setNames(mdl_nms_chr)
            paths_to_mdls_chr <- mdl_data_paths_ls %>% purrr::map_chr(~.x[endsWith(.x, 
                ".RDS")]) %>% unname()
            paths_to_mdl_plts_ls <- mdl_data_paths_ls %>% purrr::map(~{
                paths_to_all_plots_chr <- .x[endsWith(.x, ".png")]
                plt_types_chr %>% purrr::map(~paths_to_all_plots_chr[paths_to_all_plots_chr %>% 
                  stringr::str_detect(.x)]) %>% purrr::flatten_chr()
            })
            mdl_ttls_chr <- paste0(..1[1], ifelse(is.na(..1[2]), 
                "", paste(" with ", ..1[2])), " ", mdl_types_chr %>% 
                purrr::map_chr(~ready4fun::get_from_lup_obj(mdl_types_lup, 
                  match_var_nm_1L_chr = "short_name_chr", match_value_xx = .x, 
                  target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F)))
            section_ttls_chr <- paste0(section_type_1L_chr, " ", 
                mdl_ttls_chr)
            plt_nms_ls <- paths_to_mdl_plts_ls %>% purrr::map2(mdl_ttls_chr, 
                ~{
                  paths_to_mdl_plts_chr <- .x
                  mdl_ttl_1L_chr <- .y
                  transform_1L_lgl <- paths_to_mdl_plts_chr %>% 
                    endsWith("_coefs.png") %>% any()
                  if (paths_to_mdl_plts_chr %>% endsWith("_hetg.png") %>% 
                    any()) 
                    transform_1L_lgl <- F
                  plt_types_chr %>% purrr::map(~{
                    if (endsWith(paths_to_mdl_plts_chr, paste0("_", 
                      .x, ".png")) %>% any()) {
                      paste0(mdl_ttl_1L_chr, " ", ifelse(transform_1L_lgl & 
                        .x == "coefs", "population and group level effects", 
                        ready4fun::get_from_lup_obj(plt_types_lup, 
                          match_var_nm_1L_chr = "short_name_chr", 
                          match_value_xx = .x, target_var_nm_1L_chr = "long_name_chr", 
                          evaluate_lgl = F)))
                    }
                    else {
                      character(0)
                    }
                  }) %>% purrr::flatten_chr()
                })
            list(plt_nms_ls = plt_nms_ls, paths_to_mdls_chr = paths_to_mdls_chr, 
                tbl_captions_chr = mdl_ttls_chr, label_stubs_chr = paste0("lab", 
                  ..2), output_type_1L_chr = output_type_1L_chr, 
                section_ttls_chr = section_ttls_chr, paths_to_mdl_plts_ls = paths_to_mdl_plts_ls)
        }) %>% stats::setNames(predr_vars_nms_ls %>% purrr::map_chr(~paste(.x, 
        collapse = "_")))
    return(knit_pars_ls)
}
#' Make
#' @description make_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make model. The function returns Model (a model).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param covar_var_nms_chr Covariate variable names (a character vector), Default: 'NA'
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
#' @keywords internal
make_mdl <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF", 
    predr_var_nm_1L_chr, covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF", 
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
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
        transform_depnt_var_nm(depnt_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr), 
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
#' @param predr_vars_nms_ls Predictor variables names (a list)
#' @param mdl_types_chr Model types (a character vector)
#' @return Model names (a list)
#' @rdname make_mdl_nms_ls
#' @export 
#' @importFrom purrr map2
#' @keywords internal
make_mdl_nms_ls <- function (predr_vars_nms_ls, mdl_types_chr) 
{
    mdl_nms_ls <- purrr::map2(predr_vars_nms_ls, make_unique_ls_elmt_idx_int(predr_vars_nms_ls), 
        ~paste0(.x[1], "_", ifelse(is.na(.x[2]), "", paste0(.x[2], 
            "_")), .y, "_", mdl_types_chr))
    return(mdl_nms_ls)
}
#' Make model summary element table
#' @description make_mdl_smry_elmt_tbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make model summary element table. The function returns Model element sum (a tibble).
#' @param mat Matrix (a matrix)
#' @param ctg_chr Category categories (a character vector)
#' @return Model element sum (a tibble)
#' @rdname make_mdl_smry_elmt_tbl
#' @export 
#' @importFrom tibble as_tibble add_case
#' @importFrom dplyr mutate select everything filter bind_rows
#' @keywords internal
make_mdl_smry_elmt_tbl <- function (mat, ctg_chr) 
{
    tb <- mat %>% tibble::as_tibble() %>% dplyr::mutate(Parameter = rownames(mat)) %>% 
        dplyr::select(Parameter, dplyr::everything())
    mdl_elmt_sum_tb <- tb %>% dplyr::filter(F) %>% tibble::add_case(Parameter = ctg_chr) %>% 
        dplyr::bind_rows(tb)
    return(mdl_elmt_sum_tb)
}
#' Make prediction dataset with one predictor
#' @description make_predn_ds_with_one_predr() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make prediction dataset with one predictor. The function returns Prediction dataset (a tibble).
#' @param model_mdl Model (a model)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param predr_vals_dbl Predictor values (a double vector)
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @return Prediction dataset (a tibble)
#' @rdname make_predn_ds_with_one_predr
#' @export 
#' @importFrom tibble tibble
#' @importFrom rlang sym
#' @importFrom dplyr mutate
#' @importFrom stats predict
#' @keywords internal
make_predn_ds_with_one_predr <- function (model_mdl, depnt_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF", 
    predr_var_nm_1L_chr, predr_vals_dbl, predn_type_1L_chr = NULL) 
{
    predn_ds_tb <- tibble::tibble(`:=`(!!rlang::sym(predr_var_nm_1L_chr), 
        predr_vals_dbl))
    predn_ds_tb <- predn_ds_tb %>% dplyr::mutate(`:=`(!!rlang::sym(depnt_var_nm_1L_chr), 
        stats::predict(model_mdl, newdata = predn_ds_tb, type = predn_type_1L_chr) %>% 
            calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, 
                tfmn_is_outp_1L_lgl = T)))
    return(predn_ds_tb)
}
#' Make predictor values
#' @description make_predr_vals() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make predictor values. The function returns Predictor values (a double vector).
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param candidate_predrs_lup Candidate predictors (a lookup table), Default: NULL
#' @return Predictor values (a double vector)
#' @rdname make_predr_vals
#' @export 
#' @importFrom utils data
#' @importFrom purrr map_dbl
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom rlang exec
#' @keywords internal
make_predr_vals <- function (predr_var_nm_1L_chr, candidate_predrs_lup = NULL) 
{
    if (is.null(candidate_predrs_lup)) {
        utils::data("candidate_predrs_lup", envir = environment())
    }
    args_ls <- purrr::map_dbl(names(candidate_predrs_lup)[3:4], 
        ~candidate_predrs_lup %>% ready4fun::get_from_lup_obj(match_value_xx = predr_var_nm_1L_chr, 
            match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = .x, 
            evaluate_lgl = F)) %>% as.list()
    predr_vals_dbl <- rlang::exec(.fn = seq, !!!args_ls)
    return(predr_vals_dbl)
}
#' Make predictor variables names
#' @description make_predr_vars_nms_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make predictor variables names list. The function returns Predictor variables names (a list).
#' @param main_predrs_chr Main predictors (a character vector)
#' @param covars_ls Covariates (a list)
#' @return Predictor variables names (a list)
#' @rdname make_predr_vars_nms_ls
#' @export 
#' @importFrom purrr map flatten
#' @keywords internal
make_predr_vars_nms_ls <- function (main_predrs_chr, covars_ls) 
{
    predr_vars_nms_ls <- covars_ls %>% purrr::map(~{
        covars_chr <- .x
        purrr::map(main_predrs_chr, ~list(c(.x), c(.x, covars_chr))) %>% 
            purrr::flatten()
    }) %>% purrr::flatten()
    predr_vars_nms_ls <- predr_vars_nms_ls[order(sapply(predr_vars_nms_ls, 
        length))]
    return(predr_vars_nms_ls)
}
#' Make preferred models vector
#' @description make_prefd_mdls_vec() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make preferred models vector. The function returns Preferred models (a character vector).
#' @param smry_of_sngl_predr_mdls_tb Summary of single predictor models (a tibble)
#' @param choose_from_pfx_chr Choose from prefix (a character vector), Default: c("GLM", "OLS")
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @return Preferred models (a character vector)
#' @rdname make_prefd_mdls_vec
#' @export 
#' @importFrom utils data
#' @importFrom dplyr inner_join select rename pull
#' @importFrom purrr map_chr
#' @keywords internal
make_prefd_mdls_vec <- function (smry_of_sngl_predr_mdls_tb, choose_from_pfx_chr = c("GLM", 
    "OLS"), mdl_types_lup = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    ordered_mdl_types_chr <- dplyr::inner_join(smry_of_sngl_predr_mdls_tb %>% 
        dplyr::select(Model) %>% dplyr::rename(long_name_chr = Model), 
        mdl_types_lup) %>% dplyr::pull(short_name_chr)
    prefd_mdls_chr <- purrr::map_chr(choose_from_pfx_chr, ~ordered_mdl_types_chr[startsWith(ordered_mdl_types_chr, 
        .x)][1])
    return(prefd_mdls_chr)
}
#' Make shareable
#' @description make_shareable_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make shareable model. The function returns Model (a model).
#' @param data_tb Data (a tibble)
#' @param mdl_smry_tb Model summary (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
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
#' @keywords internal
make_shareable_mdl <- function (data_tb, mdl_smry_tb, depnt_var_nm_1L_chr = "utl_total_w", 
    id_var_nm_1L_chr = "fkClientID", tfmn_1L_chr = "CLL", mdl_type_1L_chr = "OLS_CLL", 
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NA_character_, 
    seed_1L_int = 12345L) 
{
    if (is.null(mdl_types_lup)) 
        utils::data(mdl_types_lup, envir = environment())
    all_var_nms_chr <- names(data_tb)
    tfd_depnt_var_nm_1L_chr <- all_var_nms_chr[all_var_nms_chr %>% 
        startsWith(depnt_var_nm_1L_chr)]
    predr_var_nms_chr <- setdiff(all_var_nms_chr, c(tfd_depnt_var_nm_1L_chr, 
        id_var_nm_1L_chr))
    if (length(predr_var_nms_chr) > 1) {
        covar_var_nms_chr <- predr_var_nms_chr[2:length(predr_var_nms_chr)]
    }
    else {
        covar_var_nms_chr <- NA_character_
    }
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = all_var_nms_chr[all_var_nms_chr != 
        id_var_nm_1L_chr], seed = seed_1L_int)
    model_mdl <- make_mdl(fk_data_ls$syn, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
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
#' Make summary of bayesian regression model
#' @description make_smry_of_brm_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make summary of bayesian regression model model. The function returns Summary of bayesian regression model model (a tibble).
#' @param mdl_ls Model list (a list of models)
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param fn Function (a function), Default: calculate_rmse
#' @param mdl_nm_1L_chr Model name (a character vector of length one), Default: 'NA'
#' @param seed_1L_dbl Seed (a double vector of length one), Default: 23456
#' @return Summary of bayesian regression model model (a tibble)
#' @rdname make_smry_of_brm_mdl
#' @export 
#' @importFrom stats predict
#' @importFrom brms bayes_R2
#' @importFrom psych describe
#' @importFrom dplyr pull mutate rename select
#' @importFrom rlang sym
#' @importFrom purrr map flatten_chr
#' @keywords internal
make_smry_of_brm_mdl <- function (mdl_ls, data_tb, depnt_var_nm_1L_chr = "utl_total_w", 
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
        dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))), quant = c(0.25, 
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
#' Make summary of model output
#' @description make_smry_of_mdl_outp() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make summary of model output. The function returns Summary of one predictor model (a tibble).
#' @param data_tb Data (a tibble)
#' @param model_mdl Model (a model)
#' @param folds_1L_int Folds (an integer vector of length one), Default: 10
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param start_1L_chr Start (a character vector of length one), Default: NULL
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param covar_var_nms_chr Covariate variable names (a character vector), Default: 'NA'
#' @param mdl_type_1L_chr Model type (a character vector of length one), Default: 'OLS_NTF'
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @return Summary of one predictor model (a tibble)
#' @rdname make_smry_of_mdl_outp
#' @export 
#' @importFrom utils data
#' @importFrom dplyr filter pull summarise_all mutate select everything
#' @importFrom rlang sym
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom purrr map_dfr
#' @importFrom stats predict
#' @importFrom tibble tibble
#' @importFrom caret R2 RMSE MAE
#' @keywords internal
make_smry_of_mdl_outp <- function (data_tb, model_mdl, folds_1L_int = 10, depnt_var_nm_1L_chr = "utl_total_w", 
    start_1L_chr = NULL, tfmn_1L_chr = "NTF", predr_var_nm_1L_chr, 
    covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF", 
    mdl_types_lup = NULL, predn_type_1L_chr = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- data_tb %>% dplyr::filter(!is.na(!!rlang::sym(predr_var_nm_1L_chr)))
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    mdl_desc_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup, 
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
        target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F)
    folds_ls <- make_folds_ls(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        folds_1L_int = folds_1L_int)
    smry_of_one_predr_mdl_tb <- purrr::map_dfr(folds_ls, ~{
        model_mdl <- make_mdl(data_tb[-.x, ], depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
            start_1L_chr = start_1L_chr, tfmn_1L_chr = tfmn_1L_chr, 
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr, 
            mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup)
        pred_old_dbl <- stats::predict(model_mdl, type = predn_type_1L_chr)
        pred_new_dbl <- stats::predict(model_mdl, newdata = data_tb[.x, 
            ], type = predn_type_1L_chr) %>% calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, 
            tfmn_is_outp_1L_lgl = T)
        tibble::tibble(Rsquared = caret::R2(pred_old_dbl, data_tb[-.x, 
            ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)), 
            form = "traditional"), RMSE = caret::RMSE(pred_old_dbl, 
            data_tb[-.x, ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))), 
            MAE = caret::MAE(pred_old_dbl, data_tb[-.x, ] %>% 
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))), 
            RsquaredP = caret::R2(pred_new_dbl, data_tb[.x, ] %>% 
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)), 
                form = "traditional"), RMSEP = caret::RMSE(pred_new_dbl, 
                data_tb[.x, ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))), 
            MAEP = caret::MAE(pred_new_dbl, data_tb[.x, ] %>% 
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))))
    }) %>% dplyr::summarise_all(mean) %>% dplyr::mutate(Model = mdl_desc_1L_chr) %>% 
        dplyr::select(Model, dplyr::everything())
    return(smry_of_one_predr_mdl_tb)
}
#' Make summary of time series model output
#' @description make_smry_of_ts_mdl_outp() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make summary of time series model output. The function returns Summary of time series (a list of models).
#' @param data_tb Data (a tibble)
#' @param fn Function (a function)
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one), Default: 'NA'
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one), Default: 'Baseline'
#' @param predictors_lup Predictors (a lookup table)
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iterations (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Summary of time series (a list of models)
#' @rdname make_smry_of_ts_mdl_outp
#' @export 
#' @importFrom purrr map_dbl
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom rlang exec
#' @keywords internal
make_smry_of_ts_mdl_outp <- function (data_tb, fn, predr_vars_nms_chr, mdl_nm_1L_chr, path_to_write_to_1L_chr = NA_character_, 
    depnt_var_nm_1L_chr = "utl_total_w", id_var_nm_1L_chr = "fkClientID", 
    round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", 
    predictors_lup, backend_1L_chr = getOption("brms.backend", 
        "rstan"), iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    scaling_fctr_dbl <- predr_vars_nms_chr %>% purrr::map_dbl(~ifelse(.x %in% 
        predictors_lup$short_name_chr, ready4fun::get_from_lup_obj(predictors_lup, 
        target_var_nm_1L_chr = "mdl_scaling_dbl", match_value_xx = .x, 
        match_var_nm_1L_chr = "short_name_chr", evaluate_lgl = F), 
        1))
    tfd_data_tb <- transform_tb_to_mdl_inp(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr, 
        scaling_fctr_dbl = scaling_fctr_dbl)
    tfd_depnt_var_nm_1L_chr <- ifelse(identical(fn, fit_clg_log_tfmn), 
        transform_depnt_var_nm_for_cll(depnt_var_nm_1L_chr), 
        depnt_var_nm_1L_chr)
    args_ls <- list(data_tb = tfd_data_tb, depnt_var_nm_1L_chr = tfd_depnt_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, iters_1L_int = iters_1L_int, 
        backend_1L_chr = backend_1L_chr, seed_1L_int = seed_1L_int)
    mdl_ls <- rlang::exec(fn, !!!args_ls)
    smry_of_ts_mdl_ls <- list(smry_of_ts_mdl_tb = make_smry_of_brm_mdl(mdl_ls, 
        data_tb = tfd_data_tb, depnt_var_nm_1L_chr = tfd_depnt_var_nm_1L_chr, 
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
            tfd_data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
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
#' Make transformation comparison
#' @description make_tfmn_cmprsn_plt() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make transformation comparison plot. The function returns Transformation comparison (a plot).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one)
#' @param dictionary_tb Dictionary (a tibble)
#' @return Transformation comparison (a plot)
#' @rdname make_tfmn_cmprsn_plt
#' @export 
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @importFrom psych logit
#' @importFrom ggplot2 ggplot aes geom_rug facet_wrap theme_bw labs
#' @importFrom ggalt geom_bkde
#' @importFrom viridis scale_fill_viridis
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
make_tfmn_cmprsn_plt <- function (data_tb, depnt_var_nm_1L_chr, dictionary_tb) 
{
    tfmn_cmprsn_plt <- tidyr::gather(data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(paste0(depnt_var_nm_1L_chr, 
        "_log")), log(!!rlang::sym(depnt_var_nm_1L_chr))), `:=`(!!rlang::sym(paste0(depnt_var_nm_1L_chr, 
        "_logit")), psych::logit(!!rlang::sym(depnt_var_nm_1L_chr))), 
        `:=`(!!rlang::sym(paste0(depnt_var_nm_1L_chr, "_loglog")), 
            -log(-log(!!rlang::sym(depnt_var_nm_1L_chr)))), `:=`(!!rlang::sym(paste0(depnt_var_nm_1L_chr, 
            "_cloglog")), log(-log(1 - !!rlang::sym(depnt_var_nm_1L_chr))))), 
        variable, value, !!rlang::sym(depnt_var_nm_1L_chr), !!rlang::sym(paste0(depnt_var_nm_1L_chr, 
            "_log")), !!rlang::sym(paste0(depnt_var_nm_1L_chr, 
            "_logit")), !!rlang::sym(paste0(depnt_var_nm_1L_chr, 
            "_loglog")), !!rlang::sym(paste0(depnt_var_nm_1L_chr, 
            "_cloglog"))) %>% dplyr::mutate(variable = factor(variable, 
        levels = paste0(depnt_var_nm_1L_chr, c("", "_log", "_logit", 
            "_loglog", "_cloglog")), labels = c("No transformation", 
            "Log", "Logit", "Log-log", "Complementary log-log"))) %>% 
        ggplot2::ggplot(ggplot2::aes(x = value, fill = variable)) + 
        ggalt::geom_bkde() + ggplot2::geom_rug() + viridis::scale_fill_viridis(guide = FALSE, 
        discrete = TRUE) + ggplot2::facet_wrap(~variable, scales = "free") + 
        ggplot2::theme_bw() + ggplot2::labs(x = paste0("Transformed ", 
        dictionary_tb %>% ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "var_nm_chr", 
            match_value_xx = depnt_var_nm_1L_chr, target_var_nm_1L_chr = "var_desc_chr", 
            evaluate_lgl = F)))
    return(tfmn_cmprsn_plt)
}
#' Make unique list element index
#' @description make_unique_ls_elmt_idx_int() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make unique list element index integer vector. The function returns Unique list element index (an integer vector).
#' @param data_ls Data (a list)
#' @return Unique list element index (an integer vector)
#' @rdname make_unique_ls_elmt_idx_int
#' @export 
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate case_when group_by row_number ungroup
#' @importFrom purrr map2_chr map flatten_int
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
make_unique_ls_elmt_idx_int <- function (data_ls) 
{
    combos_tb <- tibble::as_tibble(data_ls, .name_repair = ~paste0("r_", 
        1:length(data_ls))) %>% t() %>% as.data.frame()
    combos_tb <- combos_tb %>% tibble::as_tibble()
    combos_tb <- combos_tb %>% dplyr::mutate(V2 = dplyr::case_when(V1 == 
        V2 ~ NA_character_, T ~ V2)) %>% dplyr::mutate(combo_chr = purrr::map2_chr(V1, 
        V2, ~ifelse(ncol(combos_tb) == 1 | is.na(.y), .x, paste0(.x, 
            "_", .y))))
    combos_tb <- combos_tb %>% dplyr::group_by(combo_chr) %>% 
        dplyr::mutate(combo_id = dplyr::row_number())
    unique_ls_elmt_idx_int <- purrr::map(data_ls %>% unique(), 
        ~ready4fun::get_from_lup_obj(combos_tb %>% dplyr::ungroup(), 
            match_var_nm_1L_chr = "combo_chr", match_value_xx = paste0(.x[1], 
                ifelse(is.na(.x[2]), "", paste0("_", .x[2]))), 
            target_var_nm_1L_chr = "combo_id", evaluate_lgl = F)) %>% 
        purrr::flatten_int()
    return(unique_ls_elmt_idx_int)
}
