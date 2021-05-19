#' Make analysis dataset summary
#' @description make_analysis_ds_smry_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make analysis dataset summary list. The function returns Analysis dataset summary (a list).
#' @param ds_descvs_ls Dataset descriptives (a list)
#' @param candidate_covar_nms_chr Candidate covariate names (a character vector)
#' @param predictors_lup Predictors (a lookup table)
#' @return Analysis dataset summary (a list)
#' @rdname make_analysis_ds_smry_ls
#' @export 

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
#' Make cmpst scatter and dnsty
#' @description make_cmpst_sctr_and_dnsty_plt() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make cmpst scatter and dnsty plot. The function is called for its side effects and does not return a value.
#' @param outp_smry_ls Output summary (a list)
#' @param output_data_dir_1L_chr Output data directory (a character vector of length one)
#' @param predr_var_nms_chr Predictor variable names (a character vector)
#' @param labels_chr Labels (a character vector), Default: c("A", "B", "C", "D")
#' @param label_x_1L_dbl Label x (a double vector of length one), Default: 0.1
#' @param label_y_1L_dbl Label y (a double vector of length one), Default: 0.9
#' @param label_size_1L_dbl Label size (a double vector of length one), Default: 22
#' @return NULL
#' @rdname make_cmpst_sctr_and_dnsty_plt
#' @export 
#' @importFrom purrr map_lgl map
#' @importFrom stringr str_detect
#' @importFrom cowplot ggdraw draw_image plot_grid
#' @keywords internal
make_cmpst_sctr_and_dnsty_plt <- function (outp_smry_ls, output_data_dir_1L_chr, predr_var_nms_chr, 
    labels_chr = c("A", "B", "C", "D"), label_x_1L_dbl = 0.1, 
    label_y_1L_dbl = 0.9, label_size_1L_dbl = 22) 
{
    plot_ls <- paste0(output_data_dir_1L_chr, "/", outp_smry_ls$file_paths_chr[outp_smry_ls$file_paths_chr %>% 
        purrr::map_lgl(~stringr::str_detect(.x, paste0(predr_var_nms_chr, 
            "_1")) & (stringr::str_detect(.x, "_dnst.png") | 
            stringr::str_detect(.x, "_sctr_plt.png")))]) %>% 
        purrr::map(~cowplot::ggdraw() + cowplot::draw_image(.x))
    composite_plt <- cowplot::plot_grid(plot_ls[[1]], plot_ls[[2]], 
        plot_ls[[3]], plot_ls[[4]], nrow = 2, labels = labels_chr, 
        label_x = label_x_1L_dbl, label_y = label_y_1L_dbl, label_size = label_size_1L_dbl)
}
#' Make cohort
#' @description make_cohort_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make cohort list. The function returns Cohort (a list).
#' @param descv_tbls_ls Descriptive tables (a list)
#' @param ctgl_vars_regrouping_ls Ctgl variables regrouping (a list), Default: NULL
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @return Cohort (a list)
#' @rdname make_cohort_ls
#' @export 
#' @importFrom dplyr filter pull
#' @importFrom purrr map_dbl map map2
#' @importFrom stringr str_remove
#' @importFrom stats setNames
#' @keywords internal
make_cohort_ls <- function (descv_tbls_ls, ctgl_vars_regrouping_ls = NULL, nbr_of_digits_1L_int = 2L) 
{
    numeric_vars_chr <- descv_tbls_ls$cohort_desc_tb %>% dplyr::filter(label == 
        "Median (Q1, Q3)") %>% dplyr::pull(variable)
    ctgl_vars_chr <- unique(descv_tbls_ls$cohort_desc_tb$variable[!descv_tbls_ls$cohort_desc_tb$variable %in% 
        numeric_vars_chr])
    nbr_by_round_dbl <- c("Baseline_val_1_dbl", "Follow-up_val_1_dbl") %>% 
        purrr::map_dbl(~descv_tbls_ls$cohort_desc_tb %>% dplyr::filter(variable == 
            ctgl_vars_chr[1]) %>% dplyr::pull(.x) %>% as.numeric() %>% 
            purrr::map_dbl(~.x[[1]]) %>% sum())
    numeric_vars_smry_ls <- numeric_vars_chr %>% purrr::map(~{
        var_smry_tb <- descv_tbls_ls$cohort_desc_tb %>% dplyr::filter(variable == 
            .x)
        list(bl_min_1L_dbl = var_smry_tb %>% dplyr::filter(label == 
            "Min - Max") %>% dplyr::pull(Baseline_val_1_dbl) %>% 
            as.numeric(), bl_max_1L_dbl = var_smry_tb %>% dplyr::filter(label == 
            "Min - Max") %>% dplyr::pull(Baseline_val_2_ls) %>% 
            as.numeric(), bl_mean_1L_dbl = round(var_smry_tb %>% 
            dplyr::filter(label == "Mean (SD)") %>% dplyr::pull(Baseline_val_1_dbl) %>% 
            as.numeric(), nbr_of_digits_1L_int), bl_sd_1L_dbl = round(var_smry_tb %>% 
            dplyr::filter(label == "Mean (SD)") %>% dplyr::pull(Baseline_val_2_ls) %>% 
            stringr::str_remove("\\(") %>% stringr::str_remove("\\)") %>% 
            as.numeric(), nbr_of_digits_1L_int))
    }) %>% stats::setNames(numeric_vars_chr)
    cohort_ls <- list(n_all_1l_dbl = descv_tbls_ls$ds_descvs_ls$nbr_participants_1L_int, 
        n_inc_1L_dbl = nbr_by_round_dbl[1], n_fup_1L_dbl = nbr_by_round_dbl[2], 
        numeric_vars_smry_ls = numeric_vars_smry_ls)
    if (!is.null(ctgl_vars_regrouping_ls)) {
        append_ls <- ctgl_vars_regrouping_ls %>% purrr::map2(names(ctgl_vars_regrouping_ls), 
            ~{
                var_nm_1L_chr <- .y
                .x %>% purrr::map(~{
                  list(name_1L_chr = .x$name_1L_chr, n_in_group_1L_dbl = descv_tbls_ls$cohort_desc_tb %>% 
                    dplyr::filter(variable == var_nm_1L_chr) %>% 
                    dplyr::filter(label %in% .x$ctgs_chr) %>% 
                    dplyr::pull(Baseline_val_1_dbl) %>% as.numeric() %>% 
                    purrr::map_dbl(~.x) %>% sum(), n_msng_1L_dbl = descv_tbls_ls$cohort_desc_tb %>% 
                    dplyr::filter(variable == var_nm_1L_chr) %>% 
                    dplyr::filter(label == "Missing") %>% dplyr::pull(Baseline_val_1_dbl) %>% 
                    as.numeric())
                })
            })
        cohort_ls <- append(cohort_ls, append_ls)
    }
    return(cohort_ls)
}
#' Make eq5d dataset dictionary
#' @description make_eq5d_ds_dict() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make eq5d dataset dictionary. The function returns Dictionary (a tibble).
#' @param data_tb Data (a tibble), Default: make_fake_eq5d_ds()
#' @param predictors_lup Predictors (a lookup table), Default: make_psych_predrs_lup()
#' @return Dictionary (a tibble)
#' @rdname make_eq5d_ds_dict
#' @export 
#' @importFrom youthvars make_tfd_repln_ds_dict_r3 make_final_rpln_ds_dict
#' @importFrom dplyr filter arrange
#' @importFrom ready4use make_pt_ready4_dictionary
#' @keywords internal
make_eq5d_ds_dict <- function (data_tb = make_fake_eq5d_ds(), predictors_lup = make_psych_predrs_lup()) 
{
    dictionary_tb <- youthvars::make_tfd_repln_ds_dict_r3() %>% 
        dplyr::filter(var_nm_chr %in% names(data_tb)) %>% youthvars::make_final_rpln_ds_dict(additions_tb = ready4use::make_pt_ready4_dictionary(var_nm_chr = c("uid", 
        "Timepoint", "bl_data_collection_dtm", paste0("eq5dq_", 
            c("MO", "SC", "UA", "PD", "AD")), "EQ5D_total_dbl", 
        predictors_lup$short_name_chr), var_ctg_chr = c("Identifier", 
        rep("Temporal", 2), rep("Multi-Attribute Utility Instrument Question", 
            5), "Multi-Attribute Utility Instrument Score", rep("Clinical", 
            2)), var_desc_chr = c("Unique identifier", "Data collection round", 
        "Date of baseline data collection", "EQ5D - Mobility Domain Score", 
        "EQ5D - Self-Care Domain Score", "EQ5D - Usual Activities Domain Score", 
        "EQ5D - Pain / Discomfort Domain Score", "EQ5D - Anxiety / Depression Domain Score", 
        "EQ5D - Total weighted score", predictors_lup$long_name_chr), 
        var_type_chr = c("integer", "character", "date", rep("integer", 
            5), "double", predictors_lup$class_chr))) %>% dplyr::arrange(var_ctg_chr)
    return(dictionary_tb)
}
#' Make fake eq5d dataset
#' @description make_fake_eq5d_ds() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make fake eq5d dataset. The function returns Data (a tibble).
#' @param fl_nm_1L_chr File name (a character vector of length one), Default: 'eq5d5l_example.xlsx'
#' @param country_1L_chr Country (a character vector of length one), Default: 'UK'
#' @param version_1L_chr Version (a character vector of length one), Default: '5L'
#' @param type_1L_chr Type (a character vector of length one), Default: 'CW'
#' @return Data (a tibble)
#' @rdname make_fake_eq5d_ds
#' @export 
#' @importFrom readxl read_excel
#' @importFrom dplyr mutate left_join n group_by ungroup arrange filter pull select rename case_when across everything rename_with
#' @importFrom eq5d eq5d
#' @importFrom tibble tibble
#' @importFrom sn rsn
#' @importFrom purrr map_dbl map2_dbl map_int
#' @importFrom faux rnorm_pre
#' @importFrom youthvars transform_raw_ds_for_analysis
#' @importFrom stringr str_c
make_fake_eq5d_ds <- function (fl_nm_1L_chr = "eq5d5l_example.xlsx", country_1L_chr = "UK", 
    version_1L_chr = "5L", type_1L_chr = "CW") 
{
    require(eq5d)
    data_tb <- readxl::read_excel(system.file("extdata", fl_nm_1L_chr, 
        package = "eq5d")) %>% dplyr::mutate(total_eq5d = eq5d::eq5d(., 
        country = country_1L_chr, version = version_1L_chr, type = type_1L_chr))
    k10_lup_tb <- tibble::tibble(k10_dbl = c(sn::rsn(2500, 10.3, 
        omega = 0.1) %>% round(), sn::rsn(2500, 12, omega = 0.4) %>% 
        round(), sn::rsn(2500, 14.5, omega = 0.4) %>% round(), 
        sn::rsn(2500, 21, omega = 6, alpha = 1) %>% round() %>% 
            purrr::map_dbl(~max(.x, 10) %>% min(50))) %>% sample(10000), 
        pred_eq5d_dbl = purrr::map2_dbl(k10_dbl, rnorm(10000, 
            0, 0.075), ~predict_utl_from_k10(.x, eq5d_error_1L_dbl = .y)[2]), 
        match_idx_int = purrr::map_dbl(pred_eq5d_dbl, ~which.min(abs(data_tb$total_eq5d - 
            .x))))
    data_tb <- dplyr::left_join(k10_lup_tb, data_tb %>% dplyr::mutate(match_idx_int = 1:dplyr::n()))
    data_tb <- data_tb %>% dplyr::group_by(Group) %>% dplyr::mutate(uid = 1:dplyr::n()) %>% 
        dplyr::ungroup() %>% dplyr::arrange(uid, Group)
    bl_uids_chr <- data_tb %>% dplyr::filter(Group == "Group1") %>% 
        dplyr::pull(uid)
    data_tb <- data_tb %>% dplyr::filter(uid %in% bl_uids_chr)
    data_tb <- data_tb %>% dplyr::mutate(psych_well_int = faux::rnorm_pre(k10_dbl, 
        mu = 69.9, sd = 9.9, r = -0.56) %>% round(0) %>% purrr::map_int(~min(.x, 
        90) %>% max(18) %>% as.integer()))
    data_tb <- data_tb %>% dplyr::select(uid, Group, MO, SC, 
        UA, PD, AD, k10_dbl, psych_well_int) %>% dplyr::rename(Timepoint = Group) %>% 
        dplyr::mutate(Timepoint = dplyr::case_when(Timepoint == 
            "Group1" ~ "BL", T ~ "FUP")) %>% dplyr::mutate(dplyr::across(where(is.numeric), 
        ~as.integer(.x))) %>% dplyr::rename(k10_int = k10_dbl)
    data("replication_popl_tb", package = "youthvars", envir = environment())
    demog_data_tb <- replication_popl_tb %>% youthvars::transform_raw_ds_for_analysis() %>% 
        dplyr::filter(round == "Baseline") %>% dplyr::mutate(uid = 1:dplyr::n()) %>% 
        dplyr::select(uid, d_interview_date, d_age, Gender, d_sex_birth_s, 
            d_sexual_ori_s, d_relation_s, d_ATSI, CALD, Region, 
            d_studying_working) %>% dplyr::rename(bl_data_collection_dtm = d_interview_date)
    data_tb <- dplyr::left_join(data_tb %>% dplyr::filter(uid %in% 
        demog_data_tb$uid), demog_data_tb) %>% dplyr::select(uid, 
        Timepoint, bl_data_collection_dtm, d_age, Gender, d_sex_birth_s, 
        d_sexual_ori_s, d_relation_s, d_ATSI, CALD, Region, d_studying_working, 
        dplyr::everything()) %>% dplyr::rename_with(~stringr::str_c("eq5dq_", 
        .), .cols = c("MO", "SC", "UA", "PD", "AD"))
    return(data_tb)
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
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param folds_1L_int Folds (an integer vector of length one), Default: 10
#' @return Folds (a list)
#' @rdname make_folds_ls
#' @export 
#' @importFrom caret createFolds
#' @importFrom dplyr pull
#' @importFrom rlang sym
#' @keywords internal
make_folds_ls <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", folds_1L_int = 10L) 
{
    folds_ls <- caret::createFolds(data_tb %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)), 
        k = folds_1L_int, list = TRUE, returnTrain = FALSE)
    return(folds_ls)
}
#' Make health utility and predictors
#' @description make_hlth_utl_and_predrs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make health utility and predictors list. The function returns Health utility and predictors (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param descv_tbls_ls Descriptive tables (a list)
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @param old_nms_chr Old names (a character vector), Default: NULL
#' @param new_nms_chr New names (a character vector), Default: NULL
#' @return Health utility and predictors (a list)
#' @rdname make_hlth_utl_and_predrs_ls
#' @export 
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom dplyr filter
#' @importFrom stringr str_remove
#' @keywords internal
make_hlth_utl_and_predrs_ls <- function (outp_smry_ls, descv_tbls_ls, nbr_of_digits_1L_int = 2L, 
    old_nms_chr = NULL, new_nms_chr = NULL) 
{
    ranked_predrs_ls <- make_ranked_predrs_ls(descv_tbls_ls, 
        old_nms_chr = old_nms_chr, new_nms_chr = new_nms_chr)
    var_nm_1L_chr <- descv_tbls_ls$ds_descvs_ls$dictionary_tb %>% 
        ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "var_nm_chr", 
            match_value_xx = outp_smry_ls$depnt_var_nm_1L_chr, 
            target_var_nm_1L_chr = "var_desc_chr", evaluate_lgl = F) %>% 
        as.vector()
    hlth_utl_and_predrs_ls = list(bl_hu_mean_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>% 
        dplyr::filter(label == "Mean (SD)") %>% ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable", 
        match_value_xx = var_nm_1L_chr, target_var_nm_1L_chr = "Baseline_val_1_dbl", 
        evaluate_lgl = F) %>% as.numeric() %>% round(nbr_of_digits_1L_int), 
        bl_hu_sd_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>% 
            dplyr::filter(label == "Mean (SD)") %>% ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable", 
            match_value_xx = var_nm_1L_chr, target_var_nm_1L_chr = "Baseline_val_2_ls", 
            evaluate_lgl = F) %>% stringr::str_remove("\\(") %>% 
            stringr::str_remove("\\)") %>% as.numeric() %>% round(nbr_of_digits_1L_int), 
        fup_hu_mean_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>% 
            dplyr::filter(label == "Mean (SD)") %>% ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable", 
            match_value_xx = var_nm_1L_chr, target_var_nm_1L_chr = "Follow-up_val_1_dbl", 
            evaluate_lgl = F) %>% as.numeric() %>% round(nbr_of_digits_1L_int), 
        fup_hu_sd_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>% 
            dplyr::filter(label == "Mean (SD)") %>% ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable", 
            match_value_xx = var_nm_1L_chr, target_var_nm_1L_chr = "Follow-up_val_2_ls", 
            evaluate_lgl = F) %>% stringr::str_remove("\\(") %>% 
            stringr::str_remove("\\)") %>% as.numeric() %>% round(nbr_of_digits_1L_int), 
        predrs_nartv_seq_chr = ranked_predrs_ls$unranked_predrs_chr, 
        cor_seq_dscdng_chr = ranked_predrs_ls$ranked_predrs_chr)
    return(hlth_utl_and_predrs_ls)
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
#' Make model coefficient ratio
#' @description make_mdl_coef_ratio_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make model coefficient ratio list. The function returns Model coefficient ratios (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param predr_ctgs_ls Predictor category categoriess (a list), Default: NULL
#' @return Model coefficient ratios (a list)
#' @rdname make_mdl_coef_ratio_ls
#' @export 
#' @importFrom purrr map map2 map_dbl
#' @importFrom dplyr filter pull
#' @importFrom stats setNames
#' @keywords internal
make_mdl_coef_ratio_ls <- function (outp_smry_ls, predr_ctgs_ls = NULL) 
{
    main_mdls_ls <- outp_smry_ls$predr_cmprsn_tb$predr_chr %>% 
        purrr::map(~paste0(paste0(.x, "_1_"), outp_smry_ls$prefd_mdl_types_chr))
    ratios_ls <- main_mdls_ls %>% purrr::map2(outp_smry_ls$predr_cmprsn_tb$predr_chr, 
        ~{
            mdls_chr <- .x
            predr_1L_chr <- .y
            mdls_chr %>% purrr::map_dbl(~{
                coefs_dbl <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Model %in% 
                  .x) %>% dplyr::filter(Parameter %in% paste0(predr_1L_chr, 
                  (c(" baseline", " change")))) %>% dplyr::pull(Estimate)
                coefs_dbl[2]/coefs_dbl[1]
            })
        })
    mdl_coef_ratios_ls = list(mean_ratios_dbl = ratios_ls %>% 
        purrr::map_dbl(~mean(.x)))
    if (!is.null(predr_ctgs_ls)) {
        append_ls <- purrr::map(predr_ctgs_ls, ~mdl_coef_ratios_ls$mean_ratios_dbl[outp_smry_ls$predr_cmprsn_tb$predr_chr %in% 
            .x]) %>% stats::setNames(names(predr_ctgs_ls))
        mdl_coef_ratios_ls <- append(mdl_coef_ratios_ls, append_ls)
    }
    return(mdl_coef_ratios_ls)
}
#' Make model names
#' @description make_mdl_nms_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make model names list. The function returns Model names (a list).
#' @param predr_vars_nms_ls Predictor variables names (a list)
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
#' Make model type summary table
#' @description make_mdl_type_smry_tbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make model type summary table. The function returns Model type summary table (a tibble).
#' @param mdls_tb Models (a tibble)
#' @param mdl_nms_chr Model names (a character vector)
#' @param mdl_type_1L_chr Model type (a character vector of length one)
#' @param add_mdl_nm_sfx_1L_lgl Add model name suffix (a logical vector of length one), Default: T
#' @return Model type summary table (a tibble)
#' @rdname make_mdl_type_smry_tbl
#' @export 
#' @importFrom purrr map_dfr
#' @keywords internal
make_mdl_type_smry_tbl <- function (mdls_tb, mdl_nms_chr, mdl_type_1L_chr, add_mdl_nm_sfx_1L_lgl = T) 
{
    mdl_type_smry_tbl_tb <- mdl_nms_chr %>% purrr::map_dfr(~make_sngl_mdl_smry_tb(mdls_tb, 
        mdl_nm_1L_chr = .x, mdl_type_1L_chr = mdl_type_1L_chr, 
        add_mdl_nm_sfx_1L_lgl = add_mdl_nm_sfx_1L_lgl))
    return(mdl_type_smry_tbl_tb)
}
#' Make models list
#' @description make_mdls_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make models list. The function returns Models (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param mdls_tb Models (a tibble)
#' @return Models (a list)
#' @rdname make_mdls_ls
#' @export 
#' @importFrom purrr map
#' @keywords internal
make_mdls_ls <- function (outp_smry_ls, mdls_tb) 
{
    mdls_chr <- mdls_tb$Model %>% unique()
    mdls_ls <- outp_smry_ls$prefd_mdl_types_chr %>% purrr::map(~mdls_chr[mdls_chr %>% 
        endsWith(.x)])
    return(mdls_ls)
}
#' Make models summary tables list
#' @description make_mdls_smry_tbls_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make models summary tables list. The function returns Models summary tables (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @return Models summary tables (a list)
#' @rdname make_mdls_smry_tbls_ls
#' @export 
#' @importFrom dplyr mutate across filter
#' @importFrom purrr map flatten_chr map_lgl
#' @keywords internal
make_mdls_smry_tbls_ls <- function (outp_smry_ls, nbr_of_digits_1L_int = 2L) 
{
    mdls_smry_tb <- outp_smry_ls$mdls_smry_tb %>% dplyr::mutate(dplyr::across(c("Estimate", 
        "SE"), ~round(.x, nbr_of_digits_1L_int) %>% format(nsmall = nbr_of_digits_1L_int))) %>% 
        dplyr::mutate(`95% CI` = `95% CI` %>% transform_chr_digit_pairs())
    rownames(mdls_smry_tb) <- NULL
    indpt_predrs_mdls_tb <- mdls_smry_tb %>% dplyr::filter(Model %in% 
        (paste0(outp_smry_ls$predr_cmprsn_tb$predr_chr, "_1_") %>% 
            purrr::map(~paste0(.x, outp_smry_ls$prefd_mdl_types_chr)) %>% 
            purrr::flatten_chr()))
    covar_mdls_tb <- mdls_smry_tb %>% dplyr::filter(!Model %in% 
        indpt_predrs_mdls_tb$Model)
    prefd_predr_mdl_smry_tb <- indpt_predrs_mdls_tb %>% dplyr::filter(Model %>% 
        purrr::map_lgl(~startsWith(.x, outp_smry_ls$predr_vars_nms_ls[[1]])))
    mdls_smry_tbls_ls <- list(indpt_predrs_mdls_tb = indpt_predrs_mdls_tb, 
        covar_mdls_tb = covar_mdls_tb, prefd_predr_mdl_smry_tb = prefd_predr_mdl_smry_tb)
    return(mdls_smry_tbls_ls)
}
#' Make paths to ss plots
#' @description make_paths_to_ss_plts_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make paths to ss plots list. The function returns Paths to ss plots (a list).
#' @param output_data_dir_1L_chr Output data directory (a character vector of length one)
#' @param outp_smry_ls Output summary (a list)
#' @param additional_paths_chr Additional paths (a character vector), Default: '/dens_and_sctr.png'
#' @return Paths to ss plots (a list)
#' @rdname make_paths_to_ss_plts_ls
#' @export 
#' @importFrom purrr map_lgl
#' @importFrom stringr str_detect
#' @keywords internal
make_paths_to_ss_plts_ls <- function (output_data_dir_1L_chr, outp_smry_ls, additional_paths_chr = "/dens_and_sctr.png") 
{
    paths_to_ss_plts_ls = list(combined_utl = paste0(output_data_dir_1L_chr, 
        "/_Descriptives/combined_utl.png"), composite = paste0(output_data_dir_1L_chr, 
        additional_paths_chr[1]), items = paste0(output_data_dir_1L_chr, 
        "/_Descriptives/qstn_rspns.png"), density = paste0(output_data_dir_1L_chr, 
        "/", outp_smry_ls$file_paths_chr[outp_smry_ls$file_paths_chr %>% 
            purrr::map_lgl(~stringr::str_detect(.x, "A_TFMN_CMPRSN_DNSTY"))]), 
        importance = paste0(output_data_dir_1L_chr, "/", outp_smry_ls$file_paths_chr[outp_smry_ls$file_paths_chr %>% 
            purrr::map_lgl(~stringr::str_detect(.x, "B_PRED_CMPRSN_BORUTA_VAR_IMP"))]))
    return(paths_to_ss_plts_ls)
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
#' @param existing_predrs_ls Existing predictors (a list), Default: NULL
#' @return Predictor variables names (a list)
#' @rdname make_predr_vars_nms_ls
#' @export 
#' @importFrom purrr map discard flatten map_lgl
make_predr_vars_nms_ls <- function (main_predrs_chr, covars_ls, existing_predrs_ls = NULL) 
{
    predr_vars_nms_ls <- covars_ls %>% purrr::map(~{
        covars_chr <- .x
        purrr::map(main_predrs_chr, ~list(c(.x), c(.x, covars_chr) %>% 
            purrr::discard(is.na))) %>% purrr::flatten()
    }) %>% purrr::flatten() %>% unique()
    predr_vars_nms_ls <- predr_vars_nms_ls[order(sapply(predr_vars_nms_ls, 
        length))]
    if (!is.null(existing_predrs_ls)) {
        predr_vars_nms_ls <- predr_vars_nms_ls[predr_vars_nms_ls %>% 
            purrr::map_lgl(~{
                test_chr <- .x
                !any(existing_predrs_ls %>% purrr::map_lgl(~identical(.x, 
                  test_chr)))
            })]
    }
    return(predr_vars_nms_ls)
}
#' Make predictors for best models
#' @description make_predrs_for_best_mdls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make predictors for best models. The function returns Predictors for best models (a character vector).
#' @param outp_smry_ls Output summary (a list)
#' @param old_nms_chr Old names (a character vector), Default: NULL
#' @param new_nms_chr New names (a character vector)
#' @return Predictors for best models (a character vector)
#' @rdname make_predrs_for_best_mdls
#' @export 
#' @importFrom dplyr filter arrange desc pull
#' @importFrom purrr map flatten_chr map_lgl map2_chr
#' @importFrom stringr str_remove
#' @keywords internal
make_predrs_for_best_mdls <- function (outp_smry_ls, old_nms_chr = NULL, new_nms_chr) 
{
    ordered_mdl_nms_chr <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Parameter == 
        "R2") %>% dplyr::arrange(dplyr::desc(Estimate)) %>% dplyr::pull(Model)
    ind_predr_mdls_by_mdl_type_ls <- outp_smry_ls$prefd_mdl_types_chr %>% 
        purrr::map(~{
            mdl_type_1L_chr <- .x
            paste0(outp_smry_ls$predr_cmprsn_tb$predr_chr, "_1_") %>% 
                purrr::map(~paste0(.x, mdl_type_1L_chr)) %>% 
                purrr::flatten_chr()
        })
    ordered_mdls_by_type_ls <- ind_predr_mdls_by_mdl_type_ls %>% 
        purrr::map(~{
            ind_predr_mdls_chr <- .x
            ordered_mdl_nms_chr[ordered_mdl_nms_chr %>% purrr::map_lgl(~.x %in% 
                ind_predr_mdls_chr)]
        })
    predrs_for_best_mdls_chr <- ordered_mdls_by_type_ls %>% purrr::map2_chr(outp_smry_ls$prefd_mdl_types_chr, 
        ~{
            .x[1] %>% stringr::str_remove(paste0("_1_", .y))
        })
    if (!is.null(old_nms_chr)) {
        predrs_for_best_mdls_chr <- transform_predr_nm_part_of_phrases(predrs_for_best_mdls_chr, 
            old_nms_chr = old_nms_chr, new_nms_chr = new_nms_chr)
    }
    return(predrs_for_best_mdls_chr)
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
make_prefd_mdls_vec <- function (smry_of_sngl_predr_mdls_tb, choose_from_pfx_chr = c("GLM", 
    "OLS"), mdl_types_lup = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", envir = environment(), package = "TTU")
    ordered_mdl_types_chr <- dplyr::inner_join(smry_of_sngl_predr_mdls_tb %>% 
        dplyr::select(Model) %>% dplyr::rename(long_name_chr = Model), 
        mdl_types_lup) %>% dplyr::pull(short_name_chr)
    prefd_mdls_chr <- purrr::map_chr(choose_from_pfx_chr, ~ordered_mdl_types_chr[startsWith(ordered_mdl_types_chr, 
        .x)][1])
    return(prefd_mdls_chr)
}
#' Make psych predictors
#' @description make_psych_predrs_lup() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make psych predictors lookup table. The function returns Predictors (a lookup table).

#' @return Predictors (a lookup table)
#' @rdname make_psych_predrs_lup
#' @export 

#' @keywords internal
make_psych_predrs_lup <- function () 
{
    predictors_lup <- TTU_predictors_lup(make_pt_TTU_predictors_lup(short_name_chr = c("k10_int", 
        "psych_well_int"), long_name_chr = c("Kessler Psychological Distress - 10 Item Total Score", 
        "Overall Wellbeing Measure (Winefield et al. 2012)"), 
        min_val_dbl = c(10, 18), max_val_dbl = c(50, 90), class_chr = "integer", 
        increment_dbl = 1, class_fn_chr = "integer", mdl_scaling_dbl = 0.01, 
        covariate_lgl = F))
    return(predictors_lup)
}
#' Make ranked predictors
#' @description make_ranked_predrs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make ranked predictors list. The function returns Ranked predictors (a list).
#' @param descv_tbls_ls Descriptive tables (a list)
#' @param old_nms_chr Old names (a character vector), Default: NULL
#' @param new_nms_chr New names (a character vector), Default: NULL
#' @return Ranked predictors (a list)
#' @rdname make_ranked_predrs_ls
#' @export 
#' @importFrom purrr map_dbl map flatten_chr
#' @keywords internal
make_ranked_predrs_ls <- function (descv_tbls_ls, old_nms_chr = NULL, new_nms_chr = NULL) 
{
    unranked_predrs_chr <- rownames(descv_tbls_ls[["bl_cors_tb"]])[-1] %>% 
        transform_predr_nm_part_of_phrases(old_nms_chr = old_nms_chr, 
            new_nms_chr = new_nms_chr)
    ranks_dbl <- descv_tbls_ls[["bl_cors_tb"]][2:nrow(descv_tbls_ls[["bl_cors_tb"]]), 
        1] %>% purrr::map_dbl(~{
        vec_dbl <- regmatches(.x, gregexpr("[[:digit:]]+", .x)) %>% 
            unlist() %>% as.numeric()
        vec_dbl[1] + vec_dbl[2]/100
    }) %>% rank()
    unique_ranks_dbl <- unique(ranks_dbl) %>% sort(decreasing = T)
    ranked_predrs_chr <- purrr::map(unique_ranks_dbl, ~unranked_predrs_chr[ranks_dbl == 
        .x]) %>% purrr::flatten_chr()
    ranked_predrs_ls <- list(unranked_predrs_chr = unranked_predrs_chr, 
        ranked_predrs_chr = ranked_predrs_chr)
    return(ranked_predrs_ls)
}
#' Make results
#' @description make_results_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make results list. The function returns Results (a list).
#' @param spine_of_results_ls Spine of results (a list)
#' @param cs_ts_ratios_tb Cs time series ratios (a tibble)
#' @param ctgl_vars_regrouping_ls Ctgl variables regrouping (a list)
#' @param sig_covars_some_predrs_mdls_tb Sig covariates some predictors models (a tibble)
#' @param sig_thresh_covars_1L_chr Sig thresh covariates (a character vector of length one)
#' @return Results (a list)
#' @rdname make_results_ls
#' @export 
#' @importFrom purrr map_lgl
#' @importFrom stringr str_detect
#' @importFrom cowplot save_plot
#' @importFrom tibble tibble
#' @importFrom dplyr filter pull
make_results_ls <- function (spine_of_results_ls, cs_ts_ratios_tb, ctgl_vars_regrouping_ls, 
    sig_covars_some_predrs_mdls_tb, sig_thresh_covars_1L_chr) 
{
    mdls_smry_tbls_ls <- make_mdls_smry_tbls_ls(spine_of_results_ls$outp_smry_ls, 
        nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int)
    covars_mdls_ls <- make_mdls_ls(spine_of_results_ls$outp_smry_ls, 
        mdls_tb = mdls_smry_tbls_ls$covar_mdls_tb)
    descv_tbls_ls <- paste0(spine_of_results_ls$output_data_dir_1L_chr, 
        "/", spine_of_results_ls$outp_smry_ls$file_paths_chr[spine_of_results_ls$outp_smry_ls$file_paths_chr %>% 
            purrr::map_lgl(~stringr::str_detect(.x, "descv_tbls_ls.RDS"))]) %>% 
        readRDS()
    composite_plt <- make_cmpst_sctr_and_dnsty_plt(spine_of_results_ls$outp_smry_ls, 
        output_data_dir_1L_chr = spine_of_results_ls$output_data_dir_1L_chr, 
        predr_var_nms_chr = spine_of_results_ls$outp_smry_ls$predr_vars_nms_ls[[1]])
    cowplot::save_plot(paste0(spine_of_results_ls$output_data_dir_1L_chr, 
        "/dens_and_sctr.png"), composite_plt, base_height = 20)
    ttu_cs_ls = make_ttu_cs_ls(spine_of_results_ls$outp_smry_ls, 
        sig_covars_some_predrs_mdls_tb = sig_covars_some_predrs_mdls_tb, 
        sig_thresh_covars_1L_chr = sig_thresh_covars_1L_chr)
    ttu_lngl_ls = list(best_mdls_tb = tibble::tibble(model_type = c("GLMM", 
        "LMM"), link_and_tfmn_chr = c("Gaussian distribution with log link", 
        "c-loglog transformed)"), name_chr = make_predrs_for_best_mdls(spine_of_results_ls$outp_smry_ls, 
        old_nms_chr = spine_of_results_ls$var_nm_change_lup$old_nms_chr, 
        new_nms_chr = spine_of_results_ls$var_nm_change_lup$new_nms_chr), 
        r2_dbl = mdls_smry_tbls_ls$prefd_predr_mdl_smry_tb %>% 
            dplyr::filter(Parameter == "R2") %>% dplyr::pull(Estimate)), 
        cs_ts_ratios_tb = cs_ts_ratios_tb, incld_covars_chr = spine_of_results_ls$outp_smry_ls$prefd_covars_chr)
    results_ls <- list(cohort_ls = make_cohort_ls(descv_tbls_ls, 
        ctgl_vars_regrouping_ls = ctgl_vars_regrouping_ls, nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int), 
        hlth_utl_and_predrs_ls = make_hlth_utl_and_predrs_ls(spine_of_results_ls$outp_smry_ls, 
            descv_tbls_ls = descv_tbls_ls, nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int, 
            old_nms_chr = spine_of_results_ls$var_nm_change_lup$old_nms_chr, 
            new_nms_chr = spine_of_results_ls$var_nm_change_lup$new_nms_chr), 
        mdl_coef_ratios_ls = spine_of_results_ls$mdl_coef_ratios_ls, 
        paths_to_figs_ls = make_paths_to_ss_plts_ls(spine_of_results_ls$output_data_dir_1L_chr, 
            outp_smry_ls = spine_of_results_ls$outp_smry_ls), 
        r_version_1L_chr = paste0(spine_of_results_ls$outp_smry_ls$session_data_ls$R.version$major, 
            ".", spine_of_results_ls$outp_smry_ls$session_data_ls$R.version$minor), 
        study_descs_ls = spine_of_results_ls$study_descs_ls, 
        tables_ls = make_ss_tbls_ls(spine_of_results_ls$outp_smry_ls, 
            mdls_smry_tbls_ls = mdls_smry_tbls_ls, covars_mdls_ls = covars_mdls_ls, 
            descv_tbls_ls = descv_tbls_ls, nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int), 
        ttu_cs_ls = ttu_cs_ls, ttu_lngl_ls = ttu_lngl_ls)
    return(results_ls)
}
#' Make results list spine
#' @description make_results_ls_spine() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make results list spine. The function returns Spine of results (a list).
#' @param output_data_dir_1L_chr Output data directory (a character vector of length one)
#' @param var_nm_change_lup Variable name change (a lookup table)
#' @param study_descs_ls Study descriptions (a list)
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @return Spine of results (a list)
#' @rdname make_results_ls_spine
#' @export 

make_results_ls_spine <- function (output_data_dir_1L_chr, var_nm_change_lup, study_descs_ls, 
    nbr_of_digits_1L_int = 2L) 
{
    outp_smry_ls <- readRDS(paste0(output_data_dir_1L_chr, "/I_ALL_OUTPUT_.RDS"))
    mdl_coef_ratios_ls <- make_mdl_coef_ratio_ls(outp_smry_ls, 
        predr_ctgs_ls = study_descs_ls$predr_ctgs_ls)
    mdls_smry_tbls_ls <- make_mdls_smry_tbls_ls(outp_smry_ls, 
        nbr_of_digits_1L_int = nbr_of_digits_1L_int)
    covars_mdls_ls <- make_mdls_ls(outp_smry_ls, mdls_tb = mdls_smry_tbls_ls$covar_mdls_tb)
    spine_of_results_ls <- list(outp_smry_ls = outp_smry_ls, 
        output_data_dir_1L_chr = output_data_dir_1L_chr, mdl_coef_ratios_ls = mdl_coef_ratios_ls, 
        nbr_of_digits_1L_int = nbr_of_digits_1L_int, study_descs_ls = study_descs_ls, 
        var_nm_change_lup = var_nm_change_lup)
    return(spine_of_results_ls)
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
#' @param model_mdl Model (a model), Default: NULL
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
make_smry_of_mdl_outp <- function (data_tb, model_mdl = NULL, folds_1L_int = 10, depnt_var_nm_1L_chr = "utl_total_w", 
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
    control_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup, 
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr, 
        target_var_nm_1L_chr = "control_chr", evaluate_lgl = F)
    smry_of_one_predr_mdl_tb <- purrr::map_dfr(folds_ls, ~{
        model_mdl <- make_mdl(data_tb[-.x, ], depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
            start_1L_chr = start_1L_chr, tfmn_1L_chr = tfmn_1L_chr, 
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr, 
            mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup, 
            control_1L_chr = control_1L_chr)
        pred_old_dbl <- stats::predict(model_mdl, type = predn_type_1L_chr) %>% 
            calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, 
                tfmn_is_outp_1L_lgl = T)
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
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
        iters_1L_int = iters_1L_int, backend_1L_chr = backend_1L_chr, 
        seed_1L_int = seed_1L_int)
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
#' Make single model summary
#' @description make_sngl_mdl_smry_tb() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make single model summary tibble. The function returns New (a tibble).
#' @param mdls_tb Models (a tibble)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param mdl_type_1L_chr Model type (a character vector of length one)
#' @param add_mdl_nm_sfx_1L_lgl Add model name suffix (a logical vector of length one), Default: T
#' @return New (a tibble)
#' @rdname make_sngl_mdl_smry_tb
#' @export 
#' @importFrom dplyr filter select mutate case_when rename_with
#' @importFrom tibble add_case
#' @importFrom rlang sym
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom purrr map_chr
#' @importFrom stringr str_replace
#' @keywords internal
make_sngl_mdl_smry_tb <- function (mdls_tb, mdl_nm_1L_chr, mdl_type_1L_chr, add_mdl_nm_sfx_1L_lgl = T) 
{
    new_tb <- mdls_tb %>% dplyr::filter(Model == mdl_nm_1L_chr) %>% 
        tibble::add_case(Parameter = mdl_nm_1L_chr, .before = 1) %>% 
        dplyr::select(-Model)
    mdl_nm_sfx_1L_chr <- ifelse(add_mdl_nm_sfx_1L_lgl, paste0("_", 
        mdl_type_1L_chr), "")
    new_tb <- new_tb %>% dplyr::mutate(`:=`(!!rlang::sym(paste0("R2", 
        mdl_nm_sfx_1L_chr)), dplyr::case_when(Parameter == mdl_nm_1L_chr ~ 
        ready4fun::get_from_lup_obj(new_tb, match_value_xx = "R2", 
            match_var_nm_1L_chr = "Parameter", target_var_nm_1L_chr = "Estimate", 
            evaluate_lgl = F), T ~ NA_character_)), `:=`(!!rlang::sym(paste0("Sigma", 
        mdl_nm_sfx_1L_chr)), dplyr::case_when(Parameter == mdl_nm_1L_chr ~ 
        ready4fun::get_from_lup_obj(new_tb, match_value_xx = "Sigma", 
            match_var_nm_1L_chr = "Parameter", target_var_nm_1L_chr = "Estimate", 
            evaluate_lgl = F), T ~ NA_character_))) %>% dplyr::filter(!Parameter %in% 
        c("R2", "RMSE", "Sigma")) %>% dplyr::rename_with(~paste0(.x, 
        mdl_nm_sfx_1L_chr), .cols = c("Parameter", "Estimate", 
        "SE")) %>% dplyr::rename_with(~paste0("CI", mdl_nm_sfx_1L_chr), 
        .cols = c("95% CI")) %>% dplyr::mutate(`:=`(!!rlang::sym(paste0("Parameter", 
        mdl_nm_sfx_1L_chr)), !!rlang::sym(paste0("Parameter", 
        mdl_nm_sfx_1L_chr)) %>% purrr::map_chr(~stringr::str_replace(.x, 
        paste0("_1_", mdl_type_1L_chr), " model"))))
    rownames(new_tb) <- NULL
    return(new_tb)
}
#' Make ss tables list
#' @description make_ss_tbls_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make ss tables list. The function returns Ss tables (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param mdls_smry_tbls_ls Models summary tables (a list)
#' @param covars_mdls_ls Covariates models (a list)
#' @param descv_tbls_ls Descriptive tables (a list)
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @return Ss tables (a list)
#' @rdname make_ss_tbls_ls
#' @export 
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate across everything
#' @importFrom purrr map_dbl
#' @importFrom stringr str_replace_all
#' @keywords internal
make_ss_tbls_ls <- function (outp_smry_ls, mdls_smry_tbls_ls, covars_mdls_ls, descv_tbls_ls, 
    nbr_of_digits_1L_int = 2L) 
{
    ss_tbls_ls <- list(mdl_type_1_covar_mdls_tb = make_mdl_type_smry_tbl(mdls_tb = mdls_smry_tbls_ls$covar_mdls_tb, 
        mdl_nms_chr = covars_mdls_ls[[1]], mdl_type_1L_chr = outp_smry_ls$prefd_mdl_types_chr[1], 
        add_mdl_nm_sfx_1L_lgl = F), mdl_type_2_covar_mdls_tb = make_mdl_type_smry_tbl(mdls_tb = mdls_smry_tbls_ls$covar_mdls_tb, 
        mdl_nms_chr = covars_mdls_ls[[2]], mdl_type_1L_chr = outp_smry_ls$prefd_mdl_types_chr[2], 
        add_mdl_nm_sfx_1L_lgl = F), ind_preds_coefs_tbl = make_two_mdl_types_smry_tbl(outp_smry_ls, 
        mdls_tb = mdls_smry_tbls_ls$indpt_predrs_mdls_tb), participant_descs = descv_tbls_ls$cohort_desc_tb, 
        pred_dist_and_cors = descv_tbls_ls$predr_pars_and_cors_tb, 
        tenf_glm = outp_smry_ls[["smry_of_mdl_sngl_predrs_tb"]] %>% 
            tibble::as_tibble() %>% dplyr::mutate(dplyr::across(where(is.numeric), 
            ~.x %>% purrr::map_dbl(~min(max(.x, -1.1), 1.1)))) %>% 
            transform_tbl_to_rnd_vars(nbr_of_digits_1L_int = nbr_of_digits_1L_int) %>% 
            dplyr::mutate(dplyr::across(.cols = dplyr::everything(), 
                ~.x %>% stringr::str_replace_all("-1.10", "< -1.00") %>% 
                  stringr::str_replace_all("1.10", "> 1.00"))), 
        tenf_phq9 = make_tfd_sngl_predr_mdls_tb(outp_smry_ls, 
            nbr_of_digits_1L_int = nbr_of_digits_1L_int))
    return(ss_tbls_ls)
}
#' Make study descriptions
#' @description make_study_descs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make study descriptions list. The function returns Study descriptions (a list).
#' @param health_utl_nm_1L_chr Health utility name (a character vector of length one)
#' @param time_btwn_bl_and_fup_1L_chr Time btwn baseline and follow-up (a character vector of length one)
#' @param predr_ctgs_ls Predictor category categoriess (a list)
#' @return Study descriptions (a list)
#' @rdname make_study_descs_ls
#' @export 

make_study_descs_ls <- function (health_utl_nm_1L_chr, time_btwn_bl_and_fup_1L_chr, 
    predr_ctgs_ls) 
{
    study_descs_ls <- list(health_utl_nm_1L_chr = health_utl_nm_1L_chr, 
        time_btwn_bl_and_fup_1L_chr = time_btwn_bl_and_fup_1L_chr, 
        predr_ctgs_ls = predr_ctgs_ls)
    return(study_descs_ls)
}
#' Make transformed single predictor models
#' @description make_tfd_sngl_predr_mdls_tb() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make transformed single predictor models tibble. The function returns Transformed single predictor models (a tibble).
#' @param outp_smry_ls Output summary (a list)
#' @param nbr_of_digits_1L_int Number of digits (an integer vector of length one), Default: 2
#' @param mdl_pfx_ls Model prefix (a list), Default: list(OLS = "Ordinary Least Squares ", GLM = c("Generalised Linear Mixed Model with ", 
#'    "Beta Regression Model with Binomial "))
#' @return Transformed single predictor models (a tibble)
#' @rdname make_tfd_sngl_predr_mdls_tb
#' @export 
#' @importFrom purrr map2 map_lgl map_dfr
#' @importFrom dplyr filter mutate case_when
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom tibble add_case
#' @keywords internal
make_tfd_sngl_predr_mdls_tb <- function (outp_smry_ls, nbr_of_digits_1L_int = 2L, mdl_pfx_ls = list(OLS = "Ordinary Least Squares ", 
    GLM = c("Generalised Linear Mixed Model with ", "Beta Regression Model with Binomial "))) 
{
    tfd_sngl_predr_mdls_tb <- mdl_pfx_ls %>% purrr::map2(names(mdl_pfx_ls), 
        ~{
            pfx_chr <- .x
            mdls_tb <- outp_smry_ls$smry_of_sngl_predr_mdls_tb %>% 
                dplyr::filter(Model %>% purrr::map_lgl(~{
                  term_1L_chr <- .x
                  pfx_chr %>% purrr::map_lgl(~startsWith(term_1L_chr, 
                    .x)) %>% any()
                }))
            pfx_chr %>% purrr::map_dfr(~{
                pfx_1L_chr <- .x
                mdls_tb %>% dplyr::filter(startsWith(Model, pfx_1L_chr)) %>% 
                  dplyr::mutate(Model = dplyr::case_when(Model %>% 
                    startsWith(mdl_pfx_ls[[2]][2]) ~ stringr::str_replace_all(Model, 
                    pfx_1L_chr, "Beta "), T ~ stringr::str_remove_all(Model, 
                    pfx_1L_chr))) %>% dplyr::mutate(Model = Model %>% 
                  stringr::str_remove_all("\\(|\\)"))
            }) %>% tibble::add_case(Model = .y, .before = 1)
        }) %>% purrr::map_dfr(~.x) %>% transform_tbl_to_rnd_vars(nbr_of_digits_1L_int = nbr_of_digits_1L_int)
    return(tfd_sngl_predr_mdls_tb)
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
#' Make ttu cs
#' @description make_ttu_cs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make ttu cs list. The function returns Ttu cs (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param sig_covars_some_predrs_mdls_tb Sig covariates some predictors models (a tibble)
#' @param sig_thresh_covars_1L_chr Sig thresh covariates (a character vector of length one)
#' @return Ttu cs (a list)
#' @rdname make_ttu_cs_ls
#' @export 

#' @keywords internal
make_ttu_cs_ls <- function (outp_smry_ls, sig_covars_some_predrs_mdls_tb, sig_thresh_covars_1L_chr) 
{
    ttu_cs_ls <- list(best_mdl_types_ls = list(GLM = c("Gaussian distribution and log link"), 
        OLS = c("no transformation", "log transformation", "clog-log transformation")), 
        selected_mdls_chr = c("GLM with Gaussian distribution and log link", 
            "OLS with clog-log transformation"), cs_mdls_predrs_seq_dscdng_chr = outp_smry_ls$smry_of_mdl_sngl_predrs_tb$Predictor, 
        sig_covars_all_predrs_mdls_chr = outp_smry_ls$signt_covars_chr, 
        sig_thresh_covars_1L_chr = sig_thresh_covars_1L_chr, 
        sig_covars_some_predrs_mdls_tb = sig_covars_some_predrs_mdls_tb, 
        rf_seq_dscdng_chr = outp_smry_ls$predr_cmprsn_tb$predr_chr, 
        mdl_predrs_and_rf_seqs_cmprsn_1L_chr = "is consistent")
    return(ttu_cs_ls)
}
#' Make two model types summary table
#' @description make_two_mdl_types_smry_tbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make two model types summary table. The function returns Two model types summary table (a tibble).
#' @param outp_smry_ls Output summary (a list)
#' @param mdls_tb Models (a tibble)
#' @return Two model types summary table (a tibble)
#' @rdname make_two_mdl_types_smry_tbl
#' @export 
#' @importFrom purrr map_dfc
#' @importFrom dplyr select rename_with
#' @keywords internal
make_two_mdl_types_smry_tbl <- function (outp_smry_ls, mdls_tb) 
{
    mdls_ls <- make_mdls_ls(outp_smry_ls, mdls_tb = mdls_tb)
    two_mdl_types_smry_tbl_tb <- 1:2 %>% purrr::map_dfc(~{
        make_mdl_type_smry_tbl(mdls_tb = mdls_tb, mdl_nms_chr = mdls_ls[[.x]], 
            mdl_type_1L_chr = outp_smry_ls$prefd_mdl_types_chr[.x], 
            add_mdl_nm_sfx_1L_lgl = T)
    }) %>% dplyr::select(-paste0("Parameter_", outp_smry_ls$prefd_mdl_types_chr[2])) %>% 
        dplyr::rename_with(~"Parameter", .cols = paste0("Parameter_", 
            outp_smry_ls$prefd_mdl_types_chr[1]))
    return(two_mdl_types_smry_tbl_tb)
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
    if (ncol(combos_tb) > 1) {
        combos_tb <- combos_tb %>% dplyr::mutate(V2 = dplyr::case_when(V1 == 
            V2 ~ NA_character_, T ~ V2)) %>% dplyr::mutate(combo_chr = purrr::map2_chr(V1, 
            V2, ~ifelse(ncol(combos_tb) == 1 | is.na(.y), .x, 
                paste0(.x, "_", .y))))
        combos_tb <- combos_tb %>% dplyr::group_by(combo_chr) %>% 
            dplyr::mutate(combo_id = dplyr::row_number())
        unique_ls_elmt_idx_int <- purrr::map(data_ls %>% unique(), 
            ~ready4fun::get_from_lup_obj(combos_tb %>% dplyr::ungroup(), 
                match_var_nm_1L_chr = "combo_chr", match_value_xx = paste0(.x[1], 
                  ifelse(is.na(.x[2]), "", paste0("_", .x[2]))), 
                target_var_nm_1L_chr = "combo_id", evaluate_lgl = F)) %>% 
            purrr::flatten_int()
    }
    else {
        unique_ls_elmt_idx_int <- 1
    }
    return(unique_ls_elmt_idx_int)
}
