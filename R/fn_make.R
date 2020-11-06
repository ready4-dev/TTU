#' Make adolescent Assessment of Quality of Life Six Dimension disvalue
#' @description make_adol_aqol6d_disv_lup() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make adolescent assessment of quality of life six dimension disvalue lookup table. The function returns Adolescent Assessment of Quality of Life Six Dimension disvalue (a lookup table).

#' @return Adolescent Assessment of Quality of Life Six Dimension disvalue (a lookup table)
#' @rdname make_adol_aqol6d_disv_lup
#' @export
#' @importFrom dplyr mutate case_when
#' @keywords internal
make_adol_aqol6d_disv_lup <- function ()
{
    adol_aqol6d_disv_lup <- aqol6d_adult_disv_lup_tb %>% dplyr::mutate(Answer_4_dbl = dplyr::case_when(Question_chr ==
        "Q18" ~ 0.622, TRUE ~ Answer_4_dbl), Answer_5_dbl = dplyr::case_when(Question_chr ==
        "Q3" ~ 0.827, TRUE ~ Answer_5_dbl), Answer_6_dbl = dplyr::case_when(Question_chr ==
        "Q1" ~ 0.073, TRUE ~ Answer_5_dbl))
    return(adol_aqol6d_disv_lup)
}
#' Make Assessment of Quality of Life Six Dimension adolescent pop tibbles
#' @description make_aqol6d_adol_pop_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make assessment of quality of life six dimension adolescent pop tibbles list. The function returns Assessment of Quality of Life Six Dimension adolescent pop tibbles (a list).
#' @param aqol_items_props_tbs_ls Assessment of Quality of Life items props tibbles (a list)
#' @param aqol_scores_pars_ls Assessment of Quality of Life scores parameters (a list)
#' @param series_names_chr Series names (a character vector)
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param temporal_cors_ls Temporal correlations (a list)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param prefix_chr Prefix (a character vector), Default: c(uid = "Participant_", aqol_item = "aqol6d_q", domain_unwtd_pfx_1L_chr = "aqol6d_subtotal_c_",
#'    domain_wtd_pfx_1L_chr = "aqol6d_subtotal_w_")
#' @return Assessment of Quality of Life Six Dimension adolescent pop tibbles (a list)
#' @rdname make_aqol6d_adol_pop_tbs_ls
#' @export
#' @importFrom purrr map
#' @importFrom dplyr select starts_with everything
#' @importFrom rlang sym
#' @keywords internal
make_aqol6d_adol_pop_tbs_ls <- function (aqol_items_props_tbs_ls, aqol_scores_pars_ls, series_names_chr,
    synth_data_spine_ls, temporal_cors_ls, id_var_nm_1L_chr = "fkClientID",
    prefix_chr = c(uid = "Participant_", aqol_item = "aqol6d_q",
        domain_unwtd_pfx_1L_chr = "aqol6d_subtotal_c_", domain_wtd_pfx_1L_chr = "aqol6d_subtotal_w_"))
{
    item_pfx_1L_chr <- prefix_chr[["aqol_item"]]
    uid_pfx_1L_chr <- prefix_chr[["uid"]]
    aqol6d_adol_pop_tbs_ls <- make_synth_series_tbs_ls(synth_data_spine_ls,
        series_names_chr = series_names_chr) %>% add_cors_and_uts_to_aqol6d_tbs_ls(aqol_scores_pars_ls = aqol_scores_pars_ls,
        aqol_items_props_tbs_ls = aqol_items_props_tbs_ls, temporal_cors_ls = temporal_cors_ls,
        prefix_chr = prefix_chr, aqol_tots_var_nms_chr = synth_data_spine_ls$aqol_tots_var_nms_chr,
        id_var_nm_1L_chr = id_var_nm_1L_chr) %>% purrr::map(~{
        domain_items_ls <- make_domain_items_ls(domain_qs_lup_tb = aqol6d_domain_qs_lup_tb,
            item_pfx_1L_chr = item_pfx_1L_chr)
        domain_items_ls %>% add_unwtd_dim_tots(items_tb = .x,
            domain_pfx_1L_chr = prefix_chr[["domain_unwtd_pfx_1L_chr"]]) %>%
            add_wtd_dim_tots(domain_items_ls = domain_items_ls,
                domain_unwtd_pfx_1L_chr = prefix_chr[["domain_unwtd_pfx_1L_chr"]],
                domain_wtd_pfx_1L_chr = prefix_chr[["domain_wtd_pfx_1L_chr"]]) %>%
            add_labels_to_aqol6d_tb()
    }) %>% purrr::map(~.x %>% dplyr::select(!!rlang::sym(id_var_nm_1L_chr),
        dplyr::starts_with(item_pfx_1L_chr), dplyr::starts_with(prefix_chr[["domain_unwtd_pfx_1L_chr"]]),
        dplyr::starts_with(prefix_chr[["domain_wtd_pfx_1L_chr"]]),
        dplyr::everything()))
    return(aqol6d_adol_pop_tbs_ls)
}
#' Make Assessment of Quality of Life Six Dimension functions
#' @description make_aqol6d_fns_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make assessment of quality of life six dimension functions list. The function returns Assessment of Quality of Life Six Dimension disu (a list of functions).
#' @param domain_items_ls Domain items (a list)
#' @return Assessment of Quality of Life Six Dimension disu (a list of functions)
#' @rdname make_aqol6d_fns_ls
#' @export
#' @importFrom purrr map
#' @importFrom rlang sym
make_aqol6d_fns_ls <- function (domain_items_ls)
{
    aqol6d_disu_fn_ls <- paste0("calculate_aqol6d_dim_", 1:length(domain_items_ls),
        "_disv") %>% purrr::map(~rlang::sym(.x))
    return(aqol6d_disu_fn_ls)
}
#' Make Assessment of Quality of Life Six Dimension items
#' @description make_aqol6d_items_tb() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make assessment of quality of life six dimension items tibble. The function returns Assessment of Quality of Life Six Dimension items (a tibble).
#' @param aqol_tb Assessment of Quality of Life (a tibble)
#' @param old_pfx_1L_chr Old prefix (a character vector of length one)
#' @param new_pfx_1L_chr New prefix (a character vector of length one)
#' @return Assessment of Quality of Life Six Dimension items (a tibble)
#' @rdname make_aqol6d_items_tb
#' @export
#' @importFrom dplyr select starts_with rename_all
#' @importFrom stringr str_replace
make_aqol6d_items_tb <- function (aqol_tb, old_pfx_1L_chr, new_pfx_1L_chr)
{
    aqol6d_items_tb <- aqol_tb %>% dplyr::select(dplyr::starts_with(old_pfx_1L_chr)) %>%
        dplyr::rename_all(~{
            stringr::str_replace(., old_pfx_1L_chr, new_pfx_1L_chr)
        })
    return(aqol6d_items_tb)
}
#' Make brms mdl print list
#' @description make_brms_mdl_print_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make brms mdl print list. The function returns Brms mdl print (a list).
#' @param mdl_ls Mdl (a list)
#' @param label_stub_1L_chr Label stub (a character vector of length one)
#' @param caption_1L_chr Caption (a character vector of length one)
#' @param output_type_1L_chr Output type (a character vector of length one), Default: 'PDF'
#' @param digits_1L_dbl Digits (a double vector of length one), Default: 2
#' @param big_mark_1L_chr Big mark (a character vector of length one), Default: ' '
#' @return Brms mdl print (a list)
#' @rdname make_brms_mdl_print_ls
#' @export
#' @importFrom purrr map_dbl map_chr map2_chr discard map
#' @importFrom dplyr mutate all_of across case_when
#' @importFrom Hmisc latexTranslate
#' @importFrom stringr str_replace
#' @keywords internal
make_brms_mdl_print_ls <- function (mdl_ls, label_stub_1L_chr, caption_1L_chr, output_type_1L_chr = "PDF",
    digits_1L_dbl = 2, big_mark_1L_chr = " ")
{
    smry_mdl_ls <- summary(mdl_ls, digits = 4)
    mdl_smry_chr <- smry_mdl_ls %>% capture.output()
    idx_dbl <- c("Formula: ", "Samples: ", "Group-Level Effects: ",
        "Population-Level Effects: ", "Family Specific Parameters: ",
        "Samples were drawn using ") %>% purrr::map_dbl(~mdl_smry_chr %>%
        startsWith(.x) %>% which())
    data_tb <- make_brms_mdl_smry_tbl(smry_mdl_ls, grp_1L_chr = mdl_smry_chr[idx_dbl[3]],
        pop_1L_chr = mdl_smry_chr[idx_dbl[4]], fam_1L_chr = mdl_smry_chr[idx_dbl[5]])
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
            label = paste0("tab:", label_stub_1L_chr), merge_row_idx_int = as.integer(which(bold_lgl)),
            digits_dbl = c(ifelse(output_type_1L_chr == "PDF",
                0, NA_real_) %>% purrr::discard(is.na), names(data_tb) %>%
                purrr::map_dbl(~ifelse(.x %in% c("Bulk_ESS",
                  "Tail_ESS"), 0, digits_1L_dbl))), big_mark_1L_chr = big_mark_1L_chr,
            hline.after = c(-1, 0), sanitize_fn = force, footnotes_chr = NA_character_),
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
        brms_mdl_print_ls$part_5$addtorow <- list(pos = purrr::map(c(0,
            rep(nrow(data_tb), length(footnotes_chr) + 1)), ~.x),
            command = c(names(data_tb) %>% Hmisc::latexTranslate() %>%
                paste0(collapse = " & ") %>% paste0("\\\\\n"),
                c("\\toprule\n"), footnotes_chr %>% purrr::map_chr(~paste0("\\multicolumn{",
                  ncol(data_tb), "}{l}{", paste0("{\\footnotesize ",
                    .x, "}\n", collapse = ","), "}\\\\\n"))))
    }
    return(brms_mdl_print_ls)
}
#' Make brms mdl smry table
#' @description make_brms_mdl_smry_tbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make brms mdl smry table. The function returns Brms mdl smry (a tibble).
#' @param smry_mdl_ls Smry mdl (a list)
#' @param grp_1L_chr Grp (a character vector of length one)
#' @param pop_1L_chr Pop (a character vector of length one)
#' @param fam_1L_chr Fam (a character vector of length one)
#' @return Brms mdl smry (a tibble)
#' @rdname make_brms_mdl_smry_tbl
#' @export
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @keywords internal
make_brms_mdl_smry_tbl <- function (smry_mdl_ls, grp_1L_chr, pop_1L_chr, fam_1L_chr)
{
    brms_mdl_smry_tb <- purrr::map(1:length(smry_mdl_ls$random),
        ~make_mdl_smry_elmt_tbl(cat_chr = c(ifelse(.x == 1, grp_1L_chr,
            character(0)), paste0(names(smry_mdl_ls$ngrps)[.x],
            " (Number of levels: ", smry_mdl_ls$ngrps[.x][[1]],
            ")")), mat = smry_mdl_ls$random[.x][[1]])) %>% dplyr::bind_rows(make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$fixed,
        cat_chr = pop_1L_chr), make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$spec_pars,
        cat_chr = fam_1L_chr))
    return(brms_mdl_smry_tb)
}
#' Make complete props tibbles
#' @description make_complete_props_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make complete props tibbles list. The function returns Complete props tibbles (a list).
#' @param raw_props_tbs_ls Raw props tibbles (a list)
#' @param question_var_nm_1L_chr Question var name (a character vector of length one), Default: 'Question'
#' @return Complete props tibbles (a list)
#' @rdname make_complete_props_tbs_ls
#' @export
#' @importFrom purrr map map2_dbl
#' @importFrom dplyr mutate select mutate_if
#' @importFrom rlang sym
#' @keywords internal
make_complete_props_tbs_ls <- function (raw_props_tbs_ls, question_var_nm_1L_chr = "Question")
{
    complete_props_tbs_ls <- raw_props_tbs_ls %>% purrr::map(~{
        .x %>% dplyr::mutate(total_prop_dbl = rowSums(dplyr::select(.,
            -!!rlang::sym(question_var_nm_1L_chr)), na.rm = T) -
            100) %>% dplyr::mutate_if(is.numeric, ~purrr::map2_dbl(.,
            total_prop_dbl, ~ifelse(.x == 100, 1 - .y, .x))) %>%
            dplyr::select(-total_prop_dbl)
    })
    return(complete_props_tbs_ls)
}
#' Make correlated data
#' @description make_correlated_data_tb() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make correlated data tibble. The function returns Correlated data (a tibble).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param synth_data_idx_1L_dbl Synth data index (a double vector of length one), Default: 1
#' @return Correlated data (a tibble)
#' @rdname make_correlated_data_tb
#' @export
#' @importFrom simstudy genCorData
make_correlated_data_tb <- function (synth_data_spine_ls, synth_data_idx_1L_dbl = 1)
{
    correlated_data_tb <- simstudy::genCorData(synth_data_spine_ls$nbr_obs_dbl[synth_data_idx_1L_dbl],
        mu = synth_data_spine_ls$means_ls[[synth_data_idx_1L_dbl]],
        sigma = synth_data_spine_ls$sds_ls[[synth_data_idx_1L_dbl]],
        corMatrix = make_pdef_cor_mat_mat(synth_data_spine_ls$cor_mat_ls[[synth_data_idx_1L_dbl]]),
        cnames = synth_data_spine_ls$var_names_chr) %>% force_min_max_and_int_cnstrs(var_names_chr = synth_data_spine_ls$var_names_chr,
        min_max_ls = synth_data_spine_ls$min_max_ls, discrete_lgl = synth_data_spine_ls$discrete_lgl)
    return(correlated_data_tb)
}
#' Make corstars table
#' @description make_corstars_tbl_xx() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make corstars table output object of multiple potential types. The function is called for its side effects and does not return a value.
#' @param x PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: c("pearson", "spearman")
#' @param removeTriangle PARAM_DESCRIPTION, Default: c("upper", "lower")
#' @param result PARAM_DESCRIPTION, Default: c("none", "html", "latex")
#' @return NULL
#' @rdname make_corstars_tbl_xx
#' @export
#' @importFrom Hmisc rcorr
make_corstars_tbl_xx <- function (x, method = c("pearson", "spearman"), removeTriangle = c("upper",
    "lower"), result = c("none", "html", "latex"))
{
    x <- as.matrix(x)
    correlation_matrix <- Hmisc::rcorr(x, type = method[1])
    R <- correlation_matrix$r
    p <- correlation_matrix$P
    mystars <- ifelse(p < 1e-04, "****", ifelse(p < 0.001, "*** ",
        ifelse(p < 0.01, "**  ", ifelse(p < 0.05, "*   ", "    "))))
    R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[, -1]
    Rnew <- matrix(paste(R, mystars, sep = ""), ncol = ncol(x))
    diag(Rnew) <- paste(diag(R), " ", sep = "")
    rownames(Rnew) <- colnames(x)
    colnames(Rnew) <- paste(colnames(x), "", sep = "")
    if (removeTriangle[1] == "upper") {
        Rnew <- as.matrix(Rnew)
        Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
        Rnew <- as.data.frame(Rnew)
    }
    else if (removeTriangle[1] == "lower") {
        Rnew <- as.matrix(Rnew)
        Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
        Rnew <- as.data.frame(Rnew)
    }
    Rnew <- cbind(Rnew[1:length(Rnew) - 1])
    if (result[1] == "none")
        return(Rnew)
    else {
        if (result[1] == "html")
            print(xtable(Rnew), type = "html")
        else print(xtable(Rnew), type = "latex")
    }
}
#' Make dimension sclg constants
#' @description make_dim_sclg_cons_dbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make dimension sclg constants double vector. The function returns Dimension sclg constants (a double vector).
#' @param domains_chr Domains (a character vector)
#' @param dim_sclg_con_lup_tb Dimension sclg constant lookup table (a tibble)
#' @return Dimension sclg constants (a double vector)
#' @rdname make_dim_sclg_cons_dbl
#' @export
#' @importFrom purrr map_dbl
#' @importFrom ready4fun get_from_lup_obj
make_dim_sclg_cons_dbl <- function (domains_chr, dim_sclg_con_lup_tb)
{
    dim_sclg_cons_dbl <- purrr::map_dbl(domains_chr, ~ready4fun::get_from_lup_obj(dim_sclg_con_lup_tb,
        match_var_nm_1L_chr = "Dimension_chr", match_value_xx = .x,
        target_var_nm_1L_chr = "Constant_dbl", evaluate_lgl = F))
    return(dim_sclg_cons_dbl)
}
#' Make domain items
#' @description make_domain_items_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make domain items list. The function returns Domain items (a list).
#' @param domain_qs_lup_tb Domain questions lookup table (a tibble)
#' @param item_pfx_1L_chr Item prefix (a character vector of length one)
#' @return Domain items (a list)
#' @rdname make_domain_items_ls
#' @export
#' @importFrom purrr map
#' @importFrom dplyr filter pull
#' @importFrom stats setNames
make_domain_items_ls <- function (domain_qs_lup_tb, item_pfx_1L_chr)
{
    domains_chr <- domain_qs_lup_tb$Domain_chr %>% unique()
    q_nbrs_ls <- purrr::map(domains_chr, ~domain_qs_lup_tb %>%
        dplyr::filter(Domain_chr == .x) %>% dplyr::pull(Question_dbl))
    domain_items_ls <- purrr::map(q_nbrs_ls, ~paste0(item_pfx_1L_chr,
        .x)) %>% stats::setNames(domains_chr)
    return(domain_items_ls)
}
#' Make item wrst wghts
#' @description make_item_wrst_wghts_ls_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make item wrst wghts list list. The function returns Item wrst wghts (a list of lists).
#' @param domain_items_ls Domain items (a list)
#' @param itm_wrst_wghts_lup_tb Itm wrst wghts lookup table (a tibble)
#' @return Item wrst wghts (a list of lists)
#' @rdname make_item_wrst_wghts_ls_ls
#' @export
#' @importFrom purrr map map_dbl
#' @importFrom ready4fun get_from_lup_obj
make_item_wrst_wghts_ls_ls <- function (domain_items_ls, itm_wrst_wghts_lup_tb)
{
    item_wrst_wghts_ls_ls <- domain_items_ls %>% purrr::map(~{
        purrr::map_dbl(.x, ~{
            ready4fun::get_from_lup_obj(itm_wrst_wghts_lup_tb,
                match_var_nm_1L_chr = "Question_chr", match_value_xx = .x,
                target_var_nm_1L_chr = "Worst_Weight_dbl", evaluate_lgl = F)
        })
    })
    return(item_wrst_wghts_ls_ls)
}
#' Make knit parameters
#' @description make_knit_pars_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make knit parameters list. The function returns Knit parameters (a list).
#' @param mdl_smry_dir_1L_chr Mdl smry directory (a character vector of length one)
#' @param mdl_types_chr Mdl types (a character vector)
#' @param predictor_vars_nms_ls Predictor vars names (a list)
#' @param output_type_1L_chr Output type (a character vector of length one), Default: 'HTML'
#' @param mdl_types_lup Mdl types (a lookup table), Default: NULL
#' @param plt_types_lup Plt types (a lookup table), Default: NULL
#' @param plt_types_chr Plt types (a character vector), Default: c("coefs", "hetg", "dnst", "sctr_plt")
#' @param section_type_1L_chr Section type (a character vector of length one), Default: '#'
#' @return Knit parameters (a list)
#' @rdname make_knit_pars_ls
#' @export
#' @importFrom purrr pmap map
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
make_knit_pars_ls <- function (mdl_smry_dir_1L_chr, mdl_types_chr, predictor_vars_nms_ls,
    output_type_1L_chr = "HTML", mdl_types_lup = NULL, plt_types_lup = NULL,
    plt_types_chr = c("coefs", "hetg", "dnst", "sctr_plt"), section_type_1L_chr = "#")
{
    if (is.null(mdl_types_lup))
        data(mdl_types_lup, envir = environment())
    if (is.null(plt_types_lup))
        data(plt_types_lup, envir = environment())
    lab_idx_dbl <- 1:(length(mdl_types_chr) * length(predictor_vars_nms_ls))
    knit_pars_ls <- purrr::pmap(list(predictor_vars_nms_ls, split(lab_idx_dbl,
        ceiling(seq_along(lab_idx_dbl)/length(mdl_types_chr))),
        make_unique_ls_elmt_idx_int(predictor_vars_nms_ls), make_mdl_nms_ls(predictor_vars_nms_ls,
            mdl_types_chr = mdl_types_chr)), ~{
        mdl_nms_chr <- ..4
        path_to_mdl_stubs_chr <- paste0(mdl_smry_dir_1L_chr,
            "/", mdl_nms_chr)
        paths_to_mdls_chr <- paste0(path_to_mdl_stubs_chr, ".RDS")
        paths_to_mdl_plts_ls <- purrr::map(path_to_mdl_stubs_chr,
            ~paste0(..1, paste0("_", plt_types_chr, ".png")))
        mdl_ttls_chr <- paste0(..1[1], ifelse(is.na(..1[2]),
            "", paste(" with ", ..1[2])), " ", ready4fun::get_from_lup_obj(mdl_types_lup,
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_types_chr,
            target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F))
        section_ttls_chr <- paste0(section_type_1L_chr, " ",
            mdl_ttls_chr)
        list(plt_nms_ls = purrr::map(mdl_ttls_chr, ~{
            paste0(.x, " ", ready4fun::get_from_lup_obj(plt_types_lup,
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = plt_types_chr,
                target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F))
        }), paths_to_mdls_chr = paths_to_mdls_chr, tbl_captions_chr = mdl_ttls_chr,
            label_stubs_chr = paste0("lab", ..2), output_type_1L_chr = output_type_1L_chr,
            section_ttls_chr = section_ttls_chr, ls_elmt_idx_1L_int = ..3)
    })
    return(knit_pars_ls)
}
#' Make mdl names
#' @description make_mdl_nms_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make mdl names list. The function returns Mdl names (a list).
#' @param predictor_vars_nms_ls Predictor vars names (a list)
#' @param mdl_types_chr Mdl types (a character vector)
#' @return Mdl names (a list)
#' @rdname make_mdl_nms_ls
#' @export
#' @importFrom purrr map2
#' @keywords internal
make_mdl_nms_ls <- function (predictor_vars_nms_ls, mdl_types_chr)
{
    mdl_nms_ls <- purrr::map2(predictor_vars_nms_ls, make_unique_ls_elmt_idx_int(predictor_vars_nms_ls),
        ~paste0(.x[1], "_", ifelse(is.na(.x[2]), "", paste0(.x[2],
            "_")), .y, "_", mdl_types_chr))
    return(mdl_nms_ls)
}
#' Make mdl smry elmt table
#' @description make_mdl_smry_elmt_tbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make mdl smry elmt table. The function returns Mdl elmt sum (a tibble).
#' @param mat Matrix (a matrix)
#' @param cat_chr Cat (a character vector)
#' @return Mdl elmt sum (a tibble)
#' @rdname make_mdl_smry_elmt_tbl
#' @export
#' @importFrom tibble as_tibble add_case
#' @importFrom dplyr mutate select everything filter bind_rows
#' @keywords internal
make_mdl_smry_elmt_tbl <- function (mat, cat_chr)
{
    tb <- mat %>% tibble::as_tibble() %>% dplyr::mutate(Parameter = rownames(mat)) %>%
        dplyr::select(Parameter, dplyr::everything())
    mdl_elmt_sum_tb <- tb %>% dplyr::filter(F) %>% tibble::add_case(Parameter = cat_chr) %>%
        dplyr::bind_rows(tb)
    return(mdl_elmt_sum_tb)
}
#' Make pdef correlation matrix
#' @description make_pdef_cor_mat_mat() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make pdef correlation matrix matrix. The function returns Pdef correlation (a matrix).
#' @param lower_diag_mat Lower diag (a matrix)
#' @return Pdef correlation (a matrix)
#' @rdname make_pdef_cor_mat_mat
#' @export
#' @importFrom Matrix forceSymmetric
#' @importFrom matrixcalc is.positive.definite
#' @importFrom psych cor.smooth
make_pdef_cor_mat_mat <- function (lower_diag_mat)
{
    pdef_cor_mat <- lower_diag_mat %>% Matrix::forceSymmetric(uplo = "L") %>%
        as.matrix()
    if (!matrixcalc::is.positive.definite(pdef_cor_mat)) {
        pdef_cor_mat <- psych::cor.smooth(pdef_cor_mat)
    }
    return(pdef_cor_mat)
}
#' Make smry of brm mdl
#' @description make_smry_of_brm_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make smry of brm mdl. The function returns Smry of brm mdl (a tibble).
#' @param mdl_ls Mdl (a list)
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param predictor_vars_nms_chr Predictor vars names (a character vector)
#' @param fn Function (a function), Default: calculate_rmse
#' @param mdl_nm_1L_chr Mdl name (a character vector of length one), Default: 'NA'
#' @param seed_1L_dbl Seed (a double vector of length one), Default: 23456
#' @return Smry of brm mdl (a tibble)
#' @rdname make_smry_of_brm_mdl
#' @export
#' @importFrom brms bayes_R2
#' @importFrom psych describe
#' @importFrom dplyr pull mutate rename select
#' @importFrom rlang sym
#' @importFrom purrr map flatten_chr
#' @keywords internal
make_smry_of_brm_mdl <- function (mdl_ls, data_tb, dep_var_nm_1L_chr = "aqol6d_total_w",
    predictor_vars_nms_chr, fn = calculate_rmse, mdl_nm_1L_chr = NA_character_,
    seed_1L_dbl = 23456)
{
    if (is.na(mdl_nm_1L_chr))
        mdl_nm_1L_chr <- predictor_vars_nms_chr[1]
    set.seed(seed_1L_dbl)
    predictions <- predict(mdl_ls, summary = F)
    coef <- summary(mdl_ls, digits = 4)$fixed
    coef <- coef[1:nrow(coef), 1:4]
    R2 <- brms::bayes_R2(mdl_ls)
    RMSE <- psych::describe(apply(predictions, 1, fn, y_dbl = data_tb %>%
        dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))), quant = c(0.25,
        0.75), skew = F, ranges = F)
    RMSE <- cbind(RMSE$mean, RMSE$sd, RMSE$Q0.25, RMSE$Q0.75)
    Sigma <- summary(mdl_ls, digits = 4)$spec_par[1:4]
    smry_of_brm_mdl_tb <- data.frame(round(rbind(coef, R2, RMSE,
        Sigma), 3)) %>% dplyr::mutate(Parameter = c("Intercept",
        purrr::map(predictor_vars_nms_chr, ~paste0(.x, c(" baseline",
            " change"))) %>% purrr::flatten_chr(), "R2", "RMSE",
        "Sigma"), Model = mdl_nm_1L_chr) %>% dplyr::mutate(`95% CI` = paste(l.95..CI,
        ",", u.95..CI)) %>% dplyr::rename(SE = Est.Error) %>%
        dplyr::select(Model, Parameter, Estimate, SE, `95% CI`)
    return(smry_of_brm_mdl_tb)
}
#' Make smry of ts mdl
#' @description make_smry_of_ts_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make smry of ts mdl. The function returns Smry of ts mdl (a list).
#' @param data_tb Data (a tibble)
#' @param fn Function (a function)
#' @param predictor_vars_nms_chr Predictor vars names (a character vector)
#' @param mdl_nm_1L_chr Mdl name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one), Default: 'NA'
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @param iters_1L_int Iters (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Smry of ts mdl (a list)
#' @rdname make_smry_of_ts_mdl
#' @export
#' @importFrom rlang exec
#' @keywords internal
make_smry_of_ts_mdl <- function (data_tb, fn, predictor_vars_nms_chr, mdl_nm_1L_chr,
    path_to_write_to_1L_chr = NA_character_, dep_var_nm_1L_chr = "aqol6d_total_w",
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round",
    round_bl_val_1L_chr = "Baseline", iters_1L_int = 4000L, seed_1L_int = 1000L)
{
    tfd_data_tb <- transform_tb_to_mdl_inp(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        predictor_vars_nms_chr = predictor_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr,
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr)
    tfd_dep_var_nm_1L_chr <- ifelse(identical(fn, fit_clg_log_tfmn),
        transform_dep_var_nm_for_cll(dep_var_nm_1L_chr), dep_var_nm_1L_chr)
    args_ls <- list(data_tb = tfd_data_tb, dep_var_nm_1L_chr = tfd_dep_var_nm_1L_chr,
        predictor_vars_nms_chr = predictor_vars_nms_chr, iters_1L_int = iters_1L_int,
        seed_1L_int = seed_1L_int)
    mdl_ls <- rlang::exec(fn, !!!args_ls)
    smry_of_ts_mdl_ls <- list(smry_of_ts_mdl_tb = make_smry_of_brm_mdl(mdl_ls,
        data_tb = tfd_data_tb, dep_var_nm_1L_chr = tfd_dep_var_nm_1L_chr,
        predictor_vars_nms_chr = predictor_vars_nms_chr, fn = ifelse(identical(fn,
            fit_gsn_log_lnk), calculate_rmse, calculate_rmse_tfmn),
        mdl_nm_1L_chr = mdl_nm_1L_chr))
    if (!is.na(path_to_write_to_1L_chr)) {
        smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr <- paste0(path_to_write_to_1L_chr,
            "/", mdl_nm_1L_chr, ".RDS")
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
#' Make synth series tibbles
#' @description make_synth_series_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make synth series tibbles list. The function returns Synth series tibbles (a list).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param series_names_chr Series names (a character vector)
#' @return Synth series tibbles (a list)
#' @rdname make_synth_series_tbs_ls
#' @export
#' @importFrom purrr map
#' @importFrom stats setNames
make_synth_series_tbs_ls <- function (synth_data_spine_ls, series_names_chr)
{
    synth_series_tbs_ls <- 1:length(series_names_chr) %>% purrr::map(~make_correlated_data_tb(synth_data_spine_ls = synth_data_spine_ls,
        synth_data_idx_1L_dbl = .x) %>% replace_with_missing_vals(synth_data_spine_ls = synth_data_spine_ls,
        idx_int = .x)) %>% stats::setNames(series_names_chr)
    return(synth_series_tbs_ls)
}
#' Make unique list elmt index
#' @description make_unique_ls_elmt_idx_int() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make unique list elmt index integer vector. The function returns Unique list elmt index (an integer vector).
#' @param data_ls Data (a list)
#' @return Unique list elmt index (an integer vector)
#' @rdname make_unique_ls_elmt_idx_int
#' @export
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate group_by row_number ungroup
#' @importFrom purrr map flatten_int
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
make_unique_ls_elmt_idx_int <- function (data_ls)
{
    combos_tb <- tibble::as_tibble(data_ls, .name_repair = ~paste0("r_",
        1:length(data_ls))) %>% t() %>% as.data.frame()
    combos_tb <- combos_tb %>% tibble::as_tibble()
    combos_tb <- combos_tb %>% dplyr::mutate(combo_chr = ifelse(ncol(combos_tb) ==
        1, V1, paste0(V1, V2)))
    combos_tb <- combos_tb %>% dplyr::group_by(combo_chr) %>%
        dplyr::mutate(combo_id = dplyr::row_number())
    unique_ls_elmt_idx_int <- purrr::map(data_ls %>% unique(),
        ~ready4fun::get_from_lup_obj(combos_tb %>% dplyr::ungroup(),
            match_var_nm_1L_chr = "combo_chr", match_value_xx = paste0(.x[1],
                ifelse(is.na(.x[2]), "", .x[2])), target_var_nm_1L_chr = "combo_id",
            evaluate_lgl = F)) %>% purrr::flatten_int()
    return(unique_ls_elmt_idx_int)
}
#' Make vec with sum of
#' @description make_vec_with_sum_of_int() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make vec with sum of integer vector. The function returns Vec (an integer vector).
#' @param target_int Target (an integer vector)
#' @param start_int Start (an integer vector)
#' @param end_int End (an integer vector)
#' @param length_int Length (an integer vector)
#' @return Vec (an integer vector)
#' @rdname make_vec_with_sum_of_int
#' @export
#' @importFrom Surrogate RandVec
#' @importFrom purrr pluck
make_vec_with_sum_of_int <- function (target_int, start_int, end_int, length_int)
{
    vec_int <- Surrogate::RandVec(a = start_int, b = end_int,
        s = target_int, n = length_int, m = 1) %>% purrr::pluck("RandVecOutput") %>%
        as.vector() %>% round() %>% as.integer() %>% force_vec_to_sum_to_int(target_1L_int = target_int,
        min_max_int = c(start_int, end_int))
    return(vec_int)
}
