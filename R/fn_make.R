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
#' @keywords internal
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
#' @keywords internal
make_aqol6d_items_tb <- function (aqol_tb, old_pfx_1L_chr, new_pfx_1L_chr) 
{
    aqol6d_items_tb <- aqol_tb %>% dplyr::select(dplyr::starts_with(old_pfx_1L_chr)) %>% 
        dplyr::rename_all(~{
            stringr::str_replace(., old_pfx_1L_chr, new_pfx_1L_chr)
        })
    return(aqol6d_items_tb)
}
#' Make brms model print list
#' @description make_brms_mdl_print_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make brms model print list. The function returns Brms model print (a list).
#' @param mdl_ls Model (a list)
#' @param label_stub_1L_chr Label stub (a character vector of length one)
#' @param caption_1L_chr Caption (a character vector of length one)
#' @param output_type_1L_chr Output type (a character vector of length one), Default: 'PDF'
#' @param digits_1L_dbl Digits (a double vector of length one), Default: 2
#' @param big_mark_1L_chr Big mark (a character vector of length one), Default: ' '
#' @return Brms model print (a list)
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
#' Make brms model smry table
#' @description make_brms_mdl_smry_tbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make brms model smry table. The function returns Brms model smry (a tibble).
#' @param smry_mdl_ls Smry model (a list)
#' @param grp_1L_chr Grp (a character vector of length one)
#' @param pop_1L_chr Pop (a character vector of length one)
#' @param fam_1L_chr Fam (a character vector of length one)
#' @return Brms model smry (a tibble)
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
make_domain_items_ls <- function (domain_qs_lup_tb, item_pfx_1L_chr) 
{
    domains_chr <- domain_qs_lup_tb$Domain_chr %>% unique()
    q_nbrs_ls <- purrr::map(domains_chr, ~domain_qs_lup_tb %>% 
        dplyr::filter(Domain_chr == .x) %>% dplyr::pull(Question_dbl))
    domain_items_ls <- purrr::map(q_nbrs_ls, ~paste0(item_pfx_1L_chr, 
        .x)) %>% stats::setNames(domains_chr)
    return(domain_items_ls)
}
#' Make fake time series data
#' @description make_fake_ts_data() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make fake time series data. The function returns Fk data (a tibble).
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
#' Make folds
#' @description make_folds_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make folds list. The function returns Folds (a list).
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param n_folds_1L_int N folds (an integer vector of length one), Default: 10
#' @return Folds (a list)
#' @rdname make_folds_ls
#' @export 
#' @importFrom caret createFolds
#' @importFrom dplyr pull
#' @importFrom rlang sym
#' @keywords internal
make_folds_ls <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", n_folds_1L_int = 10L) 
{
    folds_ls <- caret::createFolds(data_tb %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr)), 
        k = n_folds_1L_int, list = TRUE, returnTrain = FALSE)
    return(folds_ls)
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
#' @keywords internal
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
#' @param rltv_path_to_data_dir_1L_chr Rltv path to data directory (a character vector of length one)
#' @param mdl_types_chr Model types (a character vector)
#' @param predr_vars_nms_ls Predr vars names (a list)
#' @param output_type_1L_chr Output type (a character vector of length one), Default: 'HTML'
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @param plt_types_lup Plt types (a lookup table), Default: NULL
#' @param plt_types_chr Plt types (a character vector), Default: 'NA'
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
#' @keywords internal
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
#' @keywords internal
make_mdl_nms_ls <- function (predr_vars_nms_ls, mdl_types_chr) 
{
    mdl_nms_ls <- purrr::map2(predr_vars_nms_ls, make_unique_ls_elmt_idx_int(predr_vars_nms_ls), 
        ~paste0(.x[1], "_", ifelse(is.na(.x[2]), "", paste0(.x[2], 
            "_")), .y, "_", mdl_types_chr))
    return(mdl_nms_ls)
}
#' Make model smry elmt table
#' @description make_mdl_smry_elmt_tbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make model smry elmt table. The function returns Model elmt sum (a tibble).
#' @param mat Matrix (a matrix)
#' @param cat_chr Cat (a character vector)
#' @return Model elmt sum (a tibble)
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
#' @keywords internal
make_pdef_cor_mat_mat <- function (lower_diag_mat) 
{
    pdef_cor_mat <- lower_diag_mat %>% Matrix::forceSymmetric(uplo = "L") %>% 
        as.matrix()
    if (!matrixcalc::is.positive.definite(pdef_cor_mat)) {
        pdef_cor_mat <- psych::cor.smooth(pdef_cor_mat)
    }
    return(pdef_cor_mat)
}
#' Make prediction dataset with one predr
#' @description make_predn_ds_with_one_predr() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make prediction dataset with one predr. The function returns Prediction dataset (a tibble).
#' @param model_mdl PARAM_DESCRIPTION
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param predr_var_nm_1L_chr Predr var name (a character vector of length one)
#' @param predr_vals_dbl Predr values (a double vector)
#' @param pred_type_1L_chr Pred type (a character vector of length one), Default: NULL
#' @return Prediction dataset (a tibble)
#' @rdname make_predn_ds_with_one_predr
#' @export 
#' @importFrom tibble tibble
#' @importFrom rlang sym
#' @importFrom dplyr mutate
#' @importFrom stats predict
#' @keywords internal
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
#' Make predr values
#' @description make_predr_vals() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make predr values. The function returns Predr values (a double vector).
#' @param predr_var_nm_1L_chr Predr var name (a character vector of length one)
#' @param candidate_predrs_lup Candidate predrs (a lookup table), Default: NULL
#' @return Predr values (a double vector)
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
#' Make predr vars names
#' @description make_predr_vars_nms_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make predr vars names list. The function returns Predr vars names (a list).
#' @param main_predrs_chr Main predrs (a character vector)
#' @param covars_ls Covars (a list)
#' @return Predr vars names (a list)
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
#' Make prefd models vec
#' @description make_prefd_mdls_vec() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make prefd models vec. The function returns Prefd models (a character vector).
#' @param smry_of_sngl_predr_mdls_tb Smry of sngl predr models (a tibble)
#' @param choose_from_pfx_chr Choose from prefix (a character vector), Default: c("GLM", "OLS")
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @return Prefd models (a character vector)
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
        model_mdl <- make_mdl(data_tb[-.x, ], dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
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
#' Make smry of time series model
#' @description make_smry_of_ts_mdl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make smry of time series model. The function returns Smry of time series model (a list).
#' @param data_tb Data (a tibble)
#' @param fn Function (a function)
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one), Default: 'NA'
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @param predictors_lup Predictors (a lookup table)
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iters (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Smry of time series model (a list)
#' @rdname make_smry_of_ts_mdl
#' @export 
#' @importFrom purrr map_dbl
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom rlang exec
#' @keywords internal
make_smry_of_ts_mdl <- function (data_tb, fn, predr_vars_nms_chr, mdl_nm_1L_chr, path_to_write_to_1L_chr = NA_character_, 
    dep_var_nm_1L_chr = "utl_total_w", id_var_nm_1L_chr = "fkClientID", 
    round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", 
    predictors_lup, backend_1L_chr = getOption("brms.backend", 
        "rstan"), iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    scaling_fctr_dbl <- predr_vars_nms_chr %>% purrr::map_dbl(~ifelse(.x %in% 
        predictors_lup$short_name_chr, ready4fun::get_from_lup_obj(predictors_lup, 
        target_var_nm_1L_chr = mdl_scaling_dbl, match_value_xx = .x, 
        match_var_nm_1L_chr = "short_name_chr", evaluate_lgl = F), 
        1))
    tfd_data_tb <- transform_tb_to_mdl_inp(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr, 
        scaling_fctr_dbl = scaling_fctr_dbl)
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
#' Make synth series tibbles
#' @description make_synth_series_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make synth series tibbles list. The function returns Synth series tibbles (a list).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param series_names_chr Series names (a character vector)
#' @return Synth series tibbles (a list)
#' @rdname make_synth_series_tbs_ls
#' @export 
#' @importFrom purrr map
#' @importFrom stats setNames
#' @keywords internal
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
#' @keywords internal
make_vec_with_sum_of_int <- function (target_int, start_int, end_int, length_int) 
{
    vec_int <- Surrogate::RandVec(a = start_int, b = end_int, 
        s = target_int, n = length_int, m = 1) %>% purrr::pluck("RandVecOutput") %>% 
        as.vector() %>% round() %>% as.integer() %>% force_vec_to_sum_to_int(target_1L_int = target_int)
    return(vec_int)
}
