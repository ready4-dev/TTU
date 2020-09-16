#' Make aqol items props tibbles
#' @description make_aqol_items_props_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make aqol items props a tibbles. The function returns aqol items props tibbles (a list).

#' @return Aqol items props tibbles (a list)
#' @rdname make_aqol_items_props_tbs_ls
#' @export 
#' @importFrom tibble tribble
#' @importFrom dplyr mutate select mutate_if
#' @importFrom purrr map2_dbl
#' @keywords internal
make_aqol_items_props_tbs_ls <- function () 
{
    bl_answer_props_tb <- tibble::tribble(~Question, ~Answer_1, 
        ~Answer_2, ~Answer_3, ~Answer_4, ~Answer_5, ~Answer_6, 
        "Q1", 0.35, 0.38, 0.16, 0.03, 100, NA_real_, "Q2", 0.28, 
        0.38, 0.18, 0.08, 0.04, 100, "Q3", 0.78, 0.18, 0.03, 
        0.01, 0, 100, "Q4", 0.64, 0.23, 0.09, 0, 100, NA_real_, 
        "Q5", 0.3, 0.48, 0.12, 0.05, 100, NA_real_, "Q6", 0.33, 
        0.48, 0.15, 100, NA_real_, NA_real_, "Q7", 0.44, 0.27, 
        0.11, 100, NA_real_, NA_real_, "Q8", 0.18, 0.29, 0.23, 
        0.21, 100, NA_real_, "Q9", 0.07, 0.27, 0.19, 0.37, 100, 
        NA_real_, "Q10", 0.04, 0.15, 0.4, 0.25, 100, NA_real_, 
        "Q11", 0.03, 0.13, 0.52, 0.25, 100, NA_real_, "Q12", 
        0.06, 0.21, 0.25, 0.34, 100, NA_real_, "Q13", 0.05, 0.25, 
        0.31, 0.28, 100, NA_real_, "Q14", 0.05, 0.3, 0.34, 0.25, 
        100, NA_real_, "Q15", 0.57, 0.25, 0.12, 100, NA_real_, 
        NA_real_, "Q16", 0.48, 0.42, 0.06, 100, NA_real_, NA_real_, 
        "Q17", 0.44, 0.3, 0.16, 0.07, 100, NA_real_, "Q18", 0.33, 
        0.38, 0.25, 0.04, 0, 100, "Q19", 0.33, 0.49, 0.16, 0.02, 
        0, 100, "Q20", 0.67, 0.21, 0.02, 100, NA_real_, NA_real_) %>% 
        dplyr::mutate(total_prop_dbl = rowSums(dplyr::select(., 
            -Question), na.rm = T) - 100) %>% dplyr::mutate_if(is.numeric, 
        ~purrr::map2_dbl(., total_prop_dbl, ~ifelse(.x == 100, 
            1 - .y, .x))) %>% dplyr::select(-total_prop_dbl)
    fup_answer_props_tb <- tibble::tribble(~Question, ~Answer_1, 
        ~Answer_2, ~Answer_3, ~Answer_4, ~Answer_5, ~Answer_6, 
        "Q1", 0.51, 0.33, 0.12, 0.02, 100, NA_real_, "Q2", 0.36, 
        0.38, 0.16, 0.06, 0.02, 100, "Q3", 0.81, 0.15, 0.04, 
        0, 0, 100, "Q4", 0.73, 0.18, 0.09, 0, 100, NA_real_, 
        "Q5", 0.36, 0.42, 0.12, 0.05, 100, NA_real_, "Q6", 0.48, 
        0.4, 0.11, 100, NA_real_, NA_real_, "Q7", 0.57, 0.25, 
        0.09, 100, NA_real_, NA_real_, "Q8", 0.31, 0.33, 0.17, 
        0.12, 100, NA_real_, "Q9", 0.13, 0.35, 0.19, 0.23, 100, 
        NA_real_, "Q10", 0.1, 0.21, 0.43, 0.16, 100, NA_real_, 
        "Q11", 0.06, 0.25, 0.48, 0.18, 100, NA_real_, "Q12", 
        0.08, 0.27, 0.26, 0.25, 100, NA_real_, "Q13", 0.07, 0.37, 
        0.31, 0.19, 100, NA_real_, "Q14", 0.08, 0.37, 0.34, 0.15, 
        100, NA_real_, "Q15", 0.62, 0.23, 0.09, 100, NA_real_, 
        NA_real_, "Q16", 0.52, 0.4, 0.06, 100, NA_real_, NA_real_, 
        "Q17", 0.51, 0.28, 0.15, 0.06, 100, NA_real_, "Q18", 
        0.37, 0.35, 0.25, 0.03, 0, 100, "Q19", 0.43, 0.4, 0.16, 
        0.01, 0, 100, "Q20", 0.77, 0.21, 0.02, 100, NA_real_, 
        NA_real_) %>% dplyr::mutate(total_prop_dbl = rowSums(dplyr::select(., 
        -Question), na.rm = T) - 100) %>% dplyr::mutate_if(is.numeric, 
        ~purrr::map2_dbl(., total_prop_dbl, ~ifelse(.x == 100, 
            1 - .y, .x))) %>% dplyr::select(-total_prop_dbl)
    aqol_items_props_tbs_ls <- list(bl_answer_props_tb, fup_answer_props_tb)
    return(aqol_items_props_tbs_ls)
}
#' Make aqol6d functions
#' @description make_aqol6d_fns_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make aqol6d a functions. The function returns aqol6d disu (a list of functions).
#' @param domain_items_ls Domain items (a list)
#' @return Aqol6d disu (a list of functions)
#' @rdname make_aqol6d_fns_ls
#' @export 
#' @importFrom purrr map
#' @importFrom rlang sym
#' @keywords internal
make_aqol6d_fns_ls <- function (domain_items_ls) 
{
    aqol6d_disu_fn_ls <- paste0("calculate_aqol6d_d", 1:length(domain_items_ls), 
        "_disu_dbl") %>% purrr::map(~rlang::sym(.x))
    return(aqol6d_disu_fn_ls)
}
#' Make aqol6d items
#' @description make_aqol6d_items_tb() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make aqol6d items. The function returns aqol6d items (a tibble).
#' @param aqol_tb Aqol (a tibble)
#' @param old_pfx_1L_chr Old prefix 1L (a character vector of length 1)
#' @param new_pfx_1L_chr New prefix 1L (a character vector of length 1)
#' @return Aqol6d items (a tibble)
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
#' Make correlated data
#' @description make_correlated_data_tb() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make correlated data. The function returns correlated data (a tibble).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param synth_data_idx_1L_dbl Synth data idx 1L (a double vector of length 1), Default: 1
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
        corMatrix = make_pdef_corr_mat_mat(synth_data_spine_ls$corr_mat_ls[[synth_data_idx_1L_dbl]]), 
        cnames = synth_data_spine_ls$var_names_chr) %>% force_min_max_and_int_cnstrs_tb(var_names_chr = synth_data_spine_ls$var_names_chr, 
        min_max_ls = synth_data_spine_ls$min_max_ls, discrete_lgl = synth_data_spine_ls$discrete_lgl)
    return(correlated_data_tb)
}
#' Make corstars table
#' @description make_corstars_tbl_xx() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make corstars a table. The function is called for its side effects and does not return a value.
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
#' Make dim sclg cons
#' @description make_dim_sclg_cons_dbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make dim sclg cons. The function returns dim sclg cons (a double vector of length 1).
#' @param domains_chr Domains (a character vector of length 1)
#' @param dim_sclg_constant_lup_tb Dim sclg constant lookup table (a tibble), Default: dim_sclg_constant_lup_tb
#' @return Dim sclg cons (a double vector of length 1)
#' @rdname make_dim_sclg_cons_dbl
#' @export 
#' @importFrom purrr map_dbl
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
make_dim_sclg_cons_dbl <- function (domains_chr, dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb) 
{
    dim_sclg_cons_dbl <- purrr::map_dbl(domains_chr, ~ready4fun::get_from_lup_obj(dim_sclg_constant_lup_tb, 
        match_var_nm_1L_chr = "Dimension_chr", match_value_xx = .x, 
        target_var_nm_1L_chr = "Constant_dbl", evaluate_lgl = F))
    return(dim_sclg_cons_dbl)
}
#' Make domain items
#' @description make_domain_items_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make domain items. The function returns domain items (a list).
#' @param domains_chr Domains (a character vector of length 1)
#' @param q_nbrs_ls Q nbrs (a list)
#' @param item_pfx_1L_chr Item prefix 1L (a character vector of length 1)
#' @return Domain items (a list)
#' @rdname make_domain_items_ls
#' @export 
#' @importFrom purrr map
#' @importFrom stats setNames
#' @keywords internal
make_domain_items_ls <- function (domains_chr, q_nbrs_ls, item_pfx_1L_chr) 
{
    domain_items_ls <- purrr::map(q_nbrs_ls, ~paste0(item_pfx_1L_chr, 
        .x)) %>% stats::setNames(domains_chr)
    return(domain_items_ls)
}
#' Make item wrst wghts
#' @description make_item_wrst_wghts_ls_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make item wrst wghts. The function returns item wrst wghts (a list of lists).
#' @param domain_items_ls Domain items (a list)
#' @param itm_wrst_wghts_lup_tb Itm wrst wghts lookup table (a tibble), Default: itm_wrst_wghts_lup_tb
#' @return Item wrst wghts (a list of lists)
#' @rdname make_item_wrst_wghts_ls_ls
#' @export 
#' @importFrom purrr map map_dbl
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
make_item_wrst_wghts_ls_ls <- function (domain_items_ls, itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb) 
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
#' Make pdef corr matrix
#' @description make_pdef_corr_mat_mat() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make pdef corr a matrix. The function returns pdef corr (a matrix).
#' @param lower_diag_mat Lower diag (a matrix)
#' @return Pdef corr (a matrix)
#' @rdname make_pdef_corr_mat_mat
#' @export 
#' @importFrom Matrix forceSymmetric
#' @importFrom matrixcalc is.positive.definite
#' @importFrom psych cor.smooth
#' @keywords internal
make_pdef_corr_mat_mat <- function (lower_diag_mat) 
{
    pdef_corr_mat <- lower_diag_mat %>% Matrix::forceSymmetric(uplo = "L") %>% 
        as.matrix()
    if (!matrixcalc::is.positive.definite(pdef_corr_mat)) {
        pdef_corr_mat <- psych::cor.smooth(pdef_corr_mat)
    }
    return(pdef_corr_mat)
}
#' Make synth series tibbles
#' @description make_synth_series_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make synth series a tibbles. The function returns synth series tibbles (a list).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param series_names_chr Series names (a character vector of length 1)
#' @return Synth series tibbles (a list)
#' @rdname make_synth_series_tbs_ls
#' @export 
#' @importFrom purrr map
#' @importFrom stats setNames
#' @keywords internal
make_synth_series_tbs_ls <- function (synth_data_spine_ls, series_names_chr) 
{
    synth_series_tbs_ls <- 1:length(series_names_chr) %>% purrr::map(~make_correlated_data_tb(synth_data_spine_ls = synth_data_spine_ls, 
        synth_data_idx_1L_dbl = .x) %>% replace_var_vals_with_missing_tbl(synth_data_spine_ls = synth_data_spine_ls, 
        idx_int = .x)) %>% stats::setNames(series_names_chr)
    return(synth_series_tbs_ls)
}
#' Make vec with sum of
#' @description make_vec_with_sum_of_int() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make vec with sum of. The function returns vec (an integer vector of length 1).
#' @param target_int Target (an integer vector of length 1)
#' @param start_int Start (an integer vector of length 1)
#' @param end_int End (an integer vector of length 1)
#' @param length_int Length (an integer vector of length 1)
#' @return Vec (an integer vector of length 1)
#' @rdname make_vec_with_sum_of_int
#' @export 
#' @importFrom Surrogate RandVec
#' @importFrom purrr pluck
#' @keywords internal
make_vec_with_sum_of_int <- function (target_int, start_int, end_int, length_int) 
{
    vec_int <- Surrogate::RandVec(a = start_int, b = end_int, 
        s = target_int, n = length_int, m = 1) %>% purrr::pluck("RandVecOutput") %>% 
        as.vector() %>% round() %>% as.integer() %>% force_vec_to_sum_to_int(target_1L_int = target_int, 
        min_max_int = c(start_int, end_int))
    return(vec_int)
}
