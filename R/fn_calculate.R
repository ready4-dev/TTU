#' Calculate adolescent Assessment of Quality of Life Six Dimension health utility
#' @description calculate_adol_aqol6d() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate adolescent assessment of quality of life six dimension health utility. The function returns Adolescent Assessment of Quality of Life Six Dimension health utility (a double vector).
#' @param unscored_aqol_tb Unscored Assessment of Quality of Life health utility (a tibble)
#' @param prefix_1L_chr Prefix (a character vector of length one), Default: 'aqol'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one)
#' @return Adolescent Assessment of Quality of Life Six Dimension health utility (a double vector)
#' @rdname calculate_adol_aqol6d
#' @export 
#' @importFrom dplyr select starts_with rename_all
#' @importFrom stringr str_replace
calculate_adol_aqol6d <- function (unscored_aqol_tb, prefix_1L_chr = "aqol", id_var_nm_1L_chr) 
{
    unscored_aqol_tb <- unscored_aqol_tb %>% dplyr::select(id_var_nm_1L_chr, 
        dplyr::starts_with(prefix_1L_chr))
    names(unscored_aqol_tb) <- c("ID", paste0("Q", 1:20))
    unscored_aqol_tb <- impute_adol_unscrd_aqol_ds(unscored_aqol_tb)
    disvals_tb <- unscored_aqol_tb %>% add_itm_disu_to_aqol6d_itms_tb_tb(disvalues_lup_tb = make_adol_aqol6d_disv_lup(), 
        pfx_1L_chr = "Q") %>% dplyr::select(ID, dplyr::starts_with("dv_")) %>% 
        dplyr::rename_all(~stringr::str_replace(.x, "dv_", "dv"))
    scored_aqol_tb <- add_aqol_dim_scrg_eqs(disvals_tb)
    adol_aqol6d_dbl <- scored_aqol_tb$uaqol
    return(adol_aqol6d_dbl)
}
#' Calculate Assessment of Quality of Life Six Dimension health utility d1 disu
#' @description calculate_aqol6d_d1_disu_dbl() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate assessment of quality of life six dimension health utility d1 disu double vector. The function returns DvD1 (a double vector).
#' @param dvQs_tb DvQs (a tibble)
#' @param kD_1L_dbl KD (a double vector of length one)
#' @param w_dbl W (a double vector)
#' @return DvD1 (a double vector)
#' @rdname calculate_aqol6d_d1_disu_dbl
#' @export 
#' @importFrom purrr pmap_dbl
calculate_aqol6d_d1_disu_dbl <- function (dvQs_tb, kD_1L_dbl, w_dbl) 
{
    dvD1_dbl <- purrr::pmap_dbl(dvQs_tb, ~{
        (1/kD_1L_dbl) * ((1 + (kD_1L_dbl * w_dbl[1] * ..1)) * 
            (1 + (kD_1L_dbl * w_dbl[2] * ..2)) * (1 + (kD_1L_dbl * 
            w_dbl[3] * ..3)) * (1 + (kD_1L_dbl * w_dbl[4] * ..4)) - 
            1)
    })
    return(dvD1_dbl)
}
#' Calculate Assessment of Quality of Life Six Dimension health utility d2 disu
#' @description calculate_aqol6d_d2_disu_dbl() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate assessment of quality of life six dimension health utility d2 disu double vector. The function returns DvD2 (a double vector).
#' @param dvQs_tb DvQs (a tibble)
#' @param kD_1L_dbl KD (a double vector of length one)
#' @param w_dbl W (a double vector)
#' @return DvD2 (a double vector)
#' @rdname calculate_aqol6d_d2_disu_dbl
#' @export 
#' @importFrom purrr pmap_dbl
calculate_aqol6d_d2_disu_dbl <- function (dvQs_tb, kD_1L_dbl, w_dbl) 
{
    dvD2_dbl <- purrr::pmap_dbl(dvQs_tb, ~{
        (1/kD_1L_dbl) * ((1 + (kD_1L_dbl * w_dbl[1] * ..1)) * 
            (1 + (kD_1L_dbl * w_dbl[2] * ..2)) * (1 + (kD_1L_dbl * 
            w_dbl[3] * ..3)) - 1)
    })
    return(dvD2_dbl)
}
#' Calculate Assessment of Quality of Life Six Dimension health utility d3 disu
#' @description calculate_aqol6d_d3_disu_dbl() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate assessment of quality of life six dimension health utility d3 disu double vector. The function returns DvD3 (a double vector).
#' @param dvQs_tb DvQs (a tibble)
#' @param kD_1L_dbl KD (a double vector of length one)
#' @param w_dbl W (a double vector)
#' @return DvD3 (a double vector)
#' @rdname calculate_aqol6d_d3_disu_dbl
#' @export 
#' @importFrom purrr pmap_dbl
calculate_aqol6d_d3_disu_dbl <- function (dvQs_tb, kD_1L_dbl, w_dbl) 
{
    dvD3_dbl <- purrr::pmap_dbl(dvQs_tb, ~{
        (1/kD_1L_dbl) * ((1 + (kD_1L_dbl * w_dbl[1] * ..1)) * 
            (1 + (kD_1L_dbl * w_dbl[2] * ..2)) * (1 + (kD_1L_dbl * 
            w_dbl[3] * ..3)) * (1 + (kD_1L_dbl * w_dbl[4] * 1 * 
            ..4)) - 1)
    })
    return(dvD3_dbl)
}
#' Calculate Assessment of Quality of Life Six Dimension health utility d4 disu
#' @description calculate_aqol6d_d4_disu_dbl() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate assessment of quality of life six dimension health utility d4 disu double vector. The function returns DvD4 (a double vector).
#' @param dvQs_tb DvQs (a tibble)
#' @param kD_1L_dbl KD (a double vector of length one)
#' @param w_dbl W (a double vector)
#' @return DvD4 (a double vector)
#' @rdname calculate_aqol6d_d4_disu_dbl
#' @export 
#' @importFrom purrr pmap_dbl
calculate_aqol6d_d4_disu_dbl <- function (dvQs_tb, kD_1L_dbl, w_dbl) 
{
    dvD4_dbl <- purrr::pmap_dbl(dvQs_tb, ~{
        (1/kD_1L_dbl) * ((1 + (kD_1L_dbl * w_dbl[1] * ..1)) * 
            (1 + (kD_1L_dbl * w_dbl[2] * ..2)) * (1 + (kD_1L_dbl * 
            w_dbl[3] * ..3)) - 1)
    })
    return(dvD4_dbl)
}
#' Calculate Assessment of Quality of Life Six Dimension health utility d5 disu
#' @description calculate_aqol6d_d5_disu_dbl() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate assessment of quality of life six dimension health utility d5 disu double vector. The function returns DvD5 (a double vector).
#' @param dvQs_tb DvQs (a tibble)
#' @param kD_1L_dbl KD (a double vector of length one)
#' @param w_dbl W (a double vector)
#' @return DvD5 (a double vector)
#' @rdname calculate_aqol6d_d5_disu_dbl
#' @export 
#' @importFrom purrr pmap_dbl
calculate_aqol6d_d5_disu_dbl <- function (dvQs_tb, kD_1L_dbl, w_dbl) 
{
    dvD5_dbl <- purrr::pmap_dbl(dvQs_tb, ~{
        (1/kD_1L_dbl) * ((1 + (kD_1L_dbl * w_dbl[1] * ..1)) * 
            (1 + (kD_1L_dbl * w_dbl[2] * ..2)) * (1 + (kD_1L_dbl * 
            w_dbl[3] * ..3)) - 1)
    })
    return(dvD5_dbl)
}
#' Calculate Assessment of Quality of Life Six Dimension health utility d6 disu
#' @description calculate_aqol6d_d6_disu_dbl() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate assessment of quality of life six dimension health utility d6 disu double vector. The function returns DvD6 (a double vector).
#' @param dvQs_tb DvQs (a tibble)
#' @param kD_1L_dbl KD (a double vector of length one)
#' @param w_dbl W (a double vector)
#' @return DvD6 (a double vector)
#' @rdname calculate_aqol6d_d6_disu_dbl
#' @export 
#' @importFrom purrr pmap_dbl
calculate_aqol6d_d6_disu_dbl <- function (dvQs_tb, kD_1L_dbl, w_dbl) 
{
    dvD6_dbl <- purrr::pmap_dbl(dvQs_tb, ~{
        (1/kD_1L_dbl) * ((1 + (kD_1L_dbl * w_dbl[1] * ..1)) * 
            (1 + (kD_1L_dbl * w_dbl[2] * ..2)) * (1 + (kD_1L_dbl * 
            w_dbl[3] * ..3)) - 1)
    })
    return(dvD6_dbl)
}
#' Calculate aqol6dU
#' @description calculate_aqol6dU_dbl() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate aqol6du double vector. The function returns Aqol6dU (a double vector).
#' @param aqol6d_items_tb Assessment of Quality of Life Six Dimension health utility items (a tibble)
#' @param prefix_1L_chr Prefix (a character vector of length one)
#' @param aqol6d_from_8d_coeffs_lup_tb Assessment of Quality of Life Six Dimension health utility from 8d coeffs lookup table (a tibble), Default: aqol6d_from_8d_coeffs_lup_tb
#' @param dim_sclg_constant_lup_tb Dimension sclg constant lookup table (a tibble), Default: dim_sclg_constant_lup_tb
#' @param disvalues_lup_tb Disvalues lookup table (a tibble), Default: disvalues_lup_tb
#' @param itm_wrst_wghts_lup_tb Itm wrst wghts lookup table (a tibble), Default: itm_wrst_wghts_lup_tb
#' @return Aqol6dU (a double vector)
#' @rdname calculate_aqol6dU_dbl
#' @export 
#' @importFrom hutils longest_prefix
calculate_aqol6dU_dbl <- function (aqol6d_items_tb, prefix_1L_chr, aqol6d_from_8d_coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb, 
    dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb, disvalues_lup_tb = disvalues_lup_tb, 
    itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb) 
{
    domains_chr <- dim_sclg_constant_lup_tb$Dimension_chr
    item_pfx_1L_chr <- hutils::longest_prefix(disvalues_lup_tb$Question_chr)
    domain_items_ls <- make_domain_items_ls(domains_chr = domains_chr, 
        q_nbrs_ls = list(1:4, 5:7, 8:11, 12:14, 15:17, 18:20), 
        item_pfx_1L_chr = item_pfx_1L_chr)
    aqol6d_items_tb <- aqol6d_items_tb %>% make_aqol6d_items_tb(old_pfx_1L_chr = prefix_1L_chr, 
        new_pfx_1L_chr = item_pfx_1L_chr) %>% impute_miss_itms_in_aqol6d_items_tb_tb(domain_items_ls = domain_items_ls) %>% 
        add_itm_disu_to_aqol6d_itms_tb_tb(disvalues_lup_tb = disvalues_lup_tb, 
            pfx_1L_chr = item_pfx_1L_chr) %>% add_dmn_disu_to_aqol6d_items_tb_tb(domain_items_ls = domain_items_ls, 
        domains_chr = domains_chr, dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb, 
        itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb) %>% add_dmn_scores_to_aqol6d_items_tb_tb(domain_items_ls = domain_items_ls) %>% 
        add_aqol6dU_to_aqol6d_items_tb_tb(aqol6d_from_8d_coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb)
    aqol6dU_dbl <- aqol6d_items_tb$aqol6dU
    return(aqol6dU_dbl)
}
