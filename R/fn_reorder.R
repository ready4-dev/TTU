#' Reorder tibble for target cors
#' @description reorder_tb_for_target_cors() is a Reorder function that reorders an object to conform to a pre-specified schema. Specifically, this function implements an algorithm to reorder a tibble for target cors. The function returns a tibbles (a list).
#' @param tbs_ls Tibbles (a list)
#' @param corr_dbl Corr (a double vector of length 1)
#' @param corr_var_1L_chr Corr var 1L (a character vector of length 1)
#' @param id_var_to_rm_1L_chr Id var to rm 1L (a character vector of length 1), Default: 'NA'
#' @param included_fup_idx_dbl Included fup idx (a double vector of length 1), Default: NA
#' @return Tibbles (a list)
#' @rdname reorder_tb_for_target_cors
#' @export 
#' @importFrom MASS mvrnorm
#' @importFrom dplyr slice arrange bind_rows select
#' @importFrom rlang sym
#' @importFrom purrr map
#' @keywords internal
reorder_tb_for_target_cors <- function (tbs_ls, corr_dbl, corr_var_1L_chr, id_var_to_rm_1L_chr = NA_character_, 
    included_fup_idx_dbl = NA_real_) 
{
    n_fup_dbl <- nrow(tbs_ls[[2]])
    if (is.na(included_fup_idx_dbl[1])) 
        included_fup_idx_dbl <- sample(1:nrow(tbs_ls[[1]]), n_fup_dbl) %>% 
            sort()
    corr_mat <- matrix(corr_dbl, ncol = 2, nrow = 2)
    diag(corr_mat) <- 1
    mvdat_mat <- MASS::mvrnorm(n = n_fup_dbl, mu = c(0, 0), Sigma = corr_mat, 
        empirical = TRUE)
    rank_x_int <- rank(mvdat_mat[, 1], ties.method = "first")
    rank_y_int <- rank(mvdat_mat[, 2], ties.method = "first")
    bl_fltd_part_1_tb <- tbs_ls[[1]] %>% dplyr::slice(included_fup_idx_dbl) %>% 
        dplyr::arrange(!!rlang::sym(corr_var_1L_chr)) %>% dplyr::slice(rank_x_int)
    tbs_ls[[1]] <- dplyr::bind_rows(bl_fltd_part_1_tb, tbs_ls[[1]] %>% 
        dplyr::slice(setdiff(1:nrow(tbs_ls[[1]]), included_fup_idx_dbl)))
    tbs_ls[[2]] <- tbs_ls[[1]] %>% dplyr::arrange(!!rlang::sym(corr_var_1L_chr)) %>% 
        dplyr::slice(rank_y_int)
    if (!is.na(id_var_to_rm_1L_chr)) 
        tbs_ls <- tbs_ls %>% purrr::map(~.x %>% dplyr::select(-!!rlang::sym(id_var_to_rm_1L_chr)))
    return(tbs_ls)
}
