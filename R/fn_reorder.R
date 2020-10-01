#' Reorder tibbles for target cors
#' @description reorder_tbs_for_target_cors() is a Reorder function that reorders an object to conform to a pre-specified schema. Specifically, this function implements an algorithm to reorder tibbles for target cors. The function returns Tibbles (a list).
#' @param tbs_ls Tibbles (a list)
#' @param corr_dbl Corr (a double vector)
#' @param corr_var_chr Corr var (a character vector)
#' @param id_var_to_rm_1L_chr Id var to rm (a character vector of length one), Default: 'NA'
#' @return Tibbles (a list)
#' @rdname reorder_tbs_for_target_cors
#' @export 
#' @importFrom MASS mvrnorm
#' @importFrom dplyr slice arrange bind_rows select
#' @importFrom rlang sym
#' @importFrom purrr map
reorder_tbs_for_target_cors <- function (tbs_ls, corr_dbl, corr_var_chr, id_var_to_rm_1L_chr = NA_character_) 
{
    n_fup_dbl <- nrow(tbs_ls[[2]])
    corr_mat <- matrix(corr_dbl, ncol = 2, nrow = 2)
    diag(corr_mat) <- 1
    mvdat_mat <- MASS::mvrnorm(n = n_fup_dbl, mu = c(0, 0), Sigma = corr_mat, 
        empirical = TRUE)
    rank_x_int <- rank(mvdat_mat[, 1], ties.method = "first")
    rank_y_int <- rank(mvdat_mat[, 2], ties.method = "first")
    matched_tb <- tbs_ls[[1]] %>% dplyr::slice(tbs_ls[[2]]$id) %>% 
        dplyr::arrange(!!rlang::sym(corr_var_chr[1])) %>% dplyr::slice(rank_x_int)
    unmatched_int <- setdiff(1:nrow(tbs_ls[[1]]), tbs_ls[[2]]$id)
    if (!identical(unmatched_int, character(0))) {
        tbs_ls[[1]] <- dplyr::bind_rows(matched_tb, tbs_ls[[1]] %>% 
            dplyr::slice(unmatched_int))
    }
    else {
        tbs_ls[[1]] <- matched_tb
    }
    tbs_ls[[2]] <- tbs_ls[[2]] %>% dplyr::arrange(!!rlang::sym(corr_var_chr[2])) %>% 
        dplyr::slice(rank_y_int)
    if (!is.na(id_var_to_rm_1L_chr)) 
        tbs_ls <- tbs_ls %>% purrr::map(~.x %>% dplyr::select(-!!rlang::sym(id_var_to_rm_1L_chr)))
    return(tbs_ls)
}
