#' Impute adolescent unscored Assessment of Quality of Life health utility dataset
#' @description impute_adol_unscrd_aqol_ds() is an Impute function that imputes data. Specifically, this function implements an algorithm to impute adolescent unscored assessment of quality of life health utility dataset. The function returns Imputed unscored Assessment of Quality of Life health utility dataset tibble (a tibble).
#' @param unscrd_aqol_ds_tb Unscored Assessment of Quality of Life health utility dataset (a tibble)
#' @return Imputed unscored Assessment of Quality of Life health utility dataset tibble (a tibble)
#' @rdname impute_adol_unscrd_aqol_ds
#' @export 
#' @importFrom dplyr mutate select filter
#' @importFrom mice mice complete
impute_adol_unscrd_aqol_ds <- function (unscrd_aqol_ds_tb) 
{
    unscrd_aqol_ds_tb <- unscrd_aqol_ds_tb %>% dplyr::mutate(missing = rowSums(is.na(dplyr::select(., 
        paste0("Q", c(1:10))))))
    aqol_cases_to_imp_tb <- unscrd_aqol_ds_tb %>% dplyr::filter(missing < 
        10) %>% dplyr::select(-missing)
    aqol_cases_not_to_imp_tb <- unscrd_aqol_ds_tb %>% dplyr::filter(missing >= 
        10) %>% dplyr::select(-missing)
    imputed_aqol_tb <- mice::mice(aqol_cases_to_imp_tb, m = 1, 
        maxit = 50, meth = "pmm", seed = 1234)
    aqol_cases_to_imp_tb <- mice::complete(imputed_aqol_tb, "long") %>% 
        dplyr::select(-.imp, -.id)
    imputed_unscrd_aqol_ds_tb_tb <- data.frame(rbind(aqol_cases_to_imp_tb, 
        aqol_cases_not_to_imp_tb))
    return(imputed_unscrd_aqol_ds_tb_tb)
}
#' Impute miss itms in Assessment of Quality of Life Six Dimension health utility items tibble
#' @description impute_miss_itms_in_aqol6d_items_tb_tb() is an Impute function that imputes data. Specifically, this function implements an algorithm to impute miss itms in assessment of quality of life six dimension health utility items tibble tibble. The function returns an Assessment of Quality of Life Six Dimension health utility items (a tibble).
#' @param aqol6d_items_tb Assessment of Quality of Life Six Dimension health utility items (a tibble)
#' @param domain_items_ls Domain items (a list)
#' @return an Assessment of Quality of Life Six Dimension health utility items (a tibble)
#' @rdname impute_miss_itms_in_aqol6d_items_tb_tb
#' @export 
#' @importFrom purrr reduce map2_dbl
#' @importFrom dplyr mutate select mutate_at vars
#' @importFrom rlang sym
impute_miss_itms_in_aqol6d_items_tb_tb <- function (aqol6d_items_tb, domain_items_ls) 
{
    aqol6d_items_tb <- 1:length(domain_items_ls) %>% purrr::reduce(.init = aqol6d_items_tb, 
        ~{
            idx_int <- .y
            tb <- .x %>% dplyr::mutate(`:=`(!!rlang::sym(paste0("D", 
                names(domain_items_ls)[idx_int], "miss")), rowMeans(dplyr::select(., 
                domain_items_ls[[idx_int]]), na.rm = TRUE)))
            tb %>% dplyr::mutate_at(dplyr::vars(domain_items_ls[[idx_int]]), 
                ~.x %>% purrr::map2_dbl(!!rlang::sym(paste0("D", 
                  names(domain_items_ls)[idx_int], "miss")), 
                  ~{
                    ifelse(is.na(.x), .y, .x)
                  }))
            tb
        })
    return(aqol6d_items_tb)
}
