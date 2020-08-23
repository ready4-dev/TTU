#' Impute miss itms in aqol6d items tibble
#' @description impute_miss_itms_in_aqol6d_items_tb_tb() is an Impute function that imputes data. Specifically, this function implements an algorithm to impute miss itms in aqol6d items a tibble. The function returns aqol6d items (a tibble).
#' @param aqol6d_items_tb Aqol6d items (a tibble)
#' @param domain_items_ls Domain items (a list)
#' @return Aqol6d items (a tibble)
#' @rdname impute_miss_itms_in_aqol6d_items_tb_tb
#' @export 
#' @importFrom purrr reduce map2_dbl
#' @importFrom dplyr mutate select mutate_at vars
#' @importFrom rlang sym
#' @keywords internal
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
