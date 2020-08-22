impute_miss_itms_in_aqol6d_items_tb_tb <- function(aqol6d_items_tb,
                                                   domain_items_ls){
  aqol6d_items_tb <- 1:length(domain_items_ls) %>%
    purrr::reduce(.init = aqol6d_items_tb,
                  ~ {
                    idx_int <- .y
                    tb <- .x %>% dplyr::mutate(!!rlang::sym(paste0("D",names(domain_items_ls)[idx_int],"miss")) := rowMeans(dplyr::select(., domain_items_ls[[idx_int]]),na.rm = TRUE))

                    tb %>% dplyr::mutate_at(dplyr::vars(domain_items_ls[[idx_int]]), ~ .x %>% purrr::map2_dbl(!!rlang::sym(paste0("D",names(domain_items_ls)[idx_int],"miss")),
                                                                                                              ~{
                                                                                                                ifelse(is.na(.x),.y,.x)
                                                                                                              }))
                    tb
                  })
  return(aqol6d_items_tb)
}
