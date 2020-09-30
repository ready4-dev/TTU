reorder_tb_for_target_cors <- function(tbs_ls,
                                       corr_dbl,
                                       corr_var_1L_chr,#aqol6d_total_w
                                       id_var_to_rm_1L_chr = NA_character_#,
                                       # included_fup_idx_dbl = NA_real_
                                       ){
  n_fup_dbl <- nrow(tbs_ls[[2]])
  # if(is.na(included_fup_idx_dbl[1]))
  #   included_fup_idx_dbl <- sample(1:nrow(tbs_ls[[1]]),n_fup_dbl) %>% sort()
  corr_mat <- matrix(corr_dbl, ncol = 2, nrow = 2)
  diag(corr_mat) <- 1
  mvdat_mat <- MASS::mvrnorm(n = n_fup_dbl, mu = c(0, 0), Sigma = corr_mat, empirical = TRUE)
  rank_x_int <- rank(mvdat_mat[ , 1], ties.method = "first")
  rank_y_int <- rank(mvdat_mat[ , 2], ties.method = "first")
  bl_fltd_part_1_tb <- tbs_ls[[1]] %>%
    dplyr::slice(tbs_ls[[2]]$id) %>%
    dplyr::arrange(!!rlang::sym(corr_var_1L_chr)) %>%
    dplyr::slice(rank_x_int)
  tbs_ls[[1]] <- dplyr::bind_rows(bl_fltd_part_1_tb,
                                 tbs_ls[[1]] %>%
                                   dplyr::slice(setdiff(1:nrow(tbs_ls[[1]]),tbs_ls[[2]]$id)))
  tbs_ls[[2]] <- tbs_ls[[2]] %>%
    dplyr::arrange(!!rlang::sym(corr_var_1L_chr)) %>%
    dplyr::slice(rank_y_int)
  if(!is.na(id_var_to_rm_1L_chr))
    tbs_ls <- tbs_ls %>% purrr::map(~.x %>% dplyr::select(-!!rlang::sym(id_var_to_rm_1L_chr)))
  return(tbs_ls)
}
