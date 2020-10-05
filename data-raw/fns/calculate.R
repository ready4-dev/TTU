calculate_adol_aqol6dU <- function(unscored_aqol_tb,
                                  prefix_1L_chr =  "aqol",
                                  id_var_nm_1L_chr){
  unscored_aqol_tb <- unscored_aqol_tb %>%
    dplyr::select(id_var_nm_1L_chr,
                  dplyr::starts_with(unname(prefix_1L_chr)))
  names(unscored_aqol_tb) <- c("ID",paste0("Q",1:20))
  unscored_aqol_tb <- impute_unscrd_adol_aqol6d_ds(unscored_aqol_tb)
  disvals_tb <- unscored_aqol_tb %>%
    add_itm_disv_to_aqol6d_itms_tb(disvalues_lup_tb = make_adol_aqol6d_disv_lup(),
                                      pfx_1L_chr = "Q") %>%
    dplyr::select(ID,
                  dplyr::starts_with("dv_")) %>%
    dplyr::rename_all(~stringr::str_replace(.x,"dv_","dv"))
  scored_aqol_tb <- add_aqol6d_adol_dim_scrg_eqs(disvals_tb)
  adol_aqol6d_dbl <- scored_aqol_tb$uaqol
  return(adol_aqol6d_dbl)
}
calculate_adult_aqol6dU <- function(aqol6d_items_tb,
                              prefix_1L_chr,
                              coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb,
                              dim_sclg_con_lup_tb = aqol6d_dim_sclg_con_lup_tb,
                              disvalues_lup_tb = aqol6d_adult_disv_lup_tb,
                              itm_wrst_wghts_lup_tb = aqol6d_adult_itm_wrst_wghts_lup_tb){
  domains_chr <- dim_sclg_con_lup_tb$Dimension_chr
  item_pfx_1L_chr <- hutils::longest_prefix(disvalues_lup_tb$Question_chr)
  domain_items_ls <- make_domain_items_ls(domain_qs_lup_tb = aqol6d_domain_qs_lup_tb,
                                          item_pfx_1L_chr = item_pfx_1L_chr)
  aqol6d_items_tb <- aqol6d_items_tb %>% make_aqol6d_items_tb(old_pfx_1L_chr = prefix_1L_chr,
                                                              new_pfx_1L_chr = item_pfx_1L_chr) %>%
    impute_adult_aqol6d_items_tb(domain_items_ls = domain_items_ls) %>%
    add_itm_disv_to_aqol6d_itms_tb(disvalues_lup_tb = disvalues_lup_tb,
                                   pfx_1L_chr = item_pfx_1L_chr) %>%
    add_dim_disv_to_aqol6d_items_tb(domain_items_ls = domain_items_ls,
                                    domains_chr = domains_chr,
                                    dim_sclg_con_lup_tb = dim_sclg_con_lup_tb,
                                    itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb) %>%
    add_dim_scores_to_aqol6d_items_tb(domain_items_ls = domain_items_ls) %>%
    add_aqol6dU_to_aqol6d_items_tb(coeffs_lup_tb = coeffs_lup_tb)
  aqol6dU_dbl <- aqol6d_items_tb$aqol6dU
  return(aqol6dU_dbl)
}
calculate_aqol6d_dim_1_disv <- function(dvQs_tb,
                                        kD_1L_dbl,
                                        w_dbl){
  dvD1_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))*(1+(kD_1L_dbl*w_dbl[4]*..4))-1)
                              })
  return(dvD1_dbl)
}
calculate_aqol6d_dim_2_disv <- function(dvQs_tb,
                                        kD_1L_dbl,
                                        w_dbl){
  dvD2_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))-1)
                              })
  return(dvD2_dbl)
}
calculate_aqol6d_dim_3_disv <- function(dvQs_tb,
                                        kD_1L_dbl,
                                        w_dbl){
  dvD3_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))*(1+(kD_1L_dbl*w_dbl[4]*1*..4))-1)
                              })
  return(dvD3_dbl)
}
calculate_aqol6d_dim_4_disv <- function(dvQs_tb,
                                        kD_1L_dbl,
                                        w_dbl){
  dvD4_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))-1)
                              })
  return(dvD4_dbl)
}
calculate_aqol6d_dim_5_disv <- function(dvQs_tb,
                                        kD_1L_dbl,
                                        w_dbl){
  dvD5_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))-1)
                              })
  return(dvD5_dbl)
}
calculate_aqol6d_dim_6_disv <- function(dvQs_tb,
                                        kD_1L_dbl,
                                        w_dbl){
  dvD6_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))-1)
                              })
  return(dvD6_dbl)
}




