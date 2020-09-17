calculate_aqol6d_d1_disu_dbl <- function(dvQs_tb,
                                    kD_1L_dbl,
                                    w_dbl){
  dvD1_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))*(1+(kD_1L_dbl*w_dbl[4]*..4))-1)
                              })
  return(dvD1_dbl)
}
calculate_aqol6d_d2_disu_dbl <- function(dvQs_tb,
                                    kD_1L_dbl,
                                    w_dbl){
  dvD2_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))-1)
                              })
  return(dvD2_dbl)
}
calculate_aqol6d_d3_disu_dbl <- function(dvQs_tb,
                                    kD_1L_dbl,
                                    w_dbl){
  dvD3_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))*(1+(kD_1L_dbl*w_dbl[4]*1*..4))-1)
                              })
  return(dvD3_dbl)
}
calculate_aqol6d_d4_disu_dbl <- function(dvQs_tb,
                                    kD_1L_dbl,
                                    w_dbl){
  dvD4_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))-1)
                              })
  return(dvD4_dbl)
}
calculate_aqol6d_d5_disu_dbl <- function(dvQs_tb,
                                    kD_1L_dbl,
                                    w_dbl){
  dvD5_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))-1)
                              })
  return(dvD5_dbl)
}
calculate_aqol6d_d6_disu_dbl <- function(dvQs_tb,
                                    kD_1L_dbl,
                                    w_dbl){
  dvD6_dbl <- purrr::pmap_dbl(dvQs_tb,
                              ~{
                                (1/kD_1L_dbl)*((1+(kD_1L_dbl*w_dbl[1]*..1))*(1+(kD_1L_dbl*w_dbl[2]*..2))*(1+(kD_1L_dbl*w_dbl[3]*..3))-1)
                              })
  return(dvD6_dbl)
}
calculate_aqol6dU_dbl <- function(aqol6d_items_tb,
                             prefix_1L_chr,
                             aqol6d_from_8d_coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb,
                             dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb,
                             disutilities_lup_tb = disutilities_lup_tb,
                             itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb){
  domains_chr <- dim_sclg_constant_lup_tb$Dimension_chr
  item_pfx_1L_chr <- hutils::longest_prefix(disutilities_lup_tb$Question_chr)
  domain_items_ls <- make_domain_items_ls(domains_chr = domains_chr,
                                          q_nbrs_ls = list(1:4,5:7,8:11,12:14,15:17,18:20),
                                          item_pfx_1L_chr = item_pfx_1L_chr)
  aqol6d_items_tb <- aqol6d_items_tb %>% make_aqol6d_items_tb(old_pfx_1L_chr = prefix_1L_chr,
                                                              new_pfx_1L_chr = item_pfx_1L_chr) %>%
    impute_miss_itms_in_aqol6d_items_tb_tb(domain_items_ls = domain_items_ls) %>%
    add_itm_disu_to_aqol6d_itms_tb_tb(disutilities_lup_tb = disutilities_lup_tb,
                                      pfx_1L_chr = item_pfx_1L_chr) %>%
    add_dmn_disu_to_aqol6d_items_tb_tb(domain_items_ls = domain_items_ls,
                                       domains_chr = domains_chr,
                                       dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb,
                                       itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb) %>%
    add_dmn_scores_to_aqol6d_items_tb_tb(domain_items_ls = domain_items_ls) %>%
    add_aqol6dU_to_aqol6d_items_tb_tb(aqol6d_from_8d_coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb)
  aqol6dU_dbl <- aqol6d_items_tb$aqol6dU
  return(aqol6dU_dbl)
}
add_scrg_eqs <- function(equations_tb,unscored_aqol_tb){
  for(var in equations_tb$Dim_scal) {
    expression=equations_tb[equations_tb$Dim_scal==var,]$Equ
    unscored_aqol_tb <- unscored_aqol_tb %>%
      mutate(!! var := !! parse_expr(expression))
    Hmisc::label(unscored_aqol_tb[,var])=equations_tb[equations_tb$Dim_scal==var,]$Label
  }
  return(unscored_aqol_tb)
}
impute_adol_unscrd_aqol_ds <- function(unscrd_aqol_ds){
  unscrd_aqol_ds<-unscrd_aqol_ds %>%
    dplyr::mutate(missing=rowSums(is.na(dplyr::select(.,
                                                      paste0("Q",
                                                             c(1:10))))))
  aqol_cases_to_imp_tb <- unscrd_aqol_ds%>%
    dplyr::filter(missing < 10) %>%
    dplyr::select(-missing)
  aqol_cases_not_to_imp_tb <- unscrd_aqol_ds%>%
    dplyr::filter(missing >= 10) %>%
    dplyr::select(-missing)
  imputed_aqol_tb <- mice::mice(aqol_cases_to_imp_tb, m=1, maxit=50, meth='pmm', seed=1234)
  aqol_cases_to_imp_tb <- mice::complete(imputed_aqol_tb, 'long') %>%
    dplyr::select(-.imp,- .id)
  imputed_unscrd_aqol_ds_tb <- data.frame(rbind(aqol_cases_to_imp_tb,
                                             aqol_cases_not_to_imp_tb))
  return(imputed_unscrd_aqol_ds_tb)
}

calc_adol_aqol6dU <- function(unscored_aqol_tb,
                              prefix_1L_chr =  "aqol",
                              id_var_nm_1L_chr){
  #IDname <- names(unscored_aqol_tb)[1]
  #rename variable
  unscored_aqol_tb <- unscored_aqol_tb %>%
    dplyr::select(id_var_nm_1L_chr,
                  dplyr::starts_with(prefix_1L_chr))
  names(unscored_aqol_tb) <- c("ID",paste0("Q",1:20))
  unscored_aqol_tb <- impute_adol_unscrd_aqol_ds(unscored_aqol_tb)
  # lookup <- read.csv("data-raw/AQoL_mapping.csv") #read.csv(here::here("Data cleaning","Data","AQoL_mapping.csv"))
  # names(lookup)[1] <- "Variable"
  # dvQ <- unscored_aqol_tb %>%
  #   tidyr::gather(key = "Variable",value="Value",-ID)%>%
  #   dplyr::left_join(lookup, by =c("Variable","Value") )%>%
  #   dplyr::select(-Value) %>%
  #   tidyr::spread(key = Variable, value = Scoring) %>%
  #   dplyr::select(ID,paste0("Q",1:20))
  disvals_tb <- unscored_aqol_tb %>%
    add_itm_disu_to_aqol6d_itms_tb_tb(disutilities_lup_tb = make_adol_aqol6du_disu_lup(),
                                    pfx_1L_chr = "Q") %>%
    dplyr::select(ID,
                  dplyr::starts_with("dv_")) %>%
    dplyr::rename_all(~stringr::str_replace(.x,"dv_","dv"))
  #names(dvQ) <- c("ID",paste0("dvQ",1:20))
  equations_tb <- read.csv("data-raw/AQoL_6D_Dim_Scaling.csv", stringsAsFactors = F)
  names(equations_tb)[1] <- "Dim_scal"
    #read.csv(here::here("Data cleaning","Data","AQoL_6D_Dim_Scaling.csv"),stringsAsFactors=F)
  scored_aqol_tb <- add_scrg_eqs(equations_tb,disvals_tb)
  #names(scored_aqol_tb)[1] <- IDname
  return(scored_aqol_tb)
}
make_adol_aqol6du_disu_lup <- function(){
  data("disutilities_lup_tb", package = "FBaqol", envir = environment())
  adol_aqol6du_dist_lup <- disutilities_lup_tb %>%
    dplyr::mutate(Answer_4_dbl = dplyr::case_when(Question_chr == "Q18" ~ 0.622,
                                              TRUE ~ Answer_4_dbl),
                  Answer_5_dbl = dplyr::case_when(Question_chr == "Q3" ~ 0.827,
                                              TRUE ~ Answer_5_dbl),
                  Answer_6_dbl = dplyr::case_when(Question_chr == "Q1" ~ 0.073,
                                              TRUE ~ Answer_5_dbl))
  return(adol_aqol6du_dist_lup)
}

