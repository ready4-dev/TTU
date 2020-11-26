predict_utility <- function(data_tb,
                           tfmn_1L_chr = "NTF",
                           mdl){
  predd_aqol6d_dbl <- predict(mdl, newdata = data_tb) %>% calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, tfmn_is_outp_1L_lgl = T)
  return(predd_aqol6d_dbl)
}
predict_from_mdl_coefs <- function(smry_of_mdl_tb,
                                   new_data_tb){
  coef_tb <- smry_of_mdl_tb %>% dplyr::filter(!Parameter %in% c("R2", "RMSE", "Sigma"))
  vecs_1_ls <- coef_tb$Parameter[-1] %>%
    purrr::map(~ new_data_tb %>%
                 dplyr::pull(.x %>% stringr::str_replace(" ","_")) * coef_tb %>%
                 dplyr::filter(Parameter == .x) %>%
                 dplyr::pull(Estimate)
    )
  pred_dbl <- exp(Reduce(`+`, vecs_1_ls) + coef_tb %>%
                    dplyr::filter(Parameter == "Intercept") %>%
                    dplyr::pull(Estimate))
  return(pred_dbl)
}
