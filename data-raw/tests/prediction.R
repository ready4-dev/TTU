#install.packages("synthpop")
library(magrittr)
library(FBaqol)

## 1. REMAKE STARTER DS WITH CORRELATED SOFAS
if(is.null(mdl_types_lup))
  data(mdl_types_lup, envir = environment())
shrble_ols_cll_mdl <- make_shareable_mdl(data_tb = PHQ9_SOFAS_1_CLL_TFN$data %>%
                              dplyr::rename(aqol6d_total_w_CLL = aqol6d_total_w_cloglog),
                            mdl_smry_tb = test_smry)
predict_aqol6d <- function(data_tb,
                           tfmn_1L_chr = "NTF",
                           mdl){
  predd_aqol6d_dbl <- predict(mdl) %>% calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, tfmn_is_outp_1L_lgl = T)
  return(predd_aqol6d_dbl)
}
 #%>% summary()
PHQ9_SOFAS_1_CLL_TFN$data$aqol6d_total_w_cloglog %>% calculate_dep_var_tfmn(tfmn_1L_chr = "CLL", tfmn_is_outp_1L_lgl = T) %>% summary()
