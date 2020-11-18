#install.packages("synthpop")
library(magrittr)
library(FBaqol)

mdl_2 <- make_shareable_mdl(data_tb = PHQ9_SOFAS_1_CLL_TFN$data %>%
                              dplyr::rename(aqol6d_total_w_CLL = aqol6d_total_w_cloglog),
                            mdl_smry_tb = test_smry)
predict(mdl_2) %>% calculate_dep_var_tfmn(tfmn_1L_chr = "CLL", tfmn_is_outp_1L_lgl = T) %>% summary()
PHQ9_SOFAS_1_CLL_TFN$data$aqol6d_total_w_cloglog %>% calculate_dep_var_tfmn(tfmn_1L_chr = "CLL", tfmn_is_outp_1L_lgl = T) %>% summary()
