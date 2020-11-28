predict_utility <- function (data_tb, tfmn_1L_chr = "NTF", model_mdl)
{
    predd_aqol6d_dbl <- predict(model_mdl, newdata = data_tb) %>% calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
        tfmn_is_outp_1L_lgl = T)
    return(predd_aqol6d_dbl)
}
