#' Predict utility
#' @description predict_utility() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict utility. The function returns Predd aqol6d (a double vector).
#' @param data_tb Data (a tibble)
#' @param tfmn_1L_chr Tfmn (a character vector of length one), Default: 'NTF'
#' @param ... Additional arguments
#' @return Predd aqol6d (a double vector)
#' @rdname predict_utility
#' @export 

#' @keywords internal
predict_utility <- function (data_tb, tfmn_1L_chr = "NTF", mdl) 
{
    predd_aqol6d_dbl <- predict(mdl, newdata = data_tb) %>% calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, 
        tfmn_is_outp_1L_lgl = T)
    return(predd_aqol6d_dbl)
}
