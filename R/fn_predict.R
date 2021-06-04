#' Predict utility
#' @description predict_utility() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict utility. The function returns Predicted utility (a double vector).
#' @param data_tb Data (a tibble)
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param model_mdl Model (a model)
#' @param force_min_max_1L_lgl Force minimum maximum (a logical vector of length one), Default: T
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: 0.03
#' @param impute_1L_lgl Impute (a logical vector of length one), Default: T
#' @param utl_cls_fn Utility class (a function), Default: NULL
#' @return Predicted utility (a double vector)
#' @rdname predict_utility
#' @export 
#' @importFrom rlang exec
predict_utility <- function (data_tb, tfmn_1L_chr = "NTF", model_mdl, force_min_max_1L_lgl = T, 
    utl_min_val_1L_dbl = 0.03, impute_1L_lgl = T, utl_cls_fn = NULL) 
{
    predd_utl_dbl <- predict(model_mdl, newdata = data_tb) %>% 
        calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, tfmn_is_outp_1L_lgl = T) %>% 
        as.vector()
    if (impute_1L_lgl) 
        predd_utl_dbl[which(is.na(predd_utl_dbl))] <- predd_utl_dbl %>% 
            na.omit() %>% mean()
    if (force_min_max_1L_lgl) {
        predd_utl_dbl[which(predd_utl_dbl > 1)] <- 1
        predd_utl_dbl[which(predd_utl_dbl < utl_min_val_1L_dbl)] <- utl_min_val_1L_dbl
    }
    if (!is.null(utl_cls_fn)) {
        predd_utl_dbl <- predd_utl_dbl %>% rlang::exec(.fn = utl_cls_fn)
    }
    return(predd_utl_dbl)
}
#' Predict utility from k10
#' @description predict_utl_from_k10() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict utility from k10. The function is called for its side effects and does not return a value.
#' @param k10_1L_dbl K10 (a double vector of length one)
#' @param b0_aqol_mdl_1L_dbl B0 Assessment of Quality of Life model (a double vector of length one), Default: 0.204665
#' @param b1_aqol_mdl_1L_dbl B1 Assessment of Quality of Life model (a double vector of length one), Default: -3.617134
#' @param b0_eq5d_mdl_1L_dbl B0 eq5d model (a double vector of length one), Default: 0.8644649
#' @param b1_eq5d_mdl_1L_dbl B1 eq5d model (a double vector of length one), Default: -2.926161
#' @param aqol_error_1L_dbl Assessment of Quality of Life error (a double vector of length one), Default: 0
#' @param eq5d_error_1L_dbl Eq5d error (a double vector of length one), Default: 0
#' @return NA ()
#' @rdname predict_utl_from_k10
#' @export 

predict_utl_from_k10 <- function (k10_1L_dbl, b0_aqol_mdl_1L_dbl = 0.204665, b1_aqol_mdl_1L_dbl = -3.617134, 
    b0_eq5d_mdl_1L_dbl = 0.8644649, b1_eq5d_mdl_1L_dbl = -2.926161, 
    aqol_error_1L_dbl = 0, eq5d_error_1L_dbl = 0) 
{
    meanaqol8dutility <- exp(b0_aqol_mdl_1L_dbl + b1_aqol_mdl_1L_dbl * 
        k10_1L_dbl * 0.01) + aqol_error_1L_dbl
    if (is.na(meanaqol8dutility)) 
        stop("Mean utility calculation is returning NAs")
    meaneq5dutility <- b0_eq5d_mdl_1L_dbl + b1_eq5d_mdl_1L_dbl * 
        (k10_1L_dbl * 0.01)^2 + eq5d_error_1L_dbl
    if (is.na(meaneq5dutility)) 
        stop("Mean EQ5D utility calculation is returning NAs")
    return(c(meanaqol8dutility, meaneq5dutility))
}
