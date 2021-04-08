#' Calculate dpnt variable transformation
#' @description calculate_dpnt_var_tfmn() is a Calculate function that performs a numeric calculation. Specifically, this function implements an algorithm to calculate dpnt variable transformation. The function returns Transformed dep variable value (a double vector).
#' @param dep_var_val_dbl Dep variable value (a double vector)
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param tfmn_is_outp_1L_lgl Transformation is output (a logical vector of length one), Default: F
#' @return Transformed dep variable value (a double vector)
#' @rdname calculate_dpnt_var_tfmn
#' @export 
#' @importFrom boot inv.logit
#' @importFrom psych logit
#' @keywords internal
calculate_dpnt_var_tfmn <- function (dep_var_val_dbl, tfmn_1L_chr = "NTF", tfmn_is_outp_1L_lgl = F) 
{
    tfd_dep_var_val_dbl <- dep_var_val_dbl
    if (tfmn_1L_chr == "LOG") {
        if (tfmn_is_outp_1L_lgl) 
            tfd_dep_var_val_dbl <- exp(dep_var_val_dbl)
        else tfd_dep_var_val_dbl <- log(dep_var_val_dbl)
    }
    if (tfmn_1L_chr == "LOGIT") {
        if (tfmn_is_outp_1L_lgl) 
            tfd_dep_var_val_dbl <- boot::inv.logit(dep_var_val_dbl)
        else tfd_dep_var_val_dbl <- psych::logit(dep_var_val_dbl)
    }
    if (tfmn_1L_chr == "LOGLOG") {
        if (tfmn_is_outp_1L_lgl) 
            tfd_dep_var_val_dbl <- exp(-exp(-dep_var_val_dbl))
        else tfd_dep_var_val_dbl <- -log(-log(dep_var_val_dbl))
    }
    if (tfmn_1L_chr == "CLL") {
        if (tfmn_is_outp_1L_lgl) 
            tfd_dep_var_val_dbl <- 1 - exp(-exp(dep_var_val_dbl))
        else tfd_dep_var_val_dbl <- log(-log(1 - dep_var_val_dbl))
    }
    return(tfd_dep_var_val_dbl)
}
#' Calculate root mean square error
#' @description calculate_rmse() is a Calculate function that performs a numeric calculation. Specifically, this function implements an algorithm to calculate root mean square error. The function returns Root mean square error (a double vector).
#' @param y_dbl Y (a double vector)
#' @param yhat_dbl Yhat (a double vector)
#' @return Root mean square error (a double vector)
#' @rdname calculate_rmse
#' @export 

#' @keywords internal
calculate_rmse <- function (y_dbl, yhat_dbl) 
{
    rmse_dbl <- sqrt(mean((yhat_dbl - y_dbl)^2))
    return(rmse_dbl)
}
#' Calculate root mean square error transformation
#' @description calculate_rmse_tfmn() is a Calculate function that performs a numeric calculation. Specifically, this function implements an algorithm to calculate root mean square error transformation. The function returns Root mean square error transformation (a double vector).
#' @param y_dbl Y (a double vector)
#' @param yhat_dbl Yhat (a double vector)
#' @return Root mean square error transformation (a double vector)
#' @rdname calculate_rmse_tfmn
#' @export 

#' @keywords internal
calculate_rmse_tfmn <- function (y_dbl, yhat_dbl) 
{
    rmse_tfmn_dbl <- sqrt(mean((1 - exp(-exp(yhat_dbl)) - y_dbl)^2))
    return(rmse_tfmn_dbl)
}
