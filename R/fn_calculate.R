#' Calculate dep var transformation
#' @description calculate_dep_var_tfmn() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate dep var transformation. Function argument dep_var_val_dbl specifies the calculate The function returns Transformed dep var value (a double vector).
#' @param dep_var_val_dbl Dep var value (a double vector)
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param tfmn_is_outp_1L_lgl Transformation is output (a logical vector of length one), Default: F
#' @return Transformed dep var value (a double vector)
#' @rdname calculate_dep_var_tfmn
#' @export 
#' @importFrom boot inv.logit
#' @importFrom psych logit
calculate_dep_var_tfmn <- function (dep_var_val_dbl, tfmn_1L_chr = "NTF", tfmn_is_outp_1L_lgl = F) 
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
#' Calculate rmse
#' @description calculate_rmse() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate rmse. Function argument y_dbl specifies the calculate The function returns Rmse (a double vector).
#' @param y_dbl Y (a double vector)
#' @param yhat_dbl Yhat (a double vector)
#' @return Rmse (a double vector)
#' @rdname calculate_rmse
#' @export 

calculate_rmse <- function (y_dbl, yhat_dbl) 
{
    rmse_dbl <- sqrt(mean((yhat_dbl - y_dbl)^2))
    return(rmse_dbl)
}
