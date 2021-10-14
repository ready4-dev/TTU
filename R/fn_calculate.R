#' Calculate dependent variable transformation
#' @description calculate_depnt_var_tfmn() is a Calculate function that performs a numeric calculation. Specifically, this function implements an algorithm to calculate dependent variable transformation. The function returns Transformed dependent variable value (a double vector).
#' @param depnt_var_val_dbl Dependent variable value (a double vector)
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param tfmn_is_outp_1L_lgl Transformation is output (a logical vector of length one), Default: F
#' @param depnt_var_max_val_1L_dbl Dependent variable maximum value (a double vector of length one), Default: NULL
#' @return Transformed dependent variable value (a double vector)
#' @rdname calculate_depnt_var_tfmn
#' @export 
#' @importFrom purrr map_dbl
#' @importFrom boot inv.logit
#' @importFrom psych logit
#' @keywords internal
calculate_depnt_var_tfmn <- function (depnt_var_val_dbl, tfmn_1L_chr = "NTF", tfmn_is_outp_1L_lgl = F, 
    depnt_var_max_val_1L_dbl = NULL) 
{
    if (!is.null(depnt_var_max_val_1L_dbl)) {
        depnt_var_val_dbl <- depnt_var_val_dbl %>% purrr::map_dbl(~min(.x, 
            depnt_var_max_val_1L_dbl))
    }
    tfd_depnt_var_val_dbl <- depnt_var_val_dbl
    if (tfmn_1L_chr == "LOG") {
        if (tfmn_is_outp_1L_lgl) 
            tfd_depnt_var_val_dbl <- exp(depnt_var_val_dbl)
        else tfd_depnt_var_val_dbl <- log(depnt_var_val_dbl)
    }
    if (tfmn_1L_chr == "LOGIT") {
        if (tfmn_is_outp_1L_lgl) 
            tfd_depnt_var_val_dbl <- boot::inv.logit(depnt_var_val_dbl)
        else tfd_depnt_var_val_dbl <- psych::logit(depnt_var_val_dbl)
    }
    if (tfmn_1L_chr == "LOGLOG") {
        if (tfmn_is_outp_1L_lgl) 
            tfd_depnt_var_val_dbl <- exp(-exp(-depnt_var_val_dbl))
        else tfd_depnt_var_val_dbl <- -log(-log(depnt_var_val_dbl))
    }
    if (tfmn_1L_chr == "CLL") {
        if (tfmn_is_outp_1L_lgl) 
            tfd_depnt_var_val_dbl <- 1 - exp(-exp(depnt_var_val_dbl))
        else tfd_depnt_var_val_dbl <- log(-log(1 - depnt_var_val_dbl))
    }
    return(tfd_depnt_var_val_dbl)
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
