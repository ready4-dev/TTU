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
calculate_rmse <- function (y_dbl, yhat_dbl)
{
    rmse_dbl <- sqrt(mean((yhat_dbl - y_dbl)^2))
    return(rmse_dbl)
}
calculate_rmse_tfmn <- function (y_dbl, yhat_dbl)
{
  rmse_tfmn_dbl <- sqrt(mean((1 - exp(-exp(yhat_dbl)) - y_dbl)^2))
  return(rmse_tfmn_dbl)
}
