calculate_depnt_var_tfmn <- function (depnt_var_val_dbl, tfmn_1L_chr = "NTF", tfmn_is_outp_1L_lgl = F, depnt_var_max_val_1L_dbl = NULL)
{
  if(!is.null(depnt_var_max_val_1L_dbl)){
    depnt_var_val_dbl <- depnt_var_val_dbl %>% purrr::map_dbl(~min(.x,depnt_var_max_val_1L_dbl))
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
calculate_rmse <- function (y_dbl, yhat_dbl)
{
    rmse_dbl <- sqrt(mean((yhat_dbl - y_dbl)^2))
    return(rmse_dbl)
}
# calculate_rmse_tfmn <- function (y_dbl, yhat_dbl)
# {
#   rmse_tfmn_dbl <- sqrt(mean((1 - exp(-exp(yhat_dbl)) - y_dbl)^2))
#   return(rmse_tfmn_dbl)
# }
