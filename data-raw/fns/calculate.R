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
calc_hrqol_from_k10_dbl <- function(k10value,
                                    b0.aqolmodel = 0.204665,
                                    b1.aqolmodel = -3.617134,
                                    b0.eq5dmodel = 0.8644649,
                                    b1.eq5dmodel = -2.926161,
                                    aqol_error_dbl = 0,
                                    eq5d_error_dbl = 0){
  meanaqol8dutility<-exp(b0.aqolmodel+b1.aqolmodel*k10value*.01) + aqol_error_dbl
  if(is.na(meanaqol8dutility))
    stop("Mean utility calculation is returning NAs")
  meaneq5dutility<-b0.eq5dmodel+b1.eq5dmodel*(k10value*.01)^2 + eq5d_error_dbl
  if(is.na(meaneq5dutility))
    stop("Mean EQ5D utility calculation is returning NAs")
  return(c(meanaqol8dutility,
           meaneq5dutility))
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
