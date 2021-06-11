predict_utl_from_k10 <- function(k10_1L_dbl,
                                     b0_aqol_mdl_1L_dbl = 0.204665,
                                     b1_aqol_mdl_1L_dbl = -3.617134,
                                     b0_eq5d_mdl_1L_dbl = 0.8644649,
                                     b1_eq5d_mdl_1L_dbl = -2.926161,
                                     aqol_error_1L_dbl = 0,
                                     eq5d_error_1L_dbl = 0){
  meanaqol8dutility<-exp(b0_aqol_mdl_1L_dbl+b1_aqol_mdl_1L_dbl*k10_1L_dbl*.01) + aqol_error_1L_dbl
  if(is.na(meanaqol8dutility))
    stop("Mean utility calculation is returning NAs")
  meaneq5dutility<-b0_eq5d_mdl_1L_dbl+b1_eq5d_mdl_1L_dbl*(k10_1L_dbl*.01)^2 + eq5d_error_1L_dbl
  if(is.na(meaneq5dutility))
    stop("Mean EQ5D utility calculation is returning NAs")
  return(c(meanaqol8dutility,
           meaneq5dutility))
}
predict_uncnstrd_utl <- function(data_tb, model_mdl,
                                new_data_is_1L_chr = "Predicted",
                                predn_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_, tfmn_1L_chr = "NTF", is_brms_mdl_1L_lgl = F){
  if(new_data_is_1L_chr == "Predicted")
    new_data_dbl <- stats::predict(model_mdl, type = predn_type_1L_chr)
  if(new_data_is_1L_chr == "Simulated"){
    if(is_brms_mdl_1L_lgl){
      new_data_dbl <- brms::posterior_predict(model_mdl,
                                              newdata = data_tb,
                                              nsamples=1) %>% as.vector()

    }else{
      if("betareg" %in% class(model_mdl)){
        new_data_dbl <- rlang::exec(enrichwith::get_simulate_function(model_mdl),
                                    coef(enrichwith::enrich(model_mdl, with = "auxiliary functions")))
      }else{
        if (!tfmn_for_bnml_1L_lgl){
          new_data_dbl <- stats::simulate(model_mdl)$sim_1
        }else{
          new_data_dbl <- (stats::predict(model_mdl) + stats::rnorm(nrow(data_tb),
                                                                    0, stats::sigma(model_mdl)))
        }
      }
    }
  }
  if(is.matrix(new_data_dbl))
    new_data_dbl <- new_data_dbl[, 1]
  new_data_dbl <- new_data_dbl %>% # Make CNDL ON BRMS???
    calculate_dpnt_var_tfmn(tfmn_1L_chr = ifelse(tfmn_for_bnml_1L_lgl & new_data_is_1L_chr == "Simulated",
                                                 ifelse(family_1L_chr == "quasibinomial(log)",
                                                        "LOG",
                                                        ifelse(family_1L_chr == "quasibinomial(logit)",
                                                               "LOGIT",
                                                               ifelse(family_1L_chr == "quasibinomial(cloglog)",
                                                                      "CLL", "NTF"))),
                                                 tfmn_1L_chr),
                            tfmn_is_outp_1L_lgl = T)
  return_new_data_dbl
}
predict_utility <- function (data_tb,
                             tfmn_1L_chr = "NTF",
                             model_mdl,
                             force_min_max_1L_lgl = T,
                             utl_min_val_1L_dbl = 0.03,
                             impute_1L_lgl = T,
                             utl_cls_fn = NULL,
                             new_data_is_1L_chr = "Predicted",
                             predn_type_1L_chr = NULL,
                             tfmn_for_bnml_1L_lgl = F,
                             family_1L_chr = NA_character_,
                             is_brms_mdl_1L_lgl = T)
{
  predd_utl_dbl <- predict_uncnstrd_utl(data_tb = data_tb,
                                        model_mdl = model_mdl,
                                        new_data_is_1L_chr = new_data_is_1L_chr,
                                        predn_type_1L_chr = predn_type_1L_chr,
                                        tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl,
                                        family_1L_chr = family_1L_chr,
                                        tfmn_1L_chr = tfmn_1L_chr,
                                        is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl)
    # predd_utl_dbl <- predict(model_mdl, newdata = data_tb) %>% calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
    #     tfmn_is_outp_1L_lgl = T) %>% as.vector()
    if(impute_1L_lgl)
      predd_utl_dbl[which(is.na(predd_utl_dbl))] <- predd_utl_dbl %>% na.omit() %>% mean()
    if(force_min_max_1L_lgl){
      predd_utl_dbl[which(predd_utl_dbl>1)] <- 1
      predd_utl_dbl[which(predd_utl_dbl<utl_min_val_1L_dbl)] <- utl_min_val_1L_dbl
    }
    if(!is.null(utl_cls_fn)){
      predd_utl_dbl <- predd_utl_dbl %>% rlang::exec(.fn = utl_cls_fn)
    }
    return(predd_utl_dbl)
}
