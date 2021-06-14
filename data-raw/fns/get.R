get_cndts_for_mxd_mdls <- function(mdl_types_lup = NULL){
  if(is.null(mdl_types_lup))
    utils::data("mdl_types_lup", package = "TTU", envir = environment())
  cndts_for_mxd_mdls_lup <- mdl_types_lup %>%
    dplyr::filter(!tfmn_for_bnml_lgl,
                  short_name_chr != "BET_LOG" )
  return(cndts_for_mxd_mdls_lup)
}
get_link_from_tfmn <- function(tfmn_1L_chr,
                               is_OLS_1L_lgl = F){
  link_1L_chr <- ifelse(is_OLS_1L_lgl,
                        "identity",
                        ifelse(tfmn_1L_chr == "LOG",
                               "log",
                               ifelse(tfmn_1L_chr == "LGT",
                                      "logit",
                                      ifelse(tfmn_1L_chr == "CLL",
                                             "cloglog",
                                             ifelse(tfmn_1L_chr == "LOGLOG",
                                                    "loglog",
                                                    ifelse(tfmn_1L_chr == "NTF",
                                                           "identity",
                                                           "ERROR"))))))
  if(link_1L_chr=="ERROR")
    stop("Link cannot be identified - incorrect transformation argument tfmn_1L_chr")
  return(link_1L_chr)
}
get_signft_covars <- function (mdls_with_covars_smry_tb, covar_var_nms_chr)
{
  signif_vars_chr <- mdls_with_covars_smry_tb$Significant %>%
    purrr::map(~strsplit(.x, " ")) %>% purrr::flatten() %>%
    purrr::flatten_chr() %>% unique()
  signt_covars_chr <- covar_var_nms_chr[covar_var_nms_chr %in%
                                          signif_vars_chr]
  if(identical(signt_covars_chr, character(0)))
    signt_covars_chr <- NA_character_
  return(signt_covars_chr)
}

