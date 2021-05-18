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
