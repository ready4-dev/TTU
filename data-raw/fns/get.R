get_cndts_for_mxd_mdls <- function(mdl_types_lup = NULL){
  if(is.null(mdl_types_lup))
    utils::data("mdl_types_lup", package = "TTU", envir = environment())
  cndts_for_mxd_mdls_lup <- mdl_types_lup %>%
    dplyr::filter(!tfmn_for_bnml_lgl,
                  short_name_chr != "BET_LOG" )
  return(cndts_for_mxd_mdls_lup)
}
get_mdl_type_from_nm <- function(mdl_nm_1L_chr,
                                 mdl_types_lup = NULL){
  if(is.null(mdl_types_lup))
    utils::data("mdl_types_lup", package = "TTU", envir = environment())
  mdl_type_1L_chr <- (mdl_types_lup %>%
                        dplyr::pull(short_name_chr))[mdl_types_lup %>%
                                                       dplyr::pull(short_name_chr) %>%
                                                       purrr::map_lgl(~endsWith(mdl_nm_1L_chr,.x))]
  return(mdl_type_1L_chr)
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
get_random_intercept <- function(mdls_smry_tb,
                                 mdl_nm_1L_chr,
                                 deterministic_1L_lgl = T){
  mdl_smry_tb <- mdls_smry_tb %>%
    dplyr::filter(Model == mdl_nm_1L_chr)
  sd_dbl <- c(mdl_smry_tb %>%
                ready4fun::get_from_lup_obj(match_value_xx = "SD (Intercept)",
                                            match_var_nm_1L_chr = "Parameter",
                                            target_var_nm_1L_chr = "Estimate",
                                            evaluate_lgl = F),
              ifelse(deterministic_1L_lgl,
                     0,
                     mdl_smry_tb %>%
                       ready4fun::get_from_lup_obj(match_value_xx = "SD (Intercept)",
                                                   match_var_nm_1L_chr = "Parameter",
                                                   target_var_nm_1L_chr = "SE",
                                                   evaluate_lgl = F)))
  return(sd_dbl)
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
get_table_predn_mdl <- function(mdl_nm_1L_chr,
                                ingredients_ls){
  mdl_type_1L_chr <- get_mdl_type_from_nm(selected_mdl,
                                          mdl_types_lup = ingredients_ls$mdl_types_lup)
  table_predn_mdl <- make_shareable_mdl(fake_ds_tb = ingredients_ls$fake_ds_tb,
                                        mdl_smry_tb = mdls_smry_ls %>% purrr::pluck(selected_mdl),
                                        depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr,
                                        id_var_nm_1L_chr = ingredients_ls$id_var_nm_1L_chr,
                                        tfmn_1L_chr = ready4fun::get_from_lup_obj(ingredients_ls$mdl_types_lup,
                                                                                  match_value_xx = mdl_type_1L_chr,
                                                                                  match_var_nm_1L_chr = "short_name_chr",
                                                                                  target_var_nm_1L_chr = "tfmn_chr",
                                                                                  evaluate_lgl = F),
                                        mdl_type_1L_chr = mdl_type_1L_chr,
                                        mdl_types_lup = ingredients_ls$mdl_types_lup,
                                        control_1L_chr = ready4fun::get_from_lup_obj(ingredients_ls$mdl_types_lup,
                                                                                     match_value_xx = mdl_type_1L_chr,
                                                                                     match_var_nm_1L_chr = "short_name_chr",
                                                                                     target_var_nm_1L_chr = "control_chr",
                                                                                     evaluate_lgl = F),
                                        start_1L_chr = NA_character_,
                                        seed_1L_int = ingredients_ls$seed_1L_int)
  return(table_predn_mdl)
}
