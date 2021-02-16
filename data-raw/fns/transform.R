transform_data_tb_for_cmprsn <- function (data_tb, model_mdl, dep_var_nm_1L_chr = "utl_total_w",
    source_data_nm_1L_chr = "Original", tf_type_1L_chr = "Predicted",
    pred_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_)
{
    if (tf_type_1L_chr == "Predicted")
        new_data_dbl <- stats::predict(model_mdl, type = pred_type_1L_chr)
    if (tf_type_1L_chr == "Simulated" & !tfmn_for_bnml_1L_lgl)
        new_data_dbl <- stats::simulate(model_mdl)$sim_1
    if (tf_type_1L_chr == "Simulated" & tfmn_for_bnml_1L_lgl)
        new_data_dbl <- (stats::predict(model_mdl) + stats::rnorm(nrow(data_tb),
            0, stats::sigma(model_mdl))) %>% calculate_dep_var_tfmn(tfmn_1L_chr = ifelse(family_1L_chr ==
            "quasibinomial(log)", "LOG", ifelse(family_1L_chr ==
            "quasibinomial(logit)", "LOGIT", ifelse(family_1L_chr ==
            "quasibinomial(cloglog)", "CLL", "NTF"))), tfmn_is_outp_1L_lgl = T)
    tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(tf_type_1L_chr),
        new_data_dbl), `:=`(!!rlang::sym(source_data_nm_1L_chr),
        !!rlang::sym(dep_var_nm_1L_chr)))
    return(tfd_data_tb)
}
transform_dep_var_nm <- function (dep_var_nm_1L_chr, tfmn_1L_chr = "NTF")
{
    tfd_dep_var_nm_1L_chr <- paste0(dep_var_nm_1L_chr, ifelse(tfmn_1L_chr ==
        "NTF", "", paste0("_", tfmn_1L_chr)))
    return(tfd_dep_var_nm_1L_chr)
}
transform_dep_var_nm_for_cll <- function (dep_var_nm_1L_chr)
{
    tfd_dep_var_nm_1L_chr <- paste0(dep_var_nm_1L_chr, "_cloglog")
    return(tfd_dep_var_nm_1L_chr)
}
transform_ds_for_tstng <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", dep_var_max_val_1L_dbl = 0.999,
                                    candidate_predrs_chr = NA_character_, covar_var_nms_chr = NA_character_,
                                    round_var_nm_1L_chr = "round", round_val_1L_chr = "Baseline",
                                    remove_all_mssng_1L_lgl = F)
{
  vars_to_keep_chr <- c(dep_var_nm_1L_chr, candidate_predrs_chr,
                        covar_var_nms_chr) %>% purrr::discard(is.na)
  tfd_data_tb <- data_tb %>% dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) ==
                                             round_val_1L_chr) %>% dplyr::select(!!!rlang::syms(vars_to_keep_chr)) %>%
    dplyr::mutate(`:=`(!!rlang::sym(dep_var_nm_1L_chr), ifelse(!!rlang::sym(dep_var_nm_1L_chr) >
                                                                 dep_var_max_val_1L_dbl, dep_var_max_val_1L_dbl, !!rlang::sym(dep_var_nm_1L_chr))))
  if (remove_all_mssng_1L_lgl)
    tfd_data_tb <- tfd_data_tb %>% stats::na.omit()
  return(tfd_data_tb)
}
transform_ds_for_mdlng <- function (data_tb, dep_var_nm_1L_chr = "utl_total_w", predr_var_nm_1L_chr,
    covar_var_nms_chr = NA_character_)
{
    mdl_vars_chr <- c(names(data_tb)[names(data_tb) %>% startsWith(dep_var_nm_1L_chr)],
        predr_var_nm_1L_chr, covar_var_nms_chr) %>% purrr::discard(is.na)
    tfd_data_tb <- data_tb %>% tidyr::drop_na(!!!rlang::syms(mdl_vars_chr)) %>%
        dplyr::select(!!!rlang::syms(mdl_vars_chr))
    return(tfd_data_tb)
}
transform_raw_aqol_tb_to_aqol6d_tb <- function (raw_aqol_tb)
{
  aqol6d_tb <- raw_aqol_tb %>% dplyr::mutate(d_agegroup = cut(d_age,
                                                              breaks = c(11, 17, 30), labels = c("Age 12-17", "Age 18-26"))) %>%
    dplyr::mutate(round = factor(round, labels = c("Baseline",
                                                   "Follow-up"))) %>% dplyr::select(fkClientID, c_p_diag_s,
                                                                                    s_centre, c_clinical_staging_s, d_age, d_agegroup, d_gender,
                                                                                    d_sex_birth_s, d_sexual_ori_s, d_country_bir_s, d_ATSI,
                                                                                    d_english_home, d_english_native, d_relation_s, d_studying_working,
                                                                                    k6_total, phq9_total, bads_total, gad7_total, oasis_total,
                                                                                    scared_total, dplyr::contains("aqol6d"), c_sofas, round) %>%
    dplyr::mutate(Gender = factor(ifelse(d_gender == "Genderqueer/gender nonconforming/agender" |
                                           d_gender == "Transgender", "Other", as.character(d_gender)))) %>%
    dplyr::mutate(Region = as.factor(ifelse(s_centre == "Canberra" |
                                              s_centre == "Southport" | s_centre == "Knox", "Metro",
                                            "Regional"))) %>% dplyr::mutate(CALD = factor(ifelse(d_country_bir_s ==
                                                                                                   "Other" | d_english_home == "No" | d_english_native ==
                                                                                                   "No" | d_ATSI == "Yes", "Yes", "No"))) %>% dplyr::rename(PHQ9 = phq9_total,
                                                                                                                                                            BADS = bads_total, GAD7 = gad7_total, OASIS = oasis_total,
                                                                                                                                                            SCARED = scared_total, K6 = k6_total, SOFAS = c_sofas)
  # Hmisc::label(aqol6d_tb$CALD) = "Culturally and linguistically diverse (CALD) background"
  # Hmisc::label(aqol6d_tb$d_agegroup) = "Age group"
  return(aqol6d_tb)
}
transform_tb_to_mdl_inp <- function (data_tb, dep_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr,
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round",
    round_bl_val_1L_chr = "Baseline", drop_all_msng_1L_lgl = T, scaling_fctr_dbl = 0.01, ungroup_1L_lgl = F,
    add_cll_tfmn_1L_lgl = T)
{
    if(length(scaling_fctr_dbl)!=length(predr_vars_nms_chr)){
      scaling_fctr_dbl <- rep(scaling_fctr_dbl[1],length(predr_vars_nms_chr))
    }
  data_tb <- data.frame(data_tb)
    tfd_for_mdl_inp_tb <- data_tb %>% dplyr::select(dplyr::all_of(id_var_nm_1L_chr), dplyr::all_of(round_var_nm_1L_chr),
        dplyr::all_of(predr_vars_nms_chr),
        dplyr::all_of(dep_var_nm_1L_chr) # Moved from first var in tb
        ) %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>%
        dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr), !!rlang::sym(round_var_nm_1L_chr))
    tfd_for_mdl_inp_tb <- purrr::reduce(1:length(predr_vars_nms_chr),
                                        .init = tfd_for_mdl_inp_tb,
                                        ~ {
                                          idx_1L_int <- as.integer(.y)
                                          .x %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr[idx_1L_int]),
                                                                             .fns = list(baseline = ~dplyr::first(.)*scaling_fctr_dbl[idx_1L_int],
                                                                                         change = ~ifelse(!!rlang::sym(round_var_nm_1L_chr) == "Baseline",
                                                                                                          0,
                                                                                                          (. - dplyr::lag(.))*scaling_fctr_dbl[idx_1L_int]))))
                                          })
    if(add_cll_tfmn_1L_lgl){
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
      dplyr::mutate(`:=`(!!rlang::sym(transform_dep_var_nm_for_cll(dep_var_nm_1L_chr)),
        log(-log(1 - ifelse(!!rlang::sym(dep_var_nm_1L_chr) ==
            1, 0.999, !!rlang::sym(dep_var_nm_1L_chr))))))
    }
    if(drop_all_msng_1L_lgl){
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
        stats::na.omit()
    }
    if(ungroup_1L_lgl){
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
        dplyr::ungroup()
    }
    return(tfd_for_mdl_inp_tb)
}
transform_ts_mdl_data <- function (mdl_ls, data_tb, dep_var_nm_1L_chr = "aqol6d_total_w",
                                   predr_vars_nms_chr, id_var_nm_1L_chr = "fkClientID", mdl_nm_1L_chr)
{
  old_data_tb <- data_tb %>% dplyr::select(c(dplyr::all_of(id_var_nm_1L_chr),
                                             dplyr::all_of(dep_var_nm_1L_chr), predr_vars_nms_chr %>%
                                               purrr::map(~paste0(.x, c("", "_baseline", "_change"))) %>%
                                               purrr::flatten_chr()))
  cnfdl_mdl_ls <- mdl_ls
  cnfdl_mdl_ls$data <- old_data_tb %>% as.data.frame() %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~sample(.x,
                                                                1)))
  return(cnfdl_mdl_ls)
}
