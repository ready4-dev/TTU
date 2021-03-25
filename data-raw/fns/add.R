
add_cors_and_uts_to_aqol6d_tbs_ls <- function (aqol6d_tbs_ls, aqol_scores_pars_ls, aqol_items_props_tbs_ls,
                                               temporal_cors_ls, prefix_chr, aqol_tots_var_nms_chr, id_var_nm_1L_chr = "fkClientID")
{
  aqol6d_tbs_ls <- reorder_tbs_for_target_cors(aqol6d_tbs_ls,
                                               cor_dbl = temporal_cors_ls[[1]], cor_var_chr = rep(names(temporal_cors_ls)[1],
                                                                                                  2), id_var_to_rm_1L_chr = "id") %>% add_uids_to_tbs_ls(prefix_1L_chr = prefix_chr[["uid"]],
                                                                                                                                                         id_var_nm_1L_chr = id_var_nm_1L_chr)
  aqol6d_tbs_ls <- aqol6d_tbs_ls %>% add_aqol6d_items_to_aqol6d_tbs_ls(aqol_items_props_tbs_ls = aqol_items_props_tbs_ls,
                                                                       prefix_chr = prefix_chr, aqol_tots_var_nms_chr = aqol_tots_var_nms_chr,
                                                                       id_var_nm_1L_chr = id_var_nm_1L_chr)
  return(aqol6d_tbs_ls)
}
add_interval_var <- function(data_tb,
                             id_var_nm_1L_chr = "fkClientID",
                             msrmnt_date_var_nm_1L_chr = "d_interview_date",
                             time_unit_1L_chr = "days",
                             bl_date_var_nm_1L_chr = "bl_date_dtm",
                             interval_var_nm_1L_chr = "interval_dbl",
                             temp_row_nbr_var_nm_1L_chr = "temp_row_nbr_int",
                             drop_bl_date_var = F){
  updated_data_tb <- data_tb %>%
    dplyr::ungroup() %>%
    dplyr::mutate(!!rlang::sym(temp_row_nbr_var_nm_1L_chr) := 1:dplyr::n()) %>%
    dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>%
    dplyr::arrange(!!rlang::sym(msrmnt_date_var_nm_1L_chr)) %>%
    dplyr::mutate(!!rlang::sym(bl_date_var_nm_1L_chr) := !!rlang::sym(msrmnt_date_var_nm_1L_chr) %>% dplyr::first()) %>%
    dplyr::mutate(interval_dbl = purrr::map2_dbl(!!rlang::sym(bl_date_var_nm_1L_chr),
                                                 !!rlang::sym(msrmnt_date_var_nm_1L_chr),
                                                 ~ lubridate::interval(.x, .y) %>%
                                                   lubridate::time_length(unit = time_unit_1L_chr))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(!!rlang::sym(temp_row_nbr_var_nm_1L_chr)) %>%
    dplyr::select(-!!rlang::sym(temp_row_nbr_var_nm_1L_chr))
  if(drop_bl_date_var)
    updated_data_tb <-  updated_data_tb %>%
      dplyr::select(-!!rlang::sym(bl_date_var_nm_1L_chr))
  return(updated_data_tb)
}
add_labels_to_aqol6d_tb <- function (aqol6d_tb, labels_chr = NA_character_)
{
  if (is.na(labels_chr))
    labels_chr <- c(fkClientID = "Unique client identifier",
                    round = "Data measurement round", d_age = "Age",
                    d_gender = "Gender", d_sexual_ori_s = "Sexual orientation",
                    d_studying_working = "Work and study", c_p_diag_s = " Primary diagnosis",
                    c_clinical_staging_s = "Clinical stage", c_sofas = "SOFAS",
                    s_centre = "Clinic", d_agegroup = "Age group", d_sex_birth_s = "Sex at birth",
                    d_country_bir_s = "Country of birth", d_ATSI = "Aboriginal and Torres Strait Islander",
                    d_english_home = "English spoken at home", d_english_native = "English is native language",
                    d_relation_s = "Relationship status", aqol6d_total_w = "AQoL health utility",
                    phq9_total = "PHQ9", bads_total = "BADS", gad7_total = "GAD7",
                    oasis_total = "OASIS", scared_total = "SCARED", k6_total = "K6",
                    aqol6d_total_c = "AQoL unweighted total", aqol6d_q1 = "Household tasks",
                    aqol6d_q2 = "Getting around", aqol6d_q3 = "Mobility",
                    aqol6d_q4 = "Self care", aqol6d_q5 = "Enjoy close rels",
                    aqol6d_q6 = "Family rels", aqol6d_q7 = "Community involvement",
                    aqol6d_q8 = "Despair", aqol6d_q9 = "Worry", aqol6d_q10 = "Sad",
                    aqol6d_q11 = "Agitated", aqol6d_q12 = "Energy level",
                    aqol6d_q13 = "Control", aqol6d_q14 = "Coping", aqol6d_q15 = "Frequency of pain",
                    aqol6d_q16 = "Degree of pain", aqol6d_q17 = "Pain interference",
                    aqol6d_q18 = "Vision", aqol6d_q19 = "Hearing", aqol6d_q20 = "Communication",
                    aqol6d_subtotal_c_IL = "Unweighted Independent Living",
                    aqol6d_subtotal_c_REL = "Unweighted Relationships",
                    aqol6d_subtotal_c_MH = "Unweighted Mental Health",
                    aqol6d_subtotal_c_COP = "Unweighted Coping", aqol6d_subtotal_c_P = "Unweighted Pain",
                    aqol6d_subtotal_c_SEN = "Unweighted Sense", aqol6d_subtotal_w_IL = "Independent Living",
                    aqol6d_subtotal_w_REL = "Relationships", aqol6d_subtotal_w_MH = "Mental Health",
                    aqol6d_subtotal_w_COP = "Coping", aqol6d_subtotal_w_P = "Pain",
                    aqol6d_subtotal_w_SEN = "Sense")
  Hmisc::label(aqol6d_tb) = as.list(labels_chr[match(names(aqol6d_tb),
                                                     names(labels_chr))])
  return(aqol6d_tb)
}
add_uids_to_tbs_ls <- function (tbs_ls, prefix_1L_chr, id_var_nm_1L_chr = "fkClientID")
{
  participant_ids <- paste0(prefix_1L_chr, 1:nrow(tbs_ls[[1]])) %>%
    sample(nrow(tbs_ls[[1]]))
  tbs_ls <- purrr::map(tbs_ls, ~{
    .x %>% dplyr::mutate(`:=`(!!rlang::sym(id_var_nm_1L_chr),
                              tidyselect::all_of(participant_ids[1:nrow(.x)]))) %>%
      dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr) %>%
                       purrr::map_chr(~stringr::str_replace(.x, prefix_1L_chr,
                                                            "")) %>% as.numeric())
  }) %>% stats::setNames(names(tbs_ls))
  return(tbs_ls)
}

add_utility_predn_to_ds <- function (data_tb, model_mdl, tfmn_1L_chr, dep_var_nm_1L_chr, predr_vars_nms_chr = NULL,
                                     force_min_max_1L_lgl = T, utl_min_val_1L_dbl = 0.03, impute_1L_lgl = T, utl_cls_fn = NULL,
                                     rmv_tfmd_dep_var_1L_lgl = F)
{
    dep_vars_chr <- c(dep_var_nm_1L_chr, transform_dep_var_nm(dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        tfmn_1L_chr = tfmn_1L_chr)) %>% unique()
    data_tb <- purrr::reduce(dep_vars_chr, .init = data_tb, ~dplyr::mutate(.x,
                                                                           !!rlang::sym(.y):= NA_real_))
    predictions_dbl <- predict_utility(data_tb = data_tb, tfmn_1L_chr = tfmn_1L_chr,
                                       model_mdl = model_mdl,
                                       force_min_max_1L_lgl = force_min_max_1L_lgl,
                                       utl_min_val_1L_dbl = utl_min_val_1L_dbl,
                                       impute_1L_lgl = impute_1L_lgl,
                                       utl_cls_fn = utl_cls_fn)
    data_tb <- data_tb %>% dplyr::mutate(!!rlang::sym(dep_var_nm_1L_chr):=predictions_dbl)
    if(!is.null(predr_vars_nms_chr)){
      data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(purrr::map(predr_vars_nms_chr, ~ paste0(.x,c("_baseline","_change"))) %>%
                                             purrr::flatten_chr()))
    }
    if(rmv_tfmd_dep_var_1L_lgl){
      data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(dep_vars_chr[dep_vars_chr!=dep_var_nm_1L_chr]))
    }
    return(data_tb)
}

