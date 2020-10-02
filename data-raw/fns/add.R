add_aqol6dU_to_aqol6d_items_tb_tb <- function(aqol6d_items_tb,
                                              aqol6d_from_8d_coeffs_lup_tb){
  coeff_dbl <- aqol6d_from_8d_coeffs_lup_tb[match(c(paste0("vD",1:6),"Constant"),
                                                  aqol6d_from_8d_coeffs_lup_tb$var_name_chr),] %>%
    dplyr::pull(coeff_dbl)
  aqol6d_items_tb <- aqol6d_items_tb %>%
    dplyr::mutate(aqol6dU = coeff_dbl[1]*vD1 + coeff_dbl[2]*vD2 + coeff_dbl[3]*vD3 + coeff_dbl[4]*vD4 +coeff_dbl[5]* vD5 +coeff_dbl[6]*vD6 + coeff_dbl[7]) %>%
    dplyr::mutate(aqol6dU = aqol6dU %>% purrr::map_dbl(~ifelse(.x>1,1,.x)))
  return(aqol6d_items_tb)
}
add_aqol6dU_to_tbs_ls <- function(tbs_ls, ## FUNCTION NEEDS TO BE FIXED - UNFINISHED
                                  prefix_1L_chr =  "aqol6d_q",
                                  id_var_nm_1L_chr,
                                  aqol6d_from_8d_coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb,
                                  dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb,
                                  disvalues_lup_tb = disvalues_lup_tb,
                                  itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb){
  tbs_ls <- tbs_ls %>% purrr::map(~.x %>% dplyr::mutate(aqol6dU = calculate_adol_aqol6d(.x, prefix_1L_chr = prefix_1L_chr,id_var_nm_1L_chr = id_var_nm_1L_chr)
                                                          # calculate_aqol6dU_dbl(aqol6d_items_tb = .x,
                                                          #                               prefix_1L_chr =  prefix_1L_chr,
                                                          #                               aqol6d_from_8d_coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb,
                                                          #                               dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb,
                                                          #                               disvalues_lup_tb = disvalues_lup_tb,
                                                          #                               itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb)
                                                          ))
  return(tbs_ls)
}
add_aqol_dim_scrg_eqs <- function(unscored_aqol_tb){
  data("adol_dim_scalg_eqs_lup", package = "FBaqol", envir = environment())
  for(var in adol_dim_scalg_eqs_lup$Dim_scal) {
    expression=adol_dim_scalg_eqs_lup[adol_dim_scalg_eqs_lup$Dim_scal==var,]$Equ
    unscored_aqol_tb <- unscored_aqol_tb %>%
      mutate(!! var := !!rlang::parse_expr(expression))
    Hmisc::label(unscored_aqol_tb[,var])=adol_dim_scalg_eqs_lup[adol_dim_scalg_eqs_lup$Dim_scal==var,]$Label
  }
  return(unscored_aqol_tb)
}
add_aqol_items_tbs_ls <- function(tbs_ls, # Needs to convert study data.
                                  aqol_items_props_tbs_ls,
                                  prefix_chr,
                                  aqol_tots_var_nms_chr,
                                  id_var_nm_1L_chr = "fkClientID",
                                  scaling_cnst_dbl = 5){
  updated_tbs_ls <- purrr::map2(tbs_ls,
                                aqol_items_props_tbs_ls,
                                ~ {
                                  nbr_obs_1L_int <- nrow(.x) * scaling_cnst_dbl
                                  transposed_items_props_tb <- .y %>% dplyr::select(-Question) %>% t()
                                  item_ranges_dbl_ls <- 1:ncol(transposed_items_props_tb) %>% purrr::map(~c(1, length(transposed_items_props_tb[,.x] %>% na.omit())))
                                  cat_probs_def_tbl <- purrr::reduce(1:ncol(transposed_items_props_tb),
                                                                     .init = NULL,
                                                                     ~ simstudy::defData(.x,
                                                                                         varname = paste0("aqol6d_q",.y),
                                                                                         formula = transposed_items_props_tb[,.y] %>% na.omit() %>% as.vector() %>% format(digits=10) %>% paste0(collapse = ";"),
                                                                                         dist = "categorical"))
                                  items_tb <- simstudy::genData(nbr_obs_1L_int, cat_probs_def_tbl) %>%
                                    dplyr::select(-id)  %>%
                                    dplyr::mutate(!!rlang::sym(unname(aqol_tots_var_nms_chr["cumulative"])) := rowSums(.,na.rm = T)) %>%
                                    dplyr::arrange(!!rlang::sym(unname(aqol_tots_var_nms_chr["cumulative"]))) %>%
                                    tibble::rowid_to_column("id")
                                  items_tb <- items_tb  %>%
                                    dplyr::mutate(aqol6dU = calculate_adol_aqol6d(items_tb,
                                                                                  prefix_1L_chr = prefix_chr[2],
                                                                                  id_var_nm_1L_chr = "id"))
                                  .x <- .x %>% dplyr::mutate(id = purrr::map_int(aqol6d_total_w,
                                                                                        ~which.min(abs(items_tb$aqol6dU - .x)))) %>%
                                    dplyr::left_join(items_tb)
                                  updated_tb <- .x %>%
                                    dplyr::mutate(!!rlang::sym(unname(aqol_tots_var_nms_chr["weighted"])) := aqol6dU) %>%
                                    dplyr::select(-aqol6dU,
                                                  -id) %>%
                                    dplyr::select(!!rlang::sym(id_var_nm_1L_chr),
                                                  dplyr::starts_with(unname(prefix_chr[2])),
                                                  !!rlang::sym(unname(aqol_tots_var_nms_chr["cumulative"])),
                                                  !!rlang::sym(unname(aqol_tots_var_nms_chr["weighted"])),
                                                  dplyr::everything())
                                  # %>%
                                  #   dplyr::arrange(aqol6dU)

                                  # .x <- .x %>% dplyr::arrange(aqol6d_total_w)
                                  #%>%
                                    ###dplyr::arrange(rank_x_int) %>%
                                    #dplyr::select(-totals_dbl)  %>% t()
                                  # target_dbl <- .x %>% dplyr::arrange(aqol6d_total_c) %>% dplyr::pull(aqol6d_total_c)
                                  # items_tb <- 1:ncol(items_tb) %>% purrr::map_dfr(~{
                                  #   force_vec_to_sum_to_int(items_tb[,.x],
                                  #                           target_1L_int = target_dbl[.x],
                                  #                           item_ranges_dbl_ls = item_ranges_dbl_ls)
                                  # })
                                  # updated_tb <-dplyr::bind_cols(.x %>% dplyr::arrange(aqol6d_total_c),
                                  #                               items_tb) %>%
                                  #   dplyr::mutate(temp_id = fkClientID %>% purrr::map_dbl(~stringr::str_replace(.x,prefix_chr[["uid"]],"") %>% as.numeric())) %>%
                                  #   dplyr::arrange(temp_id) %>%
                                  #   dplyr::select(-temp_id) %>%
                                  #   dplyr::select(fkClientID, dplyr::everything())
                                  updated_tb
                                })
  return(updated_tbs_ls)
}
add_aqol_scores_tbs_ls <- function(tbs_ls,
                                   means_dbl,
                                   sds_dbl,
                                   corr_dbl){
  tbs_ls <- purrr::pmap(list(tbs_ls,
                             means_dbl,
                             sds_dbl),
                        ~ {
                          aqol_score_dbl <- faux::rnorm_pre(..1 %>% dplyr::pull(aqol6d_total_w),
                                                            mu = ..2,
                                                            sd = ..3,
                                                            r = corr_dbl)
                          aqol_score_dbl <- aqol_score_dbl %>% purrr::map_dbl(~min(round(.x),99) %>% max(20))
                          ..1 %>% dplyr::mutate(aqol6d_total_c = tidyselect::all_of(aqol_score_dbl))
                        }
  )
  return(tbs_ls)
}
add_corrs_and_uts_to_tbs_ls_ls <- function(tbs_ls, # Based on: https://stats.stackexchange.com/questions/134164/how-to-rearrange-2d-data-to-get-given-correlation
                                           aqol_scores_pars_ls,
                                           aqol_items_props_tbs_ls,
                                           temporal_corrs_ls,
                                           prefix_chr,
                                           aqol_tots_var_nms_chr,
                                           id_var_nm_1L_chr = "fkClientID"){
  tbs_ls <- reorder_tbs_for_target_cors(tbs_ls,
                                       corr_dbl = temporal_corrs_ls[[1]],
                                       corr_var_chr = rep(names(temporal_corrs_ls)[1],2),
                                       id_var_to_rm_1L_chr = "id"
                                       ) %>%
    add_uids_to_tbs_ls(prefix_1L_chr = prefix_chr[["uid"]],
                       id_var_nm_1L_chr = id_var_nm_1L_chr)
  tbs_ls <- tbs_ls  %>%
    add_aqol_items_tbs_ls(aqol_items_props_tbs_ls = aqol_items_props_tbs_ls,
                          prefix_chr = prefix_chr,
                          aqol_tots_var_nms_chr = aqol_tots_var_nms_chr,
                          id_var_nm_1L_chr = id_var_nm_1L_chr) #%>%
    # add_aqol6dU_to_tbs_ls(prefix_1L_chr =  prefix_chr[["aqol_item"]],
    #                       aqol6d_from_8d_coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb,
    #                       id_var_nm_1L_chr = "fkClientID",
    #                       dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb,
    #                       disvalues_lup_tb = disvalues_lup_tb,
    #                       itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb)
  return(tbs_ls)
}
add_dmn_disu_to_aqol6d_items_tb_tb <- function(aqol6d_items_tb,
                                               domain_items_ls,
                                               domains_chr,
                                               dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb,
                                               itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb){
  aqol6d_disu_fn_ls <- make_aqol6d_fns_ls(domain_items_ls)
  kD_dbl <- make_dim_sclg_cons_dbl(domains_chr = domains_chr,
                                   dim_sclg_constant_lup_tb = dim_sclg_constant_lup_tb)
  w_dbl_ls <- make_item_wrst_wghts_ls_ls(domain_items_ls = domain_items_ls,
                                         itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb)
  aqol6d_items_tb <- purrr::reduce(1:length(domain_items_ls),
                                   .init = aqol6d_items_tb,
                                   ~{
                                     args_ls <- list(dvQs_tb =  .x %>%
                                                       dplyr::select(domain_items_ls[[.y]] %>%
                                                                       paste0("dv_",.)),
                                                     kD_1L_dbl = kD_dbl[.y],
                                                     w_dbl = w_dbl_ls[[.y]])
                                     .x %>% dplyr::mutate(!!rlang::sym(paste0("dvD",.y)) := rlang::exec(aqol6d_disu_fn_ls[[.y]],
                                                                                                        !!!args_ls))

                                   })
  return(aqol6d_items_tb)
}
add_dmn_scores_to_aqol6d_items_tb_tb <- function(aqol6d_items_tb,
                                                 domain_items_ls){
  aqol6d_items_tb <- aqol6d_items_tb %>% dplyr::mutate(dplyr::across(paste0("dvD",1:length(domain_items_ls)),
                                                                     .fns = list(vD = ~ 1-.x
                                                                     ),
                                                                     .names = "{fn}_{col}"
  )) %>% dplyr::rename_with(~stringr::str_replace(.,"vD_dvD", "vD"))
  return(aqol6d_items_tb)
}
add_domain_unwtd_tots_tb <- function(items_tb,
                                     domain_items_ls,
                                     domain_pfx_1L_chr){
  items_and_domains_tb <- purrr::reduce(1:length(domain_items_ls), .init = items_tb,
                                        ~ .x %>% dplyr::mutate(!!rlang::sym(paste0(domain_pfx_1L_chr,names(domain_items_ls)[.y])) := rowSums(dplyr::select(.,domain_items_ls[[.y]]))))
  return(items_and_domains_tb)
}
add_itm_disu_to_aqol6d_itms_tb_tb <- function(aqol6d_items_tb,
                                              disvalues_lup_tb = disvalues_lup_tb,
                                              pfx_1L_chr){
  aqol6d_items_tb <- purrr::reduce(1:20,
                                   .init = aqol6d_items_tb,
                                   ~ {
                                     q_1L_chr <- paste0(pfx_1L_chr,.y)
                                     disu_dbl <- disvalues_lup_tb[.y,-1] %>% as.numeric()
                                     .x %>%
                                       dplyr::mutate(dplyr::across(tidyselect::all_of(q_1L_chr),
                                                                   .fns = list(dv = ~ disu_dbl[.x]
                                                                   ),
                                                                   .names = "{fn}_{col}" ) )
                                   })
  return(aqol6d_items_tb)
}
add_labels_to_aqol6d_tb <- function(aqol6d_tb,
                                    labels_chr = NA_character_){
  if(is.na(labels_chr))
    labels_chr <- c(
      fkClientID ="Unique client identifier",
      round = "Data measurement round",
      d_age = "Age",
      d_gender = "Gender",
      d_sexual_ori_s = "Sexual orientation",
      d_studying_working = "Work and study",
      c_p_diag_s = " Primary diagnosis",
      c_clinical_staging_s = "Clinical stage",
      c_sofas = "SOFAS",
      s_centre = "Clinic",
      d_agegroup = "Age group",
      d_sex_birth_s = "Sex at birth",
      d_country_bir_s = "Country of birth",
      d_ATSI = "Aboriginal and Torres Strait Islander",
      d_english_home = "English spoken at home",
      d_english_native = "English is native language",
      d_relation_s = "Relationship status",
      aqol6d_total_w ="AQoL health utility",
      phq9_total = "PHQ9",
      bads_total = "BADS",
      gad7_total = "GAD7",
      oasis_total = "OASIS",
      scared_total = "SCARED",
      k6_total = "K6",
      aqol6d_total_c = "AQoL unweighted total",
      aqol6d_q1 ="Household tasks",
      aqol6d_q2 = "Getting around",
      aqol6d_q3 = "Mobility",
      aqol6d_q4 = "Self care",
      aqol6d_q5 = "Enjoy close rels", # omitted regex
      aqol6d_q6 = "Family rels", # omitted regex
      aqol6d_q7 = "Community involvement",
      aqol6d_q8 = "Despair",
      aqol6d_q9 = "Worry",
      aqol6d_q10 = "Sad",
      aqol6d_q11 = "Agitated",
      aqol6d_q12 = "Energy level",
      aqol6d_q13 = "Control",
      aqol6d_q14 = "Coping",
      aqol6d_q15 = "Frequency of pain",
      aqol6d_q16 = "Degree of pain",
      aqol6d_q17 = "Pain interference",
      aqol6d_q18 = "Vision",
      aqol6d_q19 = "Hearing",
      aqol6d_q20 = "Communication",
      aqol6d_subtotal_w_IL = "Independent Living",
      aqol6d_subtotal_w_Rel = "Relationships",
      aqol6d_subtotal_w_MH = "Mental Health",
      aqol6d_subtotal_w_COP = "Coping",
      aqol6d_subtotal_w_P = "Pain",
      aqol6d_subtotal_w_SEN = "Sense")
  Hmisc::label(aqol6d_tb) = as.list(labels_chr[match(names(aqol6d_tb), names(labels_chr))])
  return(aqol6d_tb)
}
add_uids_to_tbs_ls <- function(tbs_ls,
                               prefix_1L_chr,
                               id_var_nm_1L_chr = "fkClientID"){
  participant_ids <- paste0(prefix_1L_chr,1:nrow(tbs_ls$bl_part_1_tb)) %>% sample(nrow(tbs_ls$bl_part_1_tb))
  tbs_ls <- purrr::map(tbs_ls,
                       ~ {
                         .x %>%
                           dplyr::mutate(!!rlang::sym(id_var_nm_1L_chr) := tidyselect::all_of(participant_ids[1:nrow(.x)])) %>%
                           dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr) %>% purrr::map_chr(~stringr::str_replace(.x,prefix_1L_chr,"")) %>% as.numeric())
                       }) %>%
    stats::setNames(names(tbs_ls))
  return(tbs_ls)
}










