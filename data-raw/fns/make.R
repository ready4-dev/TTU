make_adol_aqol6d_disv_lup <- function ()
{
  adol_aqol6d_disv_lup <- aqol6d_adult_disv_lup_tb %>% dplyr::mutate(Answer_4_dbl = dplyr::case_when(Question_chr ==
                                                                                                       "Q18" ~ 0.622, TRUE ~ Answer_4_dbl), Answer_5_dbl = dplyr::case_when(Question_chr ==
                                                                                                                                                                              "Q3" ~ 0.827, TRUE ~ Answer_5_dbl), Answer_6_dbl = dplyr::case_when(Question_chr ==
                                                                                                                                                                                                                                                    "Q1" ~ 0.073, TRUE ~ Answer_5_dbl))
  return(adol_aqol6d_disv_lup)
}
make_aqol6d_adol_pop_tbs_ls <- function (aqol_items_props_tbs_ls, aqol_scores_pars_ls, series_names_chr,
                                         synth_data_spine_ls, temporal_cors_ls, id_var_nm_1L_chr = "fkClientID",
                                         prefix_chr = c(uid = "Participant_", aqol_item = "aqol6d_q",
                                                        domain_unwtd_pfx_1L_chr = "aqol6d_subtotal_c_", domain_wtd_pfx_1L_chr = "aqol6d_subtotal_w_"))
{
  item_pfx_1L_chr <- prefix_chr[["aqol_item"]]
  uid_pfx_1L_chr <- prefix_chr[["uid"]]
  aqol6d_adol_pop_tbs_ls <- make_synth_series_tbs_ls(synth_data_spine_ls,
                                                     series_names_chr = series_names_chr) %>% add_cors_and_uts_to_aqol6d_tbs_ls(aqol_scores_pars_ls = aqol_scores_pars_ls,
                                                                                                                                aqol_items_props_tbs_ls = aqol_items_props_tbs_ls, temporal_cors_ls = temporal_cors_ls,
                                                                                                                                prefix_chr = prefix_chr, aqol_tots_var_nms_chr = synth_data_spine_ls$aqol_tots_var_nms_chr,
                                                                                                                                id_var_nm_1L_chr = id_var_nm_1L_chr) %>% purrr::map(~{
                                                                                                                                  domain_items_ls <- make_domain_items_ls(domain_qs_lup_tb = aqol6d_domain_qs_lup_tb,
                                                                                                                                                                          item_pfx_1L_chr = item_pfx_1L_chr)
                                                                                                                                  domain_items_ls %>% add_unwtd_dim_tots(items_tb = .x,
                                                                                                                                                                         domain_pfx_1L_chr = prefix_chr[["domain_unwtd_pfx_1L_chr"]]) %>%
                                                                                                                                    add_wtd_dim_tots(domain_items_ls = domain_items_ls,
                                                                                                                                                     domain_unwtd_pfx_1L_chr = prefix_chr[["domain_unwtd_pfx_1L_chr"]],
                                                                                                                                                     domain_wtd_pfx_1L_chr = prefix_chr[["domain_wtd_pfx_1L_chr"]]) %>%
                                                                                                                                    add_labels_to_aqol6d_tb()
                                                                                                                                }) %>% purrr::map(~.x %>% dplyr::select(!!rlang::sym(id_var_nm_1L_chr),
                                                                                                                                                                        dplyr::starts_with(item_pfx_1L_chr), dplyr::starts_with(prefix_chr[["domain_unwtd_pfx_1L_chr"]]),
                                                                                                                                                                        dplyr::starts_with(prefix_chr[["domain_wtd_pfx_1L_chr"]]),
                                                                                                                                                                        dplyr::everything()))
  return(aqol6d_adol_pop_tbs_ls)
}
make_aqol6d_fns_ls <- function (domain_items_ls)
{
  aqol6d_disu_fn_ls <- paste0("calculate_aqol6d_dim_", 1:length(domain_items_ls),
                              "_disv") %>% purrr::map(~rlang::sym(.x))
  return(aqol6d_disu_fn_ls)
}
make_aqol6d_items_tb <- function (aqol_tb, old_pfx_1L_chr, new_pfx_1L_chr)
{
  aqol6d_items_tb <- aqol_tb %>% dplyr::select(dplyr::starts_with(old_pfx_1L_chr)) %>%
    dplyr::rename_all(~{
      stringr::str_replace(., old_pfx_1L_chr, new_pfx_1L_chr)
    })
  return(aqol6d_items_tb)
}
make_brms_mdl_print_ls <- function (mdl_ls, label_stub_1L_chr, caption_1L_chr, output_type_1L_chr = "PDF",
                                    digits_1L_dbl = 2, big_mark_1L_chr = " ")
{
    smry_mdl_ls <- summary(mdl_ls, digits = 4)
    mdl_smry_chr <- smry_mdl_ls %>% utils::capture.output()
    idx_dbl <- c("Formula: ", "Samples: ", "Group-Level Effects: ",
                 "Population-Level Effects: ", "Family Specific Parameters: ",
                 "Samples were drawn using ") %>% purrr::map_dbl(~mdl_smry_chr %>%
                                                                     startsWith(.x) %>% which())
    data_tb <- make_brms_mdl_smry_tbl(smry_mdl_ls, grp_1L_chr = mdl_smry_chr[idx_dbl[3]],
                                      pop_1L_chr = mdl_smry_chr[idx_dbl[4]], fam_1L_chr = mdl_smry_chr[idx_dbl[5]])
    bold_lgl <- data_tb$Parameter %in% c(mdl_smry_chr[idx_dbl[3]],
                                         mdl_smry_chr[idx_dbl[4]], mdl_smry_chr[idx_dbl[5]])
    if (output_type_1L_chr == "PDF") {
        data_tb <- data_tb %>% dplyr::mutate(Parameter = purrr::map_chr(Parameter,
                                                                        ~.x %>% Hmisc::latexTranslate()))
    }
    data_tb <- data_tb %>% dplyr::mutate(Parameter = Parameter %>%
                                             purrr::map2_chr(dplyr::all_of(bold_lgl), ~ifelse(.y &
                                                                                                  output_type_1L_chr == "PDF", paste0("\\textbf{",
                                                                                                                                      .x, "}"), .x)))
    if (output_type_1L_chr != "PDF") {
        data_tb <- data_tb %>% dplyr::mutate(dplyr::across(c(Bulk_ESS,
                                                             Tail_ESS), ~format(., big.mark = big_mark_1L_chr)))
    }
    if (output_type_1L_chr == "HTML") {
        data_tb <- data_tb %>% dplyr::mutate(dplyr::across(where(is.numeric),
                                                           ~format(round(., digits = digits_1L_dbl), digits = digits_1L_dbl,
                                                                   nsmall = digits_1L_dbl)))
    }
    data_tb <- data_tb %>% dplyr::mutate(dplyr::across(where(is.character),
                                                       ~dplyr::case_when(is.na(.) ~ "", . == "NA" ~ "", endsWith(.,
                                                                                                                 " NA") ~ "", TRUE ~ .)))
    end_matter_1L_chr <- trimws(mdl_smry_chr[idx_dbl[6]:length(mdl_smry_chr)]) %>%
        paste0(collapse = " ")
    brms_mdl_print_ls <- list(part_1 = mdl_smry_chr[idx_dbl[1]],
                              part_2 = "\n\n", part_3 = c(trimws(mdl_smry_chr[1:(idx_dbl[2] -
                                                                                     1)][-idx_dbl[1]]), paste0(trimws(mdl_smry_chr[idx_dbl[2]]),
                                                                                                               " ", trimws(mdl_smry_chr[idx_dbl[2] + 1]), collapse = " ")) %>%
                                  paste0(collapse = ifelse(output_type_1L_chr == "PDF",
                                                           "\n\n", "\n")), part_4 = "\n\n", part_5 = list(data_tb = data_tb,
                                                                                                          output_type_1L_chr = output_type_1L_chr, caption_1L_chr = caption_1L_chr,
                                                                                                          mkdn_tbl_ref_1L_chr = paste0("tab:", label_stub_1L_chr), merge_row_idx_int = as.integer(which(bold_lgl)),
                                                                                                          digits_dbl = c(ifelse(output_type_1L_chr == "PDF",
                                                                                                                                0, NA_real_) %>% purrr::discard(is.na), names(data_tb) %>%
                                                                                                                             purrr::map_dbl(~ifelse(.x %in% c("Bulk_ESS",
                                                                                                                                                              "Tail_ESS"), 0, digits_1L_dbl))), big_mark_1L_chr = big_mark_1L_chr,
                                                                                                          hline_after_ls = c(-1, 0), sanitize_fn = force, footnotes_chr = NA_character_),
                              part_6 = end_matter_1L_chr)
    if (output_type_1L_chr != "PDF") {
        brms_mdl_print_ls$part_5$footnotes_chr <- c(paste0(brms_mdl_print_ls$part_1,
                                                           ifelse(output_type_1L_chr == "Word", "", "\n")),
                                                    brms_mdl_print_ls$part_3, brms_mdl_print_ls$part_6)
        brms_mdl_print_ls$part_6 <- NULL
    }
    else {
        footnotes_chr <- c(mdl_smry_chr[idx_dbl[1]], trimws(mdl_smry_chr[1:(idx_dbl[2] -
                                                                                1)][-idx_dbl[1]]), trimws(mdl_smry_chr[idx_dbl[2]]),
                           trimws(mdl_smry_chr[idx_dbl[2] + 1]), trimws(mdl_smry_chr[idx_dbl[6]:length(mdl_smry_chr)])) %>%
            Hmisc::latexTranslate()
        footnotes_chr[1] <- footnotes_chr[1] %>% stringr::str_replace("~",
                                                                      "\\\\textasciitilde")
        brms_mdl_print_ls$part_5$add_to_row_ls <- list(pos = purrr::map(c(0,
                                                                     rep(nrow(data_tb), length(footnotes_chr) + 1)), ~.x),
                                                  command = c(names(data_tb) %>% Hmisc::latexTranslate() %>%
                                                                  paste0(collapse = " & ") %>% paste0("\\\\\n"),
                                                              c("\\toprule\n"), footnotes_chr %>% purrr::map_chr(~paste0("\\multicolumn{",
                                                                                                                         ncol(data_tb), "}{l}{", paste0("{\\footnotesize ",
                                                                                                                                                        .x, "}\n", collapse = ","), "}\\\\\n"))))
    }
    return(brms_mdl_print_ls)
}
make_brms_mdl_smry_tbl <- function (smry_mdl_ls, grp_1L_chr, pop_1L_chr, fam_1L_chr)
{
    brms_mdl_smry_tb <- purrr::map(1:length(smry_mdl_ls$random),
                                   ~make_mdl_smry_elmt_tbl(cat_chr = c(ifelse(.x == 1, grp_1L_chr,
                                                                              character(0)), paste0(names(smry_mdl_ls$ngrps)[.x],
                                                                                                    " (Number of levels: ", smry_mdl_ls$ngrps[.x][[1]],
                                                                                                    ")")), mat = smry_mdl_ls$random[.x][[1]])) %>% dplyr::bind_rows(make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$fixed,
                                                                                                                                                                                           cat_chr = pop_1L_chr), make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$spec_pars,
                                                                                                                                                                                                                                         cat_chr = fam_1L_chr))
    return(brms_mdl_smry_tb)
}
make_complete_props_tbs_ls <- function (raw_props_tbs_ls, question_var_nm_1L_chr = "Question")
{
  complete_props_tbs_ls <- raw_props_tbs_ls %>% purrr::map(~{
    .x %>% dplyr::mutate(total_prop_dbl = rowSums(dplyr::select(.,
                                                                -!!rlang::sym(question_var_nm_1L_chr)), na.rm = T) -
                           100) %>% dplyr::mutate_if(is.numeric, ~purrr::map2_dbl(.,
                                                                                  total_prop_dbl, ~ifelse(.x == 100, 1 - .y, .x))) %>%
      dplyr::select(-total_prop_dbl)
  })
  return(complete_props_tbs_ls)
}
make_correlated_data_tb <- function (synth_data_spine_ls, synth_data_idx_1L_dbl = 1)
{
  correlated_data_tb <- simstudy::genCorData(synth_data_spine_ls$nbr_obs_dbl[synth_data_idx_1L_dbl],
                                             mu = synth_data_spine_ls$means_ls[[synth_data_idx_1L_dbl]],
                                             sigma = synth_data_spine_ls$sds_ls[[synth_data_idx_1L_dbl]],
                                             corMatrix = make_pdef_cor_mat_mat(synth_data_spine_ls$cor_mat_ls[[synth_data_idx_1L_dbl]]),
                                             cnames = synth_data_spine_ls$var_names_chr) %>% force_min_max_and_int_cnstrs(var_names_chr = synth_data_spine_ls$var_names_chr,
                                                                                                                          min_max_ls = synth_data_spine_ls$min_max_ls, discrete_lgl = synth_data_spine_ls$discrete_lgl)
  return(correlated_data_tb)
}
make_corstars_tbl_xx <- function (x, method = c("pearson", "spearman"), removeTriangle = c("upper",
                                                                                           "lower"), result = c("none", "html", "latex"))
{
    x <- as.matrix(x)
    correlation_matrix <- Hmisc::rcorr(x, type = method[1])
    R <- correlation_matrix$r
    p <- correlation_matrix$P
    mystars <- ifelse(p < 1e-04, "****", ifelse(p < 0.001, "*** ",
                                                ifelse(p < 0.01, "**  ", ifelse(p < 0.05, "*   ", "    "))))
    R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[, -1]
    Rnew <- matrix(paste(R, mystars, sep = ""), ncol = ncol(x))
    diag(Rnew) <- paste(diag(R), " ", sep = "")
    rownames(Rnew) <- colnames(x)
    colnames(Rnew) <- paste(colnames(x), "", sep = "")
    if (removeTriangle[1] == "upper") {
        Rnew <- as.matrix(Rnew)
        Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
        Rnew <- as.data.frame(Rnew)
    }
    else if (removeTriangle[1] == "lower") {
        Rnew <- as.matrix(Rnew)
        Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
        Rnew <- as.data.frame(Rnew)
    }
    Rnew <- cbind(Rnew[1:length(Rnew) - 1])
    if (result[1] == "none")
        return(Rnew)
    else {
        if (result[1] == "html")
            print(xtable(Rnew), type = "html")
        else print(xtable(Rnew), type = "latex")
    }
}
make_dim_sclg_cons_dbl <- function (domains_chr, dim_sclg_con_lup_tb)
{
  dim_sclg_cons_dbl <- purrr::map_dbl(domains_chr, ~ready4fun::get_from_lup_obj(dim_sclg_con_lup_tb,
                                                                                match_var_nm_1L_chr = "Dimension_chr", match_value_xx = .x,
                                                                                target_var_nm_1L_chr = "Constant_dbl", evaluate_lgl = F))
  return(dim_sclg_cons_dbl)
}
make_domain_items_ls <- function (domain_qs_lup_tb, item_pfx_1L_chr)
{
  domains_chr <- domain_qs_lup_tb$Domain_chr %>% unique()
  q_nbrs_ls <- purrr::map(domains_chr, ~domain_qs_lup_tb %>%
                            dplyr::filter(Domain_chr == .x) %>% dplyr::pull(Question_dbl))
  domain_items_ls <- purrr::map(q_nbrs_ls, ~paste0(item_pfx_1L_chr,
                                                   .x)) %>% stats::setNames(domains_chr)
  return(domain_items_ls)
}
make_fake_ts_data <- function (outp_smry_ls)
{
    data_tb <- outp_smry_ls$scored_data_tb %>% transform_tb_to_mdl_inp(dep_var_nm_1L_chr = outp_smry_ls$dep_var_nm_1L_chr,
        predr_vars_nms_chr = outp_smry_ls$predr_cmprsns_tb$predr_chr,
        id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
        round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr)
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = names(data_tb)[names(data_tb) !=
        outp_smry_ls$id_var_nm_1L_chr], seed = outp_smry_ls$seed_1L_int)
    dep_vars_chr <- names(fk_data_ls$syn)[names(fk_data_ls$syn) %>%
        purrr::map_lgl(~startsWith(.x, outp_smry_ls$dep_var_nm_1L_chr))]
    fk_data_tb <- fk_data_ls$syn %>% dplyr::mutate(dplyr::across(dplyr::all_of(dep_vars_chr),
        ~NA_real_))
    return(fk_data_tb)
}
make_folds_ls <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", n_folds_1L_int = 10L)
{
    folds_ls <- caret::createFolds(data_tb %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr)),
                                   k = n_folds_1L_int, list = TRUE, returnTrain = FALSE)
    return(folds_ls)
}
make_item_wrst_wghts_ls_ls <- function (domain_items_ls, itm_wrst_wghts_lup_tb)
{
  item_wrst_wghts_ls_ls <- domain_items_ls %>% purrr::map(~{
    purrr::map_dbl(.x, ~{
      ready4fun::get_from_lup_obj(itm_wrst_wghts_lup_tb,
                                  match_var_nm_1L_chr = "Question_chr", match_value_xx = .x,
                                  target_var_nm_1L_chr = "Worst_Weight_dbl", evaluate_lgl = F)
    })
  })
  return(item_wrst_wghts_ls_ls)
}
make_knit_pars_ls <- function (rltv_path_to_data_dir_1L_chr, mdl_types_chr, predr_vars_nms_ls,
                               output_type_1L_chr = "HTML", mdl_types_lup = NULL, plt_types_lup = NULL,
                               plt_types_chr = NA_character_, section_type_1L_chr = "#")
{
    if (is.null(mdl_types_lup))
        utils::data(mdl_types_lup, envir = environment())
    if (is.null(plt_types_lup))
        utils::data(plt_types_lup, envir = environment())
  if(is.na(plt_types_chr)){
    plt_types_chr <- plt_types_lup$short_name_chr
  }
   paths_to_all_data_fls_chr <- list.files(rltv_path_to_data_dir_1L_chr,full.names = T)
    lab_idx_dbl <- 1:(length(mdl_types_chr) * length(predr_vars_nms_ls))
    knit_pars_ls <- purrr::pmap(list(predr_vars_nms_ls,
                                     split(lab_idx_dbl,
                                           ceiling(seq_along(lab_idx_dbl)/length(mdl_types_chr))),
                                     #make_unique_ls_elmt_idx_int(predr_vars_nms_ls), # POSSIBLY REMOVE AND UPDATE PURRR INDICES (..4 to ..3)
                                     make_mdl_nms_ls(predr_vars_nms_ls,
                                                     mdl_types_chr = mdl_types_chr)),
                                ~{
                                  mdl_nms_chr <- ..3#..4
                                  mdl_data_paths_ls <- mdl_nms_chr %>%
                                    purrr::map(~paths_to_all_data_fls_chr[stringr::str_detect(paths_to_all_data_fls_chr,.x)]) %>%
                                    stats::setNames(mdl_nms_chr)
                                  paths_to_mdls_chr <- mdl_data_paths_ls %>%
                                    purrr::map_chr(~.x[endsWith(.x,".RDS")]) %>% unname()
                                  paths_to_mdl_plts_ls <- mdl_data_paths_ls %>%
                                    purrr::map(~{
                                      paths_to_all_plots_chr <- .x[endsWith(.x,".png")]
                                      plt_types_chr %>% purrr::map(~paths_to_all_plots_chr[paths_to_all_plots_chr %>%
                                                                                                stringr::str_detect(.x)]) %>% purrr::flatten_chr()
                                      })
                                  mdl_ttls_chr <- paste0(..1[1], ifelse(is.na(..1[2]),
                                                                        "",
                                                                        paste(" with ", ..1[2])),
                                                         " ",
                                                         mdl_types_chr %>% purrr::map_chr(~ready4fun::get_from_lup_obj(mdl_types_lup,
                                                                                     match_var_nm_1L_chr = "short_name_chr",
                                                                                     match_value_xx = .x,
                                                                                     target_var_nm_1L_chr = "long_name_chr",
                                                                                     evaluate_lgl = F)))
                                  section_ttls_chr <- paste0(section_type_1L_chr, " ", mdl_ttls_chr)
                                  plt_nms_ls <- paths_to_mdl_plts_ls %>% purrr::map2(mdl_ttls_chr,~{
                                    paths_to_mdl_plts_chr <- .x
                                    mdl_ttl_1L_chr <- .y
                                    transform_1L_lgl <- paths_to_mdl_plts_chr %>% endsWith("_coefs.png") %>% any()
                                    if(paths_to_mdl_plts_chr %>% endsWith("_hetg.png") %>% any())
                                      transform_1L_lgl <- F # ALL OF THIS TFMN LOGIC CAN BE BINNED WHEN NEW CATEGORY OF PLOT TYPE IS ADDED.
                                    plt_types_chr %>% purrr::map(~{
                                      if(endsWith(paths_to_mdl_plts_chr,paste0("_",.x,".png")) %>% any()){
                                        paste0(mdl_ttl_1L_chr,
                                               " ",
                                               ifelse(transform_1L_lgl & .x == "coefs",
                                                      "population and group level effects",
                                                      ready4fun::get_from_lup_obj(plt_types_lup,
                                                                           match_var_nm_1L_chr = "short_name_chr",
                                                                           match_value_xx = .x,
                                                                           target_var_nm_1L_chr = "long_name_chr",
                                                                           evaluate_lgl = F)))

                                      }else{
                                        character(0)
                                      }}) %>% purrr::flatten_chr()
                                  })
                                  list(plt_nms_ls = plt_nms_ls,
                                       paths_to_mdls_chr = paths_to_mdls_chr,
                                       tbl_captions_chr = mdl_ttls_chr,
                                       label_stubs_chr = paste0("lab", ..2),
                                       output_type_1L_chr = output_type_1L_chr,
                                       section_ttls_chr = section_ttls_chr,
                                       #ls_elmt_idx_1L_int = ..3,
                                       paths_to_mdl_plts_ls = paths_to_mdl_plts_ls)
                                  }) %>%
      stats::setNames(predr_vars_nms_ls %>% purrr::map_chr(~paste(.x,collapse="_")))
    return(knit_pars_ls)
}
make_mdl <- function (data_tb, dep_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF",
    predr_var_nm_1L_chr, covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF",
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- transform_ds_for_mdlng(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    if (is.null(start_1L_chr)) {
        start_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
            target_var_nm_1L_chr = "start_chr", evaluate_lgl = F)
    }
    if (!is.na(control_1L_chr)) {
        idx_1L_int <- 1 + stringi::stri_locate_last_fixed(mdl_type_1L_chr,
            "_")[1, 1] %>% as.vector()
        link_1L_chr <- stringr::str_sub(mdl_type_1L_chr, start = idx_1L_int)
        link_1L_chr <- ifelse(link_1L_chr == "LOG", "log", ifelse(link_1L_chr ==
            "LGT", "logit", ifelse(link_1L_chr == "CLL", "cloglog",
            "ERROR")))
    }
    mdl_1L_chr <- paste0(ready4fun::get_from_lup_obj(mdl_types_lup,
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
        target_var_nm_1L_chr = "fn_chr", evaluate_lgl = F), "(",
        transform_dep_var_nm(dep_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr),
        " ~ ", predr_var_nm_1L_chr, ifelse(is.na(covar_var_nms_chr[1]),
            "", paste0(" + ", paste0(covar_var_nms_chr, collapse = " + "))),
        ", data = data_tb", ifelse(!is.na(ready4fun::get_from_lup_obj(mdl_types_lup,
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
            target_var_nm_1L_chr = "family_chr", evaluate_lgl = F)),
            paste0(", family = ", ready4fun::get_from_lup_obj(mdl_types_lup,
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                target_var_nm_1L_chr = "family_chr", evaluate_lgl = F)),
            ""), ifelse(!is.na(start_1L_chr), ", ", ""), ifelse(!is.na(control_1L_chr),
            paste0("link=\"", link_1L_chr, "\",control=", control_1L_chr,
                "("), ""), ifelse(!is.na(start_1L_chr), paste0("start=c(",
            start_1L_chr, ")"), ""), ifelse(!is.na(control_1L_chr),
            ")", ""), ")")
    model_mdl <- eval(parse(text = mdl_1L_chr))
    return(model_mdl)
}
make_mdl_nms_ls <- function (predr_vars_nms_ls, mdl_types_chr)
{
    mdl_nms_ls <- purrr::map2(predr_vars_nms_ls, make_unique_ls_elmt_idx_int(predr_vars_nms_ls),
        ~paste0(.x[1], "_", ifelse(is.na(.x[2]), "", paste0(.x[2],
            "_")), .y, "_", mdl_types_chr))
    return(mdl_nms_ls)
}
make_mdl_smry_elmt_tbl <- function (mat, cat_chr)
{
    tb <- mat %>% tibble::as_tibble() %>% dplyr::mutate(Parameter = rownames(mat)) %>%
        dplyr::select(Parameter, dplyr::everything())
    mdl_elmt_sum_tb <- tb %>% dplyr::filter(F) %>% tibble::add_case(Parameter = cat_chr) %>%
        dplyr::bind_rows(tb)
    return(mdl_elmt_sum_tb)
}
make_pdef_cor_mat_mat <- function (lower_diag_mat)
{
  pdef_cor_mat <- lower_diag_mat %>% Matrix::forceSymmetric(uplo = "L") %>%
    as.matrix()
  if (!matrixcalc::is.positive.definite(pdef_cor_mat)) {
    pdef_cor_mat <- psych::cor.smooth(pdef_cor_mat)
  }
  return(pdef_cor_mat)
}
make_predn_ds_with_one_predr <- function (model_mdl, dep_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF",
    predr_var_nm_1L_chr, predr_vals_dbl, pred_type_1L_chr = NULL)
{
    predn_ds_tb <- tibble::tibble(`:=`(!!rlang::sym(predr_var_nm_1L_chr),
        predr_vals_dbl))
    predn_ds_tb <- predn_ds_tb %>% dplyr::mutate(`:=`(!!rlang::sym(dep_var_nm_1L_chr),
        stats::predict(model_mdl, newdata = predn_ds_tb, type = pred_type_1L_chr) %>%
            calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
                tfmn_is_outp_1L_lgl = T)))
    return(predn_ds_tb)
}
make_predr_vals <- function (predr_var_nm_1L_chr, candidate_predrs_lup = NULL)
{
    if (is.null(candidate_predrs_lup)) {
        utils::data("candidate_predrs_lup", envir = environment())
    }
    args_ls <- purrr::map_dbl(names(candidate_predrs_lup)[3:4],
                              ~candidate_predrs_lup %>% ready4fun::get_from_lup_obj(match_value_xx = predr_var_nm_1L_chr,
                                                                                    match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = .x,
                                                                                    evaluate_lgl = F)) %>% as.list()
    predr_vals_dbl <- rlang::exec(.fn = seq, !!!args_ls)
    return(predr_vals_dbl)
}
make_predr_vars_nms_ls <- function (main_predrs_chr, covars_ls)
{
    predr_vars_nms_ls <- covars_ls %>% purrr::map(~{
        covars_chr <- .x
        purrr::map(main_predrs_chr, ~list(c(.x), c(.x, covars_chr))) %>%
            purrr::flatten()
    }) %>% purrr::flatten()
    predr_vars_nms_ls <- predr_vars_nms_ls[order(sapply(predr_vars_nms_ls,
                                                        length))]
    return(predr_vars_nms_ls)
}
make_prefd_mdls_vec <- function (smry_of_sngl_predr_mdls_tb, choose_from_pfx_chr = c("GLM",
                                                                                     "OLS"), mdl_types_lup = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    ordered_mdl_types_chr <- dplyr::inner_join(smry_of_sngl_predr_mdls_tb %>%
                                                   dplyr::select(Model) %>% dplyr::rename(long_name_chr = Model),
                                               mdl_types_lup) %>% dplyr::pull(short_name_chr)
    prefd_mdls_chr <- purrr::map_chr(choose_from_pfx_chr, ~ordered_mdl_types_chr[startsWith(ordered_mdl_types_chr,
                                                                                            .x)][1])
    return(prefd_mdls_chr)
}
make_shareable_mdl <- function (data_tb, mdl_smry_tb, dep_var_nm_1L_chr = "utl_total_w",
    id_var_nm_1L_chr = "fkClientID", tfmn_1L_chr = "CLL", mdl_type_1L_chr = "OLS_CLL",
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NA_character_,
    seed_1L_int = 12345L)
{
    if (is.null(mdl_types_lup))
        utils::data(mdl_types_lup, envir = environment())
    all_var_nms_chr <- names(data_tb)
    tfd_dep_var_nm_1L_chr <- all_var_nms_chr[all_var_nms_chr %>%
        startsWith(dep_var_nm_1L_chr)]
    predr_var_nms_chr <- setdiff(all_var_nms_chr, c(tfd_dep_var_nm_1L_chr,
        id_var_nm_1L_chr))
    if (length(predr_var_nms_chr) > 1) {
        covar_var_nms_chr <- predr_var_nms_chr[2:length(predr_var_nms_chr)]
    }else{
        covar_var_nms_chr <- NA_character_
    }
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = all_var_nms_chr[all_var_nms_chr !=
        id_var_nm_1L_chr], seed = seed_1L_int)
    model_mdl <- make_mdl(fk_data_ls$syn, dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        predr_var_nm_1L_chr = predr_var_nms_chr[1], covar_var_nms_chr = covar_var_nms_chr,
        tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr,
        mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr,
        start_1L_chr = start_1L_chr)
    par_nms_chr <- model_mdl$coefficients %>% names()
    mdl_smry_tb <- mdl_smry_tb %>% dplyr::mutate(Parameter = dplyr::case_when(Parameter ==
        "Intercept" ~ "(Intercept)", TRUE ~ purrr::map_chr(Parameter,
        ~stringr::str_replace_all(.x, " ", "_")))) %>% dplyr::filter(Parameter %in%
        par_nms_chr) %>% dplyr::slice(match(par_nms_chr, Parameter))
    assertthat::assert_that(all(par_nms_chr == mdl_smry_tb$Parameter),
        msg = "Parameter names mismatch between data and model summary table")
    model_mdl$coefficients <- mdl_smry_tb$Estimate
    names(model_mdl$coefficients) <- par_nms_chr
    return(model_mdl)
}
make_smry_of_brm_mdl <- function (mdl_ls, data_tb, dep_var_nm_1L_chr = "utl_total_w",
    predr_vars_nms_chr, fn = calculate_rmse, mdl_nm_1L_chr = NA_character_,
    seed_1L_dbl = 23456)
{
    if (is.na(mdl_nm_1L_chr))
        mdl_nm_1L_chr <- predr_vars_nms_chr[1]
    set.seed(seed_1L_dbl)
    predictions <- stats::predict(mdl_ls, summary = F)
    coef <- summary(mdl_ls, digits = 4)$fixed
    coef <- coef[1:nrow(coef), 1:4]
    R2 <- brms::bayes_R2(mdl_ls)
    RMSE <- psych::describe(apply(predictions, 1, fn, y_dbl = data_tb %>%
        dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))), quant = c(0.25,
        0.75), skew = F, ranges = F)
    RMSE <- cbind(RMSE$mean, RMSE$sd, RMSE$Q0.25, RMSE$Q0.75) %>%
        as.vector()
    Sigma <- summary(mdl_ls, digits = 4)$spec_par[1:4]
    smry_of_brm_mdl_tb <- data.frame(round(rbind(coef, R2, RMSE,
        Sigma), 3)) %>% dplyr::mutate(Parameter = c("Intercept",
        purrr::map(predr_vars_nms_chr, ~paste0(.x, c(" baseline",
            " change"))) %>% purrr::flatten_chr(), "R2", "RMSE",
        "Sigma"), Model = mdl_nm_1L_chr) %>% dplyr::mutate(`95% CI` = paste(l.95..CI,
        ",", u.95..CI)) %>% dplyr::rename(SE = Est.Error) %>%
        dplyr::select(Model, Parameter, Estimate, SE, `95% CI`)
    return(smry_of_brm_mdl_tb)
}
make_smry_of_mdl <- function (data_tb, model_mdl, n_folds_1L_int = 10, dep_var_nm_1L_chr = "utl_total_w",
    start_1L_chr = NULL, tfmn_1L_chr = "NTF", predr_var_nm_1L_chr,
    covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF",
    mdl_types_lup = NULL, pred_type_1L_chr = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- data_tb %>% dplyr::filter(!is.na(!!rlang::sym(predr_var_nm_1L_chr)))
    data_tb <- transform_ds_for_mdlng(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    mdl_desc_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
        target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F)
    folds_ls <- make_folds_ls(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        n_folds_1L_int = n_folds_1L_int)
    smry_of_one_predr_mdl_tb <- purrr::map_dfr(folds_ls, ~{
        model_mdl <- make_mdl(data_tb[-.x,], dep_var_nm_1L_chr = dep_var_nm_1L_chr,
            start_1L_chr = start_1L_chr, tfmn_1L_chr = tfmn_1L_chr,
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr,
            mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup)
        pred_old_dbl <- stats::predict(model_mdl, type = pred_type_1L_chr)
        pred_new_dbl <- stats::predict(model_mdl, newdata = data_tb[.x, ],
            type = pred_type_1L_chr) %>% calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
            tfmn_is_outp_1L_lgl = T)
        tibble::tibble(Rsquared = caret::R2(pred_old_dbl, data_tb[-.x,
            ] %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr)),
            form = "traditional"), RMSE = caret::RMSE(pred_old_dbl,
            data_tb[-.x, ] %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))),
            MAE = caret::MAE(pred_old_dbl, data_tb[-.x, ] %>%
                dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))),
            RsquaredP = caret::R2(pred_new_dbl, data_tb[.x, ] %>%
                dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr)),
                form = "traditional"), RMSEP = caret::RMSE(pred_new_dbl,
                data_tb[.x, ] %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))),
            MAEP = caret::MAE(pred_new_dbl, data_tb[.x, ] %>%
                dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr))))
    }) %>% dplyr::summarise_all(mean) %>% dplyr::mutate(Model = mdl_desc_1L_chr) %>%
        dplyr::select(Model, dplyr::everything())
    return(smry_of_one_predr_mdl_tb)
}
make_smry_of_ts_mdl <- function (data_tb, fn, predr_vars_nms_chr, mdl_nm_1L_chr, path_to_write_to_1L_chr = NA_character_,
    dep_var_nm_1L_chr = "utl_total_w", id_var_nm_1L_chr = "fkClientID",
    round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline",
    backend_1L_chr = getOption("brms.backend", "rstan"), iters_1L_int = 4000L,
    seed_1L_int = 1000L)
{
    tfd_data_tb <- transform_tb_to_mdl_inp(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr,
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr)
    tfd_dep_var_nm_1L_chr <- ifelse(identical(fn, fit_clg_log_tfmn),
        transform_dep_var_nm_for_cll(dep_var_nm_1L_chr), dep_var_nm_1L_chr)
    args_ls <- list(data_tb = tfd_data_tb, dep_var_nm_1L_chr = tfd_dep_var_nm_1L_chr,
        predr_vars_nms_chr = predr_vars_nms_chr, iters_1L_int = iters_1L_int,
        backend_1L_chr = backend_1L_chr, seed_1L_int = seed_1L_int)
    mdl_ls <- rlang::exec(fn, !!!args_ls)
    smry_of_ts_mdl_ls <- list(smry_of_ts_mdl_tb = make_smry_of_brm_mdl(mdl_ls,
        data_tb = tfd_data_tb, dep_var_nm_1L_chr = tfd_dep_var_nm_1L_chr,
        predr_vars_nms_chr = predr_vars_nms_chr, fn = ifelse(identical(fn,
            fit_gsn_log_lnk), calculate_rmse, calculate_rmse_tfmn),
        mdl_nm_1L_chr = mdl_nm_1L_chr))
    if (!is.na(path_to_write_to_1L_chr)) {
        smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr <- paste0(path_to_write_to_1L_chr,
            "/", mdl_nm_1L_chr, ".RDS")
        if (file.exists(smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr))
            file.remove(smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)
        saveRDS(mdl_ls, smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)
        smry_of_ts_mdl_ls$paths_to_mdl_plts_chr <- write_brm_model_plts(mdl_ls,
            tfd_data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr,
            mdl_nm_1L_chr = mdl_nm_1L_chr, path_to_write_to_1L_chr = path_to_write_to_1L_chr,
            round_var_nm_1L_chr = round_var_nm_1L_chr, tfmn_fn = ifelse(identical(fn,
                fit_gsn_log_lnk), function(x) {
                x
            }, function(x) {
                1 - exp(-exp(x))
            }))
    }
    return(smry_of_ts_mdl_ls)
}
make_synth_series_tbs_ls <- function (synth_data_spine_ls, series_names_chr)
{
  synth_series_tbs_ls <- 1:length(series_names_chr) %>% purrr::map(~make_correlated_data_tb(synth_data_spine_ls = synth_data_spine_ls,
                                                                                            synth_data_idx_1L_dbl = .x) %>% replace_with_missing_vals(synth_data_spine_ls = synth_data_spine_ls,
                                                                                                                                                      idx_int = .x)) %>% stats::setNames(series_names_chr)
  return(synth_series_tbs_ls)
}
make_vec_with_sum_of_int <- function (target_int, start_int, end_int, length_int)
{
  vec_int <- Surrogate::RandVec(a = start_int, b = end_int,
                                s = target_int, n = length_int, m = 1) %>% purrr::pluck("RandVecOutput") %>%
    as.vector() %>% round() %>% as.integer() %>% force_vec_to_sum_to_int(target_1L_int = target_int)
  return(vec_int)
}
make_unique_ls_elmt_idx_int <- function (data_ls)
{
  combos_tb <- tibble::as_tibble(data_ls, .name_repair = ~paste0("r_",
                                                                 1:length(data_ls))) %>% t() %>% as.data.frame()
  combos_tb <- combos_tb %>% tibble::as_tibble()
  combos_tb <- combos_tb %>% dplyr::mutate(V2 = dplyr::case_when(V1 ==
                                                                   V2 ~ NA_character_, T ~ V2)) %>% dplyr::mutate(combo_chr = purrr::map2_chr(V1,
                                                                                                                                              V2, ~ifelse(ncol(combos_tb) == 1 | is.na(.y), .x, paste0(.x,
                                                                                                                                                                                                       "_", .y))))
  combos_tb <- combos_tb %>% dplyr::group_by(combo_chr) %>%
    dplyr::mutate(combo_id = dplyr::row_number())
  unique_ls_elmt_idx_int <- purrr::map(data_ls %>% unique(),
                                       ~ready4fun::get_from_lup_obj(combos_tb %>% dplyr::ungroup(),
                                                                    match_var_nm_1L_chr = "combo_chr", match_value_xx = paste0(.x[1],
                                                                                                                               ifelse(is.na(.x[2]), "", paste0("_", .x[2]))),
                                                                    target_var_nm_1L_chr = "combo_id", evaluate_lgl = F)) %>%
    purrr::flatten_int()
  return(unique_ls_elmt_idx_int)
}
