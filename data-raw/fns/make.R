make_analysis_ds_smry_ls <- function(ds_descvs_ls,
                                     candidate_covar_nms_chr,
                                     predictors_lup){
  analysis_ds_smry_ls <- list(candidate_predrs_chr = ds_descvs_ls$candidate_predrs_chr,
                              candidate_covar_nms_chr = candidate_covar_nms_chr,
                              depnt_var_nm_1L_chr = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                              id_var_nm_1L_chr = ds_descvs_ls$id_var_nm_1L_chr,
                              predictors_lup = predictors_lup,
                              round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                              round_bl_val_1L_chr = ds_descvs_ls$round_vals_chr[1],
                              dictionary_tb = ds_descvs_ls$dictionary_tb)

  return(analysis_ds_smry_ls)
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
                                      popl_1L_chr = mdl_smry_chr[idx_dbl[4]], fam_1L_chr = mdl_smry_chr[idx_dbl[5]])
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
make_brms_mdl_smry_tbl <- function (smry_mdl_ls, grp_1L_chr, popl_1L_chr, fam_1L_chr)
{
    brms_mdl_smry_tb <- purrr::map(1:length(smry_mdl_ls$random),
                                   ~make_mdl_smry_elmt_tbl(ctg_chr = c(ifelse(.x == 1, grp_1L_chr,
                                                                              character(0)), paste0(names(smry_mdl_ls$ngrps)[.x],
                                                                                                    " (Number of levels: ", smry_mdl_ls$ngrps[.x][[1]],
                                                                                                    ")")), mat = smry_mdl_ls$random[.x][[1]])) %>% dplyr::bind_rows(make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$fixed,
                                                                                                                                                                                           ctg_chr = popl_1L_chr), make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$spec_pars,
                                                                                                                                                                                                                                         ctg_chr = fam_1L_chr))
    return(brms_mdl_smry_tb)
}
make_cmpst_sctr_and_dnsty_plt <- function(outp_smry_ls,
                                          predr_var_nms_chr,
                                          labels_chr = c("A","B","C","D"),
                                          label_x_1L_dbl = 0.1,
                                          label_y_1L_dbl = 0.9,
                                          label_size_1L_dbl = 22){
  plot_ls <- paste0(output_data_dir_1L_chr,"/",outp_smry_ls$file_paths_chr[outp_smry_ls$file_paths_chr %>% purrr::map_lgl(~stringr::str_detect(.x,paste0(predr_var_nms_chr,"_1")) & (stringr::str_detect(.x,"_dnst.png") | stringr::str_detect(.x,"_sctr_plt.png")))])  %>% purrr::map(~cowplot::ggdraw() + cowplot::draw_image(.x))
  composite_plt <- cowplot::plot_grid(plot_ls[[1]],plot_ls[[2]],plot_ls[[3]],plot_ls[[4]],nrow = 2, labels = labels_chr, label_x = label_x_1L_dbl,label_y = label_y_1L_dbl, label_size = label_size_1L_dbl)
}
make_fake_ts_data <- function (outp_smry_ls)
{
    data_tb <- outp_smry_ls$scored_data_tb %>% transform_tb_to_mdl_inp(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
        predr_vars_nms_chr = outp_smry_ls$predr_cmprsn_tb$predr_chr,
        id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
        round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr)
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = names(data_tb)[names(data_tb) !=
        outp_smry_ls$id_var_nm_1L_chr], seed = outp_smry_ls$seed_1L_int)
    dep_vars_chr <- names(fk_data_ls$syn)[names(fk_data_ls$syn) %>%
        purrr::map_lgl(~startsWith(.x, outp_smry_ls$depnt_var_nm_1L_chr))]
    fk_data_tb <- fk_data_ls$syn %>% dplyr::mutate(dplyr::across(dplyr::all_of(dep_vars_chr),
        ~NA_real_))
    return(fk_data_tb)
}
make_folds_ls <- function (data_tb, depnt_var_nm_1L_chr = "aqol6d_total_w", folds_1L_int = 10L)
{
    folds_ls <- caret::createFolds(data_tb %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
                                   k = folds_1L_int, list = TRUE, returnTrain = FALSE)
    return(folds_ls)
}
make_hlth_utl_and_predrs_ls <- function(outp_smry_ls,
                                        descv_tbls_ls,
                                        nbr_of_digits_1L_int = 2L,
                                        old_nms_chr = NULL,
                                        new_nms_chr = NULL){
  ranked_predrs_ls <- make_ranked_predrs_ls(descv_tbls_ls,
                                            old_nms_chr = old_nms_chr,
                                            new_nms_chr = new_nms_chr)
  var_nm_1L_chr <- descv_tbls_ls$ds_descvs_ls$dictionary_tb %>%
    ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "var_nm_chr",
                                match_value_xx = outp_smry_ls$depnt_var_nm_1L_chr,
                                target_var_nm_1L_chr = "var_desc_chr",
                                evaluate_lgl = F) %>% as.vector()


  hlth_utl_and_predrs_ls = list(bl_hu_mean_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                              match_value_xx = var_nm_1L_chr,
                                                              target_var_nm_1L_chr = "Baseline_val_1_dbl",
                                                              evaluate_lgl = F) %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                bl_hu_sd_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                              match_value_xx = var_nm_1L_chr,
                                                              target_var_nm_1L_chr = "Baseline_val_2_ls",
                                                              evaluate_lgl = F) %>%
                                  stringr::str_remove("\\(") %>%
                                  stringr::str_remove("\\)") %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                fup_hu_mean_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                              match_value_xx = var_nm_1L_chr,
                                                              target_var_nm_1L_chr = "Follow-up_val_1_dbl",
                                                              evaluate_lgl = F) %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                fup_hu_sd_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                              match_value_xx = var_nm_1L_chr,
                                                              target_var_nm_1L_chr = "Follow-up_val_2_ls",
                                                              evaluate_lgl = F) %>%
                                  stringr::str_remove("\\(") %>%
                                  stringr::str_remove("\\)") %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                predrs_nartv_seq_chr = ranked_predrs_ls$unranked_predrs_chr,
                                cor_seq_dscdng_chr =  ranked_predrs_ls$ranked_predrs_chr)
  return(hlth_utl_and_predrs_ls)

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
make_mdl <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF",
    predr_var_nm_1L_chr, covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF",
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
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
        transform_depnt_var_nm(depnt_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr),
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
make_mdl_coef_ratio_ls <- function(outp_smry_ls,
                                   predr_ctgs_ls = NULL){
  main_mdls_ls <- outp_smry_ls$predr_cmprsn_tb$predr_chr %>% purrr::map(~paste0(paste0(.x,"_1_"),
                                                                                outp_smry_ls$prefd_mdl_types_chr))
  ratios_ls <- main_mdls_ls %>%
    purrr::map2(outp_smry_ls$predr_cmprsn_tb$predr_chr, ~{
      mdls_chr <- .x
      predr_1L_chr <- .y
      mdls_chr %>% purrr::map_dbl(~{
        coefs_dbl <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Model %in% .x) %>%
          dplyr::filter(Parameter %in% paste0(predr_1L_chr,(c(" baseline"," change")))) %>%
          dplyr::pull(Estimate)
        coefs_dbl[2]/coefs_dbl[1]
      })
    }
    )
  mdl_coef_ratios_ls = list(mean_ratios_dbl = ratios_ls %>% purrr::map_dbl(~mean(.x)))
  if(!is.null(predr_ctgs_ls)){
    append_ls <- purrr::map(predr_ctgs_ls,
                            ~ mdl_coef_ratios_ls$mean_ratios_dbl[outp_smry_ls$predr_cmprsn_tb$predr_chr %in% .x]) %>%
      stats::setNames(names(predr_ctgs_ls))
    mdl_coef_ratios_ls <- append(mdl_coef_ratios_ls, append_ls)
  }
  return(mdl_coef_ratios_ls)
}
make_mdl_nms_ls <- function (predr_vars_nms_ls, mdl_types_chr)
{
    mdl_nms_ls <- purrr::map2(predr_vars_nms_ls, make_unique_ls_elmt_idx_int(predr_vars_nms_ls),
        ~paste0(.x[1], "_", ifelse(is.na(.x[2]), "", paste0(.x[2],
            "_")), .y, "_", mdl_types_chr))
    return(mdl_nms_ls)
}
make_mdl_smry_elmt_tbl <- function (mat, ctg_chr)
{
    tb <- mat %>% tibble::as_tibble() %>% dplyr::mutate(Parameter = rownames(mat)) %>%
        dplyr::select(Parameter, dplyr::everything())
    mdl_elmt_sum_tb <- tb %>% dplyr::filter(F) %>% tibble::add_case(Parameter = ctg_chr) %>%
        dplyr::bind_rows(tb)
    return(mdl_elmt_sum_tb)
}
make_mdls_ls <- function(outp_smry_ls,
                         mdls_tb){
  mdls_chr <- mdls_tb$Model %>% unique()
  mdls_ls <- outp_smry_ls$prefd_mdl_types_chr %>%
    purrr::map(~mdls_chr[mdls_chr %>% endsWith(.x)])
  return(mdls_ls)
}
make_mdls_smry_tbls_ls <- function(outp_smry_ls,
                                   nbr_of_digits_1L_int = 2L){
  mdls_smry_tb <- outp_smry_ls$mdls_smry_tb  %>%
    dplyr::mutate(dplyr::across(c("Estimate","SE"),
                                ~ round(.x,nbr_of_digits_1L_int) %>%
                                  format(nsmall = nbr_of_digits_1L_int))) %>%
    dplyr::mutate(`95% CI` = `95% CI` %>% transform_chr_digit_pairs())
  rownames(mdls_smry_tb) <- NULL
  indpt_predrs_mdls_tb <- mdls_smry_tb %>% dplyr::filter(Model %in% (paste0(outp_smry_ls$predr_cmprsn_tb$predr_chr,"_1_") %>% purrr::map(~paste0(.x,
                                                                                                                                                 outp_smry_ls$prefd_mdl_types_chr)) %>% purrr::flatten_chr()))
  covar_mdls_tb <- mdls_smry_tb %>% dplyr::filter(!Model %in% indpt_predrs_mdls_tb$Model)
  prefd_predr_mdl_smry_tb <- indpt_predrs_mdls_tb %>%
    dplyr::filter(Model %>% purrr::map_lgl(~startsWith(.x, outp_smry_ls$predr_vars_nms_ls[[1]])))
  mdls_smry_tbls_ls <- list(indpt_predrs_mdls_tb = indpt_predrs_mdls_tb,
                            covar_mdls_tb = covar_mdls_tb,
                            prefd_predr_mdl_smry_tb = prefd_predr_mdl_smry_tb)
  return(mdls_smry_tbls_ls)
}
make_mdl_type_smry_tbl <- function(mdls_tb,
                                   mdl_nms_chr,
                                   mdl_type_1L_chr,
                                   add_mdl_nm_sfx_1L_lgl = T){
  mdl_type_smry_tbl_tb <- mdl_nms_chr %>%
    purrr::map_dfr(~make_sngl_mdl_smry_tb(mdls_tb,
                                          mdl_nm_1L_chr = .x,
                                          mdl_type_1L_chr = mdl_type_1L_chr,
                                          add_mdl_nm_sfx_1L_lgl = add_mdl_nm_sfx_1L_lgl))
  return(mdl_type_smry_tbl_tb)
}
make_predn_ds_with_one_predr <- function (model_mdl, depnt_var_nm_1L_chr = "utl_total_w", tfmn_1L_chr = "NTF",
    predr_var_nm_1L_chr, predr_vals_dbl, predn_type_1L_chr = NULL)
{
    predn_ds_tb <- tibble::tibble(`:=`(!!rlang::sym(predr_var_nm_1L_chr),
        predr_vals_dbl))
    predn_ds_tb <- predn_ds_tb %>% dplyr::mutate(`:=`(!!rlang::sym(depnt_var_nm_1L_chr),
        stats::predict(model_mdl, newdata = predn_ds_tb, type = predn_type_1L_chr) %>%
            calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
                tfmn_is_outp_1L_lgl = T)))
    return(predn_ds_tb)
}
make_predrs_for_best_mdls <- function(outp_smry_ls,
                                      old_nms_chr = NULL,
                                      new_nms_chr){
  ordered_mdl_nms_chr <- outp_smry_ls$mdls_smry_tb %>%
    dplyr::filter(Parameter=="R2") %>%
    dplyr::arrange(dplyr::desc(Estimate)) %>% dplyr::pull(Model)
  ind_predr_mdls_by_mdl_type_ls <- outp_smry_ls$prefd_mdl_types_chr %>%
    purrr::map(~{
    mdl_type_1L_chr <- .x
    paste0(outp_smry_ls$predr_cmprsn_tb$predr_chr,"_1_") %>% purrr::map(~paste0(.x,
                                                                                mdl_type_1L_chr)) %>% purrr::flatten_chr()
  })
  ordered_mdls_by_type_ls <- ind_predr_mdls_by_mdl_type_ls %>%
    purrr::map(~{
      ind_predr_mdls_chr <- .x
      ordered_mdl_nms_chr[ordered_mdl_nms_chr %>% purrr::map_lgl(~.x %in% ind_predr_mdls_chr)]
    })
  predrs_for_best_mdls_chr <- ordered_mdls_by_type_ls %>% purrr::map2_chr(outp_smry_ls$prefd_mdl_types_chr,
                                                                          ~{
                                                                            .x[1] %>% stringr::str_remove(paste0("_1_",.y))
                                                                          })
  if(!is.null(old_nms_chr)){
    predrs_for_best_mdls_chr <- transform_predr_nm_part_of_phrases(predrs_for_best_mdls_chr,
                                                                   old_nms_chr = old_nms_chr,
                                                                   new_nms_chr = new_nms_chr)
  }
  return(predrs_for_best_mdls_chr)
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
make_predrs_for_best_mdls <- function(outp_smry_ls,
                                      old_nms_chr = NULL,
                                      new_nms_chr){
  ordered_mdl_nms_chr <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Parameter=="R2") %>% dplyr::arrange(dplyr::desc(Estimate)) %>% dplyr::pull(Model)
  ind_predr_mdls_by_mdl_type_ls <- outp_smry_ls$prefd_mdl_types_chr %>% purrr::map(~{
    mdl_type_1L_chr <- .x
    paste0(outp_smry_ls$predr_cmprsn_tb$predr_chr,"_1_") %>% purrr::map(~paste0(.x,
                                                                                mdl_type_1L_chr)) %>% purrr::flatten_chr()
  })
  ordered_mdls_by_type_ls <- ind_predr_mdls_by_mdl_type_ls %>%
    purrr::map(~{
      ind_predr_mdls_chr <- .x
      ordered_mdl_nms_chr[ordered_mdl_nms_chr %>% purrr::map_lgl(~.x %in% ind_predr_mdls_chr)]
    })
  predrs_for_best_mdls_chr <- ordered_mdls_by_type_ls %>% purrr::map2_chr(outp_smry_ls$prefd_mdl_types_chr,
                                                                          ~{
                                                                            .x[1] %>% stringr::str_remove(paste0("_1_",.y))
                                                                          })
  if(!is.null(old_nms_chr)){
    predrs_for_best_mdls_chr <- transform_predr_nm_part_of_phrases(predrs_for_best_mdls_chr,
                                                                   old_nms_chr = old_nms_chr,
                                                                   new_nms_chr = new_nms_chr)
  }
  return(predrs_for_best_mdls_chr)
}
make_predr_vars_nms_ls <- function (main_predrs_chr, covars_ls)
{
    predr_vars_nms_ls <- covars_ls %>% purrr::map(~{
        covars_chr <- .x
        purrr::map(main_predrs_chr, ~list(c(.x), c(.x, covars_chr))) %>%
            purrr::flatten()
    }) %>% purrr::flatten() %>% unique()
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
make_ranked_predrs_ls <- function(descv_tbls_ls,
                                  old_nms_chr = NULL,
                                  new_nms_chr = NULL){
  unranked_predrs_chr <- rownames(descv_tbls_ls[["bl_cors_tb"]])[-1] %>%
    transform_predr_nm_part_of_phrases(old_nms_chr = old_nms_chr,
                                       new_nms_chr = new_nms_chr)
  ranks_dbl <- descv_tbls_ls[["bl_cors_tb"]][2:nrow(descv_tbls_ls[["bl_cors_tb"]]),1] %>%
    purrr::map_dbl(~{
      vec_dbl <- regmatches(.x, gregexpr("[[:digit:]]+", .x)) %>% unlist() %>% as.numeric()
      vec_dbl[1]+vec_dbl[2]/100
    }) %>% rank()
  unique_ranks_dbl <- unique(ranks_dbl) %>% sort(decreasing = T)
  ranked_predrs_chr <- purrr::map(unique_ranks_dbl, ~ unranked_predrs_chr[ranks_dbl == .x]) %>% purrr::flatten_chr()
  ranked_predrs_ls <- list(unranked_predrs_chr = unranked_predrs_chr,
                           ranked_predrs_chr = ranked_predrs_chr)
  return(ranked_predrs_ls)

}
make_shareable_mdl <- function (data_tb, mdl_smry_tb, depnt_var_nm_1L_chr = "utl_total_w",
    id_var_nm_1L_chr = "fkClientID", tfmn_1L_chr = "CLL", mdl_type_1L_chr = "OLS_CLL",
    mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NA_character_,
    seed_1L_int = 12345L)
{
    if (is.null(mdl_types_lup))
        utils::data(mdl_types_lup, envir = environment())
    all_var_nms_chr <- names(data_tb)
    tfd_depnt_var_nm_1L_chr <- all_var_nms_chr[all_var_nms_chr %>%
        startsWith(depnt_var_nm_1L_chr)]
    predr_var_nms_chr <- setdiff(all_var_nms_chr, c(tfd_depnt_var_nm_1L_chr,
        id_var_nm_1L_chr))
    if (length(predr_var_nms_chr) > 1) {
        covar_var_nms_chr <- predr_var_nms_chr[2:length(predr_var_nms_chr)]
    }else{
        covar_var_nms_chr <- NA_character_
    }
    fk_data_ls <- synthpop::syn(data_tb, visit.sequence = all_var_nms_chr[all_var_nms_chr !=
        id_var_nm_1L_chr], seed = seed_1L_int)
    model_mdl <- make_mdl(fk_data_ls$syn, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
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
make_sngl_mdl_smry_tb <- function(mdls_tb,
                                  mdl_nm_1L_chr,
                                  mdl_type_1L_chr,
                                  add_mdl_nm_sfx_1L_lgl = T){
  new_tb <- mdls_tb %>%
    dplyr::filter(Model == mdl_nm_1L_chr) %>%
    tibble::add_case(Parameter = mdl_nm_1L_chr, .before = 1) %>%
    dplyr::select(-Model)
  mdl_nm_sfx_1L_chr <- ifelse(add_mdl_nm_sfx_1L_lgl,
                              paste0("_",mdl_type_1L_chr),
                              "")
  new_tb <- new_tb %>%
    dplyr::mutate(!!rlang::sym(paste0("R2",mdl_nm_sfx_1L_chr)) := dplyr::case_when(Parameter == mdl_nm_1L_chr ~ ready4fun::get_from_lup_obj(new_tb,
                                                                                                                                            match_value_xx = "R2",
                                                                                                                                            match_var_nm_1L_chr = "Parameter",
                                                                                                                                            target_var_nm_1L_chr = "Estimate",
                                                                                                                                            evaluate_lgl = F),
                                                                                   T ~ NA_character_),
                  !!rlang::sym(paste0("Sigma",mdl_nm_sfx_1L_chr)) := dplyr::case_when(Parameter == mdl_nm_1L_chr ~ ready4fun::get_from_lup_obj(new_tb,
                                                                                                                                               match_value_xx = "Sigma",
                                                                                                                                               match_var_nm_1L_chr = "Parameter",
                                                                                                                                               target_var_nm_1L_chr = "Estimate",
                                                                                                                                               evaluate_lgl = F),
                                                                                      T ~ NA_character_)) %>%
    dplyr::filter(!Parameter %in% c("R2","RMSE","Sigma")) %>%
    dplyr::rename_with(~paste0(.x,mdl_nm_sfx_1L_chr),.cols = c("Parameter","Estimate","SE")) %>%
    dplyr::rename_with(~paste0("CI",mdl_nm_sfx_1L_chr),.cols = c("95% CI")) %>%
    dplyr::mutate(!!rlang::sym(paste0("Parameter",mdl_nm_sfx_1L_chr)) := !!rlang::sym(paste0("Parameter",mdl_nm_sfx_1L_chr)) %>% purrr::map_chr(~stringr::str_replace(.x,
                                                                                                                                                                      paste0("_1_",mdl_type_1L_chr), " model")))
  rownames(new_tb)  <- NULL
  return(new_tb)
}
make_smry_of_brm_mdl <- function (mdl_ls, data_tb, depnt_var_nm_1L_chr = "utl_total_w",
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
        dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))), quant = c(0.25,
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
make_smry_of_mdl_outp <- function (data_tb, model_mdl, folds_1L_int = 10, depnt_var_nm_1L_chr = "utl_total_w",
    start_1L_chr = NULL, tfmn_1L_chr = "NTF", predr_var_nm_1L_chr,
    covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF",
    mdl_types_lup = NULL, predn_type_1L_chr = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    data_tb <- data_tb %>% dplyr::filter(!is.na(!!rlang::sym(predr_var_nm_1L_chr)))
    data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
    mdl_desc_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
        match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
        target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F)
    folds_ls <- make_folds_ls(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        folds_1L_int = folds_1L_int)
    smry_of_one_predr_mdl_tb <- purrr::map_dfr(folds_ls, ~{
        model_mdl <- make_mdl(data_tb[-.x,], depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
            start_1L_chr = start_1L_chr, tfmn_1L_chr = tfmn_1L_chr,
            predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr,
            mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup)
        pred_old_dbl <- stats::predict(model_mdl, type = predn_type_1L_chr)
        pred_new_dbl <- stats::predict(model_mdl, newdata = data_tb[.x, ],
            type = predn_type_1L_chr) %>% calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
            tfmn_is_outp_1L_lgl = T)
        tibble::tibble(Rsquared = caret::R2(pred_old_dbl, data_tb[-.x,
            ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
            form = "traditional"), RMSE = caret::RMSE(pred_old_dbl,
            data_tb[-.x, ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
            MAE = caret::MAE(pred_old_dbl, data_tb[-.x, ] %>%
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
            RsquaredP = caret::R2(pred_new_dbl, data_tb[.x, ] %>%
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
                form = "traditional"), RMSEP = caret::RMSE(pred_new_dbl,
                data_tb[.x, ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
            MAEP = caret::MAE(pred_new_dbl, data_tb[.x, ] %>%
                dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))))
    }) %>% dplyr::summarise_all(mean) %>% dplyr::mutate(Model = mdl_desc_1L_chr) %>%
        dplyr::select(Model, dplyr::everything())
    return(smry_of_one_predr_mdl_tb)
}
make_smry_of_ts_mdl_outp <- function (data_tb, fn, predr_vars_nms_chr, mdl_nm_1L_chr, path_to_write_to_1L_chr = NA_character_,
    depnt_var_nm_1L_chr = "utl_total_w", id_var_nm_1L_chr = "fkClientID",
    round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", predictors_lup,
    backend_1L_chr = getOption("brms.backend", "rstan"), iters_1L_int = 4000L,
    seed_1L_int = 1000L)
{
  scaling_fctr_dbl <- predr_vars_nms_chr %>% purrr::map_dbl(~
      ifelse(.x %in% predictors_lup$short_name_chr,
             ready4fun::get_from_lup_obj(predictors_lup,
                                         target_var_nm_1L_chr = "mdl_scaling_dbl",
                                         match_value_xx = .x,
                                         match_var_nm_1L_chr = "short_name_chr",
                                         evaluate_lgl = F),
             1)
    )
    tfd_data_tb <- transform_tb_to_mdl_inp(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr,
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr,
        scaling_fctr_dbl = scaling_fctr_dbl)
    tfd_depnt_var_nm_1L_chr <- ifelse(identical(fn, fit_clg_log_tfmn),
        transform_depnt_var_nm_for_cll(depnt_var_nm_1L_chr), depnt_var_nm_1L_chr)
    args_ls <- list(data_tb = tfd_data_tb, depnt_var_nm_1L_chr = tfd_depnt_var_nm_1L_chr,
        predr_vars_nms_chr = predr_vars_nms_chr, iters_1L_int = iters_1L_int,
        backend_1L_chr = backend_1L_chr, seed_1L_int = seed_1L_int)
    mdl_ls <- rlang::exec(fn, !!!args_ls)
    smry_of_ts_mdl_ls <- list(smry_of_ts_mdl_tb = make_smry_of_brm_mdl(mdl_ls,
        data_tb = tfd_data_tb, depnt_var_nm_1L_chr = tfd_depnt_var_nm_1L_chr,
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
            tfd_data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
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
make_tfmn_cmprsn_plt <-  function(data_tb,
                                  depnt_var_nm_1L_chr,
                                  dictionary_tb){
  tfmn_cmprsn_plt <- tidyr::gather(data_tb %>%
                                     dplyr::mutate(!!rlang::sym(paste0(depnt_var_nm_1L_chr,"_log")) := log(!!rlang::sym(depnt_var_nm_1L_chr)),
                                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_logit")) := psych::logit(!!rlang::sym(depnt_var_nm_1L_chr)),
                                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_loglog")) := -log(-log(!!rlang::sym(depnt_var_nm_1L_chr))),
                                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_cloglog")) := log(-log(1-!!rlang::sym(depnt_var_nm_1L_chr)))), variable, value,
                                   !!rlang::sym(depnt_var_nm_1L_chr),
                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_log")),
                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_logit")),
                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_loglog")),
                                   !!rlang::sym(paste0(depnt_var_nm_1L_chr,"_cloglog"))) %>%
    dplyr::mutate(variable = factor(variable, levels = paste0(depnt_var_nm_1L_chr,c("",
                                                                                    "_log", "_logit", "_loglog", "_cloglog")),
                                    labels= c("No transformation",
                                              "Log", "Logit","Log-log", "Complementary log-log"))) %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = variable)) +
    ggalt::geom_bkde() +
    ggplot2::geom_rug() +
    viridis::scale_fill_viridis(guide = FALSE, discrete = TRUE) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::theme_bw() + ggplot2::labs(x = paste0("Transformed ",
                                                   dictionary_tb %>%
                                                     ready4fun::get_from_lup_obj(match_var_nm_1L_chr = "var_nm_chr",
                                                                                 match_value_xx = depnt_var_nm_1L_chr,
                                                                                 target_var_nm_1L_chr = "var_desc_chr",
                                                                                 evaluate_lgl = F)))
  return(tfmn_cmprsn_plt)
}
make_two_mdl_types_smry_tbl <- function(outp_smry_ls,
                                        mdls_tb){
  mdls_ls <- make_mdls_ls(outp_smry_ls,
                          mdls_tb = mdls_tb)
  two_mdl_types_smry_tbl_tb <- 1:2 %>% purrr::map_dfc(
    ~ {
      make_mdl_type_smry_tbl(mdls_tb = mdls_tb,
                             mdl_nms_chr = mdls_ls[[.x]],
                             mdl_type_1L_chr = outp_smry_ls$prefd_mdl_types_chr[.x],
                             add_mdl_nm_sfx_1L_lgl = T)
    }
  ) %>%
    dplyr::select(-paste0("Parameter_",outp_smry_ls$prefd_mdl_types_chr[2])) %>%
    dplyr::rename_with(~"Parameter",.cols = paste0("Parameter_",outp_smry_ls$prefd_mdl_types_chr[1]))
  return(two_mdl_types_smry_tbl_tb)
}
make_unique_ls_elmt_idx_int <- function (data_ls)
{
  combos_tb <- tibble::as_tibble(data_ls, .name_repair = ~paste0("r_",
                                                                 1:length(data_ls))) %>% t() %>% as.data.frame()
  combos_tb <- combos_tb %>% tibble::as_tibble()
  if(ncol(combos_tb)>1){
    combos_tb <- combos_tb %>% dplyr::mutate(V2 = dplyr::case_when(V1 == V2 ~ NA_character_,
                                                                   T ~ V2)) %>%
      dplyr::mutate(combo_chr = purrr::map2_chr(V1,
                                                V2, ~ifelse(ncol(combos_tb) == 1 | is.na(.y),
                                                            .x,
                                                            paste0(.x,"_", .y))))
    combos_tb <- combos_tb %>% dplyr::group_by(combo_chr) %>%
      dplyr::mutate(combo_id = dplyr::row_number())
    unique_ls_elmt_idx_int <- purrr::map(data_ls %>% unique(),
                                         ~ready4fun::get_from_lup_obj(combos_tb %>% dplyr::ungroup(),
                                                                      match_var_nm_1L_chr = "combo_chr", match_value_xx = paste0(.x[1],
                                                                                                                                 ifelse(is.na(.x[2]), "", paste0("_", .x[2]))),
                                                                      target_var_nm_1L_chr = "combo_id", evaluate_lgl = F)) %>%
      purrr::flatten_int()
  }else{
    unique_ls_elmt_idx_int <- 1
  }
  return(unique_ls_elmt_idx_int)
}
