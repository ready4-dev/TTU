write_box_cox_tfmn <- function(data_tb,
                               predr_var_nm_1L_chr,
                               path_to_write_to_1L_chr,
                               dep_var_nm_1L_chr = "aqol6d_total_w",
                               covar_var_nms_chr = NA_character_,
                               fl_nm_pfx_1L_chr = "A_RT",
                               height_1L_dbl = 6,
                               width_1L_dbl = 6,
                               start_1L_chr = NULL,
                               mdl_types_lup = NULL){

  if(is.null(mdl_types_lup))
    mdl_types_lup <- data("mdl_types_lup", package = "FBaqol", envir = environment())
  mdl <- make_mdl(data_tb,
                  dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                  predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                  covar_var_nms_chr = covar_var_nms_chr,
                  mdl_type_1L_chr = "OLS_NTF",
                  mdl_types_lup = mdl_types_lup,
                  start_1L_chr = start_1L_chr)
  path_to_plot_1L_chr <- write_brm_mdl_plt_fl(plt_fn = MASS::boxcox,
                                              fn_args_ls = list(mdl, plotit=T),
                                              path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                              plt_nm_1L_chr = paste0(fl_nm_pfx_1L_chr,
                                                                     "_",
                                                                     predr_var_nm_1L_chr,
                                                                     "_",
                                                                     "BOXCOX"),
                                              height_1L_dbl = height_1L_dbl,
                                              width_1L_dbl = width_1L_dbl)
  return(path_to_plot_1L_chr)

}
write_brm_model_plts <- function(mdl_ls,
                                 tfd_data_tb,
                                 mdl_nm_1L_chr,
                                 path_to_write_to_1L_chr,
                                 dep_var_nm_1L_chr = "aqol6d_total_w",
                                 dep_var_desc_1L_chr = "AQoL-6D utility score",
                                 round_var_nm_1L_chr = "round",
                                 tfmn_fn = function(x){x},
                                 units_1L_chr = "in",
                                 height_dbl = c(rep(6,2),rep(5,2)),
                                 width_dbl = c(rep(6,2),rep(6,2)),
                                 rsl_dbl = rep(300,4),
                                 args_ls = NULL,
                                 seed_1L_dbl = 23456){
  set.seed(seed_1L_dbl)
  tfd_data_tb$Predicted <- predict(mdl_ls)[,1] %>% tfmn_fn()
  plt_nms_chr <- paste0(mdl_nm_1L_chr, "_", c("coefs","hetg", "dnst","sctr_plt"))
  mdl_plts_paths_ls <- purrr::map(1:4,
                                  ~  {
                                    plt_fn <- fn_args_ls <- NULL
                                    if(.x %in% c(1,2)){
                                      plt <- plot(mdl_ls, ask=F, plot = F)
                                      if(length(plt)>=.x){
                                        fn_args_ls <- list(mdl_ls = mdl_ls,
                                                           idx_1L_int = as.integer(.x))
                                        plt_fn <- function(mdl_ls, idx_1L_int){plot(mdl_ls, ask=F, plot = F)[idx_1L_int]}
                                      }
                                    }else{
                                      if(.x == 3){
                                        plt_fn <- plot_obsd_predd_dnst
                                        fn_args_ls <- list(tfd_data_tb = tfd_data_tb)
                                      }else{
                                        plt_fn <- plot_obsd_predd_sctr
                                        fn_args_ls <- list(tfd_data_tb = tfd_data_tb,
                                                           dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                                                           dep_var_desc_1L_chr = dep_var_desc_1L_chr,
                                                           round_var_nm_1L_chr = round_var_nm_1L_chr,
                                                           #plot_fn = sctr_plt_fn,
                                                           args_ls = args_ls)
                                      }
                                    }
                                    write_brm_mdl_plt_fl(plt_fn,
                                                         fn_args_ls = fn_args_ls,
                                                         path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                                         plt_nm_1L_chr = plt_nms_chr[.x],
                                                         units_1L_chr = units_1L_chr,
                                                         width_1L_dbl =  width_dbl[.x],
                                                         height_1L_dbl = height_dbl[.x],
                                                         rsl_1L_dbl = rsl_dbl[.x])
                                    }) %>%
    stats::setNames(plt_nms_chr) %>%
    purrr::discard(is.na)
  return(mdl_plts_paths_ls)
}
write_brm_mdl_plt_fl <- function(plt_fn = NULL,
                                 fn_args_ls = NULL,
                                 path_to_write_to_1L_chr,
                                 plt_nm_1L_chr,
                                 grpx_fn = grDevices::png,
                                 units_1L_chr = "in",
                                 width_1L_dbl = 6,
                                 height_1L_dbl = 6,
                                 rsl_1L_dbl = 300){
  if(!is.null(plt_fn)){
    path_to_plot_1L_chr <- paste0(path_to_write_to_1L_chr,
                                  "/",
                                  plt_nm_1L_chr,
                                  ifelse(identical(grpx_fn,grDevices::png),".png",".tiff"))
    rlang::exec(grpx_fn,
                !!!list(path_to_plot_1L_chr,
                        units = units_1L_chr,
                        width = width_1L_dbl,
                        height = height_1L_dbl,
                        res = rsl_1L_dbl))
    plt <- rlang::exec(plt_fn,!!!fn_args_ls)
    print(plt)
    dev.off()
  }else{
    path_to_plot_1L_chr <- NA_character_
  }
  return(path_to_plot_1L_chr)
}
write_mdl_plts <- function(data_tb,
                           mdl,
                           mdl_fl_nm_1L_chr = "OLS_NTF",
                           dep_var_nm_1L_chr = "aqol6d_total_w",
                           tfmn_1L_chr = "NTF",
                           predr_var_nm_1L_chr,
                           predr_var_desc_1L_chr,
                           predr_vals_dbl,
                           covar_var_nms_chr = NA_character_,
                           path_to_write_to_1L_chr,
                           pred_type_1L_chr = NULL,
                           tfmn_for_bnml_1L_lgl = F,
                           family_1L_chr = NA_character_,
                           plt_idcs_int = 1:5){
  data_tb <- transform_ds_for_mdlng(data_tb,
                                    dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                                    predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                    covar_var_nms_chr = covar_var_nms_chr)
  tfd_data_tb <- transform_data_tb_for_cmprsn(data_tb,
                                              mdl = mdl,
                                              pred_type_1L_chr = pred_type_1L_chr,
                                              tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl,
                                              family_1L_chr = family_1L_chr)
  if(1 %in% plt_idcs_int){
    predn_ds_tb <- make_predn_ds_with_one_predr(mdl,
                                                dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                                                tfmn_1L_chr = tfmn_1L_chr,
                                                predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                                predr_vals_dbl = predr_vals_dbl,
                                                pred_type_1L_chr = pred_type_1L_chr)
  }else{
    predn_ds_tb <- NULL
  }
  purrr::pwalk(list(plt_fn_ls = list(plot_lnr_cmprsn, plot_auto_lm, plot_obsd_predd_dnst, plot_obsd_predd_dnst, plot_lnr_cmprsn_sctr_plt)[plt_idcs_int],
                    fn_args_ls_ls = list(list(data_tb = data_tb,
                                              predn_ds_tb = predn_ds_tb,
                                              predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                              predr_var_desc_1L_chr = predr_var_desc_1L_chr),
                                         list(mdl,
                                              which_dbl = 1:6,
                                              ncol_1L_int = 3L,
                                              label_size_1L_int = 3),
                                         list(tfd_data_tb = tfd_data_tb),
                                         list(tfd_data_tb = transform_data_tb_for_cmprsn(data_tb,
                                                                                         mdl = mdl,
                                                                                         tf_type_1L_chr = ifelse(!4 %in% plt_idcs_int,"Predicted","Simulated"),
                                                                                         pred_type_1L_chr = NULL,
                                                                                         tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl,
                                                                                         family_1L_chr = family_1L_chr),
                                              predd_val_var_nm_1L_chr = "Simulated"),
                                         list(tfd_data_tb = tfd_data_tb))[plt_idcs_int],
                    plt_nm_sfx_chr = c("_LNR_CMPRSN","_AUTOPLT","_PRED_DNSTY","_SIM_DNSTY","_PRED_SCTR")[plt_idcs_int],
                    size_ls = list(c(6,6),c(4,7),c(6,6),c(6,6),c(6,6))[plt_idcs_int]),
               ~ write_brm_mdl_plt_fl(plt_fn = ..1,
                                      fn_args_ls = ..2,
                                      path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                      plt_nm_1L_chr = paste0(mdl_fl_nm_1L_chr,
                                                             ifelse(!is.na(covar_var_nms_chr[1]),
                                                                    paste("_",paste0(covar_var_nms_chr[1:min(length(covar_var_nms_chr),3)], collapse="")),
                                                                    ""),
                                                             ..3),
                                      height_1L_dbl = ..4[1],
                                      width_1L_dbl = ..4[2]))

}
write_mdl_type_covars_mdls <- function(data_tb,
                                       predrs_var_nms_chr,
                                       covar_var_nms_chr,
                                       mdl_type_1L_chr,
                                       path_to_write_to_1L_chr,
                                       fl_nm_pfx_1L_chr = "D_CT",
                                       mdl_types_lup = NULL,
                                       start_1L_chr = NA_character_){
  if(is.null(mdl_types_lup))
    mdl_types_lup <- data("mdl_types_lup", package = "FBaqol", envir = environment())
  smry_of_mdls_with_covars_tb <- purrr::map_dfr(predrs_var_nms_chr,
                                                ~{
                                                  mdl <- make_mdl(data_tb,
                                                                  predr_var_nm_1L_chr = .x,
                                                                  covar_var_nms_chr = covar_var_nms_chr,
                                                                  tfmn_1L_chr = "NTF",
                                                                  mdl_type_1L_chr = mdl_type_1L_chr,
                                                                  control_1L_chr = NA_character_,
                                                                  mdl_types_lup = mdl_types_lup,
                                                                  start_1L_chr = start_1L_chr)
                                                  mdl_fl_nm_1L_chr <- paste0(fl_nm_pfx_1L_chr,"_",.x,
                                                                             "_",mdl_type_1L_chr)
                                                  saveRDS(mdl,paste0(path_to_write_to_1L_chr,"/",mdl_fl_nm_1L_chr,".RDS"))

                                                  tibble::tibble(variable = .x,
                                                                 Rsquare=caret::R2(data_tb %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr)), predict(mdl), form = "traditional"),
                                                                 AIC=AIC(mdl), BIC=BIC(mdl),
                                                                 Significant=paste(names(which(summary(mdl)$coefficients[,4]<0.01)),
                                                                                   collapse = " "))
                                                })
  smry_of_mdls_with_covars_tb <- smry_of_mdls_with_covars_tb %>%
    dplyr::arrange(dplyr::desc(AIC))
  saveRDS(smry_of_mdls_with_covars_tb,paste0(path_to_write_to_1L_chr,"/",
                                             paste0(fl_nm_pfx_1L_chr,"_","SMRY",
                                                    "_",mdl_type_1L_chr),
                                             ".RDS"))
  return(smry_of_mdls_with_covars_tb)

}
write_mdl_type_multi_outps <- function(data_tb,
                                       n_folds_1L_int = 10,
                                       predrs_var_nms_chr,
                                       covar_var_nms_chr = NA_character_,
                                       start_1L_chr = NULL,
                                       mdl_type_1L_chr,
                                       dep_var_nm_1L_chr = "aqol6d_total_w",
                                       path_to_write_to_1L_chr,
                                       mdl_types_lup = NULL,
                                       fl_nm_pfx_1L_chr = "C_PREDR",
                                       plt_idcs_int = c(3,5)){
  if(is.null(mdl_types_lup))
    mdl_types_lup <- data("mdl_types_lup", package = "FBaqol", envir = environment())
  smry_of_mdl_sngl_predrs_tb <- purrr::map_dfr(predrs_var_nms_chr,
                                               ~ {
                                                 tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
                                                                                            match_var_nm_1L_chr = "short_name_chr",
                                                                                            match_value_xx = mdl_type_1L_chr,
                                                                                            target_var_nm_1L_chr = "tfmn_chr",
                                                                                            evaluate_lgl = F)
                                                 mdl_smry_tb <- write_mdl_type_sngl_outps(data_tb,
                                                                                          n_folds_1L_int = n_folds_1L_int,
                                                                                          dep_var_nm_1L_chr  = dep_var_nm_1L_chr ,
                                                                                          tfmn_1L_chr = tfmn_1L_chr,
                                                                                          start_1L_chr = start_1L_chr,
                                                                                          predr_var_nm_1L_chr = .x,
                                                                                          predr_var_desc_1L_chr = NA_character_,
                                                                                          predr_vals_dbl = NA_real_,
                                                                                          covar_var_nms_chr = covar_var_nms_chr,
                                                                                          mdl_type_1L_chr = mdl_type_1L_chr,
                                                                                          path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                                                                          mdl_types_lup = mdl_types_lup,
                                                                                          mdl_fl_nm_1L_chr = paste0(fl_nm_pfx_1L_chr,
                                                                                                                    "_",
                                                                                                                    .x,
                                                                                                                    "_",
                                                                                                                    mdl_type_1L_chr#"_MDL"
                                                                                          ),
                                                                                          plt_idcs_int = plt_idcs_int)
                                                 if(!is.null(n_folds_1L_int)){
                                                   mdl_smry_tb <- mdl_smry_tb %>%
                                                     dplyr::select((-Model)) %>%
                                                     dplyr::mutate(Predictor = .x) %>%
                                                     dplyr::select(Predictor, dplyr::everything())
                                                 }
                                                 mdl_smry_tb
                                               })
  if(!is.null(n_folds_1L_int)){
    smry_of_mdl_sngl_predrs_tb <- smry_of_mdl_sngl_predrs_tb %>%
      dplyr::arrange(dplyr::desc(RsquaredP))
  }
  return(smry_of_mdl_sngl_predrs_tb)

}
write_mdl_type_sngl_outps <- function(data_tb,
                                      n_folds_1L_int = 10,
                                      dep_var_nm_1L_chr = "aqol6d_total_w",
                                      start_1L_chr = NULL,
                                      tfmn_1L_chr = "NTF",
                                      predr_var_nm_1L_chr,
                                      predr_var_desc_1L_chr,
                                      predr_vals_dbl,
                                      covar_var_nms_chr = NA_character_,
                                      mdl_type_1L_chr = "OLS_NTF",
                                      mdl_types_lup = NULL,
                                      path_to_write_to_1L_chr,
                                      mdl_fl_nm_1L_chr,
                                      plt_idcs_int = NA_integer_){
  if(is.null(mdl_types_lup))
    mdl_types_lup <- data("mdl_types_lup", package = "FBaqol", envir = environment())
  arg_vals_chr <- c("control_chr","family_chr","pred_type_chr") %>%
    purrr::map_chr(~ready4fun::get_from_lup_obj(mdl_types_lup,
                                                match_var_nm_1L_chr = "short_name_chr",
                                                match_value_xx = mdl_type_1L_chr,
                                                target_var_nm_1L_chr = .x,
                                                evaluate_lgl = F))
  control_1L_chr <- arg_vals_chr[1]
  family_1L_chr <- arg_vals_chr[2]
  pred_type_1L_chr <- arg_vals_chr[3]
  if(is.na(pred_type_1L_chr))
    pred_type_1L_chr <- NULL
  if(is.na(plt_idcs_int[1])){
    plt_idcs_int <- 1:5
    if(!is.na(control_1L_chr)){
      if(control_1L_chr %>% startsWith("betareg"))
        plt_idcs_int <- c(1,3,5)
    }
  }
  tfmn_for_bnml_1L_lgl <- ready4fun::get_from_lup_obj(mdl_types_lup,
                                                      match_var_nm_1L_chr = "short_name_chr",
                                                      match_value_xx = mdl_type_1L_chr,
                                                      target_var_nm_1L_chr = "tfmn_for_bnml_lgl",
                                                      evaluate_lgl = F)
  data_tb <- data_tb %>% dplyr::mutate(!!rlang::sym(transform_dep_var_nm(dep_var_nm_1L_chr,
                                                                         tfmn_1L_chr = tfmn_1L_chr)) := !!rlang::sym(dep_var_nm_1L_chr) %>% calculate_dep_var_tfmn())
  mdl <- make_mdl(data_tb,
                  dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                  tfmn_1L_chr = tfmn_1L_chr,
                  predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                  covar_var_nms_chr = covar_var_nms_chr,
                  mdl_type_1L_chr = mdl_type_1L_chr,
                  mdl_types_lup = mdl_types_lup,
                  control_1L_chr = control_1L_chr,
                  start_1L_chr = start_1L_chr)
  write_mdl_plts(data_tb,
                 mdl = mdl,
                 mdl_fl_nm_1L_chr = mdl_fl_nm_1L_chr,
                 dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                 tfmn_1L_chr = tfmn_1L_chr,
                 predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                 predr_var_desc_1L_chr = predr_var_desc_1L_chr,
                 predr_vals_dbl = predr_vals_dbl,
                 covar_var_nms_chr = covar_var_nms_chr,
                 path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                 pred_type_1L_chr = pred_type_1L_chr,
                 tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl,
                 family_1L_chr = family_1L_chr,
                 plt_idcs_int = plt_idcs_int)
  if(!is.null(n_folds_1L_int)) {
    smry_of_one_predr_mdl_tb <- make_smry_of_mdl(data_tb,
                                                 mdl = mdl,
                                                 n_folds_1L_int = n_folds_1L_int,
                                                 dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                                                 tfmn_1L_chr = tfmn_1L_chr,
                                                 predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                                 mdl_type_1L_chr = mdl_type_1L_chr,
                                                 mdl_types_lup = mdl_types_lup,
                                                 pred_type_1L_chr = pred_type_1L_chr)
  }else{
    smry_of_one_predr_mdl_tb <- tibble::tibble()
  }
  saveRDS(mdl,paste0(path_to_write_to_1L_chr,"/",mdl_fl_nm_1L_chr,".RDS"))
  return(smry_of_one_predr_mdl_tb)
}
write_predr_cmprsn_outps <- function(data_tb,
                                     path_to_write_to_1L_chr,
                                     dep_var_nm_1L_chr = "aqol6d_total_w",
                                     candidate_predrs_chr,
                                     max_nbr_of_boruta_mdl_runs_int = 300L){
  if(length(candidate_predrs_chr)>1){
    covar_var_nms_chr <- candidate_predrs_chr[2:length(candidate_predrs_chr)]
  }else{
    covar_var_nms_chr <- NA_character_
  }
  data_tb <- transform_ds_for_mdlng(data_tb,
                                    dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                                    predr_var_nm_1L_chr = candidate_predrs_chr[1],
                                    covar_var_nms_chr = covar_var_nms_chr)
  rf_mdl <- randomForest::randomForest(as.formula(paste0(dep_var_nm_1L_chr," ~ .")), data=data_tb, importance = TRUE)
  boruta_mdl <- Boruta::Boruta(as.formula(paste0(dep_var_nm_1L_chr," ~ .")), data = data_tb, maxRuns = max_nbr_of_boruta_mdl_runs_int)
  purrr::pwalk(list(fn_ls = list(randomForest::varImpPlot,plot),
                    fn_args_ls_ls = list(list(rf_mdl,
                                              main = ""),
                                         list(boruta_mdl,
                                              cex=1.5,
                                              cex.axis=0.8,#pin=c(6,4),
                                              las=2,
                                              xlab="",
                                              main="")),
                    plt_nm_sfx_chr = c("_RF_VAR_IMP","_BORUTA_VAR_IMP"),
                    size_ls = list(c(6,6),c(4,6))),
               ~ write_brm_mdl_plt_fl(plt_fn = ..1,
                                      fn_args_ls = ..2,
                                      path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                      plt_nm_1L_chr = paste0("B_PRED_CMPRSN",
                                                             ..3),
                                      height_1L_dbl = ..4[1],
                                      width_1L_dbl = ..4[2]))
  confirmed_predrs_chr <- names(boruta_mdl$finalDecision)[boruta_mdl$finalDecision=="Confirmed"]
  confirmed_predrs_tb <- rf_mdl$importance %>% tibble::as_tibble(rownames = "predr_chr") %>%
    dplyr::arrange(dplyr::desc(`%IncMSE`)) %>%
    dplyr::filter(predr_chr %in% confirmed_predrs_chr)
  return(confirmed_predrs_tb)
}
write_rndrd_rprt <- function(path_to_RMD_dir_1L_chr,
                             nm_of_RMD_1L_chr = "report.RMD",
                             params_ls = list(output_type_1L_chr = "HTML"),
                             rltv_path_to_outpt_yaml_1L_chr = "output_yml", # must be relative to path_to_RMD_dir_1L_chr
                             paths_to_fls_to_copy_chr = NA_character_, # main RMD must be at top level
                             path_to_write_fls_to_1L_chr = NA_character_,
                             nm_of_rprt_dir_1L_chr = "Reports",
                             path_to_outpt_rtrp_1L_chr = "./",
                             file_nm_1L_chr,
                             overwrite_1L_lgl = T){
  if(!is.na(path_to_write_fls_to_1L_chr)){
    path_to_rprt_dir_1L_chr <- paste0(path_to_write_fls_to_1L_chr,"/",nm_of_rprt_dir_1L_chr)
    if(!dir.exists(path_to_rprt_dir_1L_chr))
      dir.create(path_to_rprt_dir_1L_chr)
    if(is.na(paths_to_fls_to_copy_chr[1]))
      paths_to_fls_to_copy_chr <- list.files(path_to_RMD_dir_1L_chr, full.names = T)
    file.copy(paths_to_fls_to_copy_chr, path_to_rprt_dir_1L_chr, overwrite = overwrite_1L_lgl)
    path_to_wd_1L_chr <- path_to_rprt_dir_1L_chr
  }else{
    path_to_wd_1L_chr <- path_to_RMD_dir_1L_chr
  }
  if(!dir.exists(path_to_outpt_rtrp_1L_chr))
    dir.create(path_to_outpt_rtrp_1L_chr)
  path_to_RMD_1L_chr <- paste0(path_to_wd_1L_chr,"/",nm_of_RMD_1L_chr)
  rmarkdown::render(path_to_RMD_1L_chr,
                    switch(params_ls$output_type_1L_chr,
                           PDF = "bookdown::pdf_book",
                           HTML = "bookdown::html_document2",
                           Word = "officedown::rdocx_document"),
                    output_yaml = paste0(path_to_wd_1L_chr,"/",rltv_path_to_outpt_yaml_1L_chr),
                    params = params_ls,
                    envir = new.env(),
                    output_file = paste0(file_nm_1L_chr,
                      ".",
                      ifelse(params_ls$output_type_1L_chr=="Word",
                             "docx",
                             tolower(params_ls$output_type_1L_chr))),
                    output_dir = path_to_outpt_rtrp_1L_chr)
}
write_results_to_csv <- function(synth_data_spine_ls,
                              output_dir_1L_chr = "."){
  measurements_tb <- tibble::tibble(timepoint_nms_chr = synth_data_spine_ls$timepoint_nms_chr,
                                    nbr_obs_dbl = synth_data_spine_ls$nbr_obs_dbl)
  var_summ_res_tb <- suppressMessages(purrr::map_dfr(1:length(synth_data_spine_ls$timepoint_nms_chr),
                                    ~ {
                                      idx_dbl <- .x
                                      suppressWarnings({synth_data_spine_ls[c(4:6)]  %>%
                                          purrr::map_dfc(~.x[idx_dbl])}) %>%
                                        stats::setNames(c("Mean","SD","N_Missing")) %>%
                                        dplyr::mutate(var_names_chr = synth_data_spine_ls$var_names_chr,
                                                      timepoint_nms_chr = synth_data_spine_ls$timepoint_nms_chr[idx_dbl]) %>%
                                        dplyr::select(timepoint_nms_chr,var_names_chr,dplyr::everything())
                                    }))
  cor_tb_ls <- synth_data_spine_ls$cor_mat_ls %>% purrr::map(~tibble::as_tibble(.x) %>% stats::setNames(synth_data_spine_ls$var_names_chr) %>% dplyr::mutate(var_names_chr = synth_data_spine_ls$var_names_chr) %>%
                                                                 dplyr::select(var_names_chr,dplyr::everything())) %>%
    stats::setNames(paste0(synth_data_spine_ls$timepoint_nms_chr,"_correlations_tb"))
  var_class_pars_tb <- synth_data_spine_ls[7:9] %>% tibble::as_tibble() %>% dplyr::mutate(min_dbl = purrr::map_dbl(min_max_ls,~.x[1]),
                                                                                          max_dbl = purrr::map_dbl(min_max_ls,~.x[2])) %>%
    dplyr::select(var_names_chr,dplyr::everything(),-min_max_ls)
  output_ls <- list(measurements_tb = measurements_tb,
                    var_summ_res_tb = var_summ_res_tb,
                    var_class_pars_tb = var_class_pars_tb) %>% append(cor_tb_ls)
  dss_tb <- tibble::tibble(ds_obj_nm_chr = names(output_ls),
                           title_chr = c("Brief summary table of the number of observations for which data was collected at each study timepoint.",
                                         "Summary statistics (Mean, SD and Number Missing) for AQoL6D health utility and six mental health outcome measures for each study timepoint.",
                                         "Brief information about the data structure (whether discrete and allowable range) of AQoL6D health utility and six mental health outcome variables.",
                                         paste0("Correlation matrix for AQoL6D health utility and six mental health outcome measures at the ",synth_data_spine_ls$timepoint_nms_chr," study timepoint.")))
  purrr::walk2(output_ls, names(output_ls), ~ write.csv(.x, file = paste0(output_dir_1L_chr,"/",.y,".csv"),
                                                        row.names = F))
  return(dss_tb)
}
write_sngl_predr_multi_mdls_outps <- function(data_tb,
                                              mdl_types_chr,
                                              predr_var_nm_1L_chr,
                                              predr_var_desc_1L_chr,
                                              predr_vals_dbl,
                                              path_to_write_to_1L_chr,
                                              start_1L_chr = NULL,
                                              covar_var_nms_chr = NA_character_,
                                              dep_var_nm_1L_chr = "aqol6d_total_w",
                                              n_folds_1L_int = 10,
                                              mdl_types_lup = NULL,
                                              fl_nm_pfx_1L_chr = "A_RT_",
                                              plt_idcs_int = NA_integer_){
  if(is.null(mdl_types_lup))
    mdl_types_lup <- data("mdl_types_lup", package = "FBaqol", envir = environment())
  data_tb <- transform_ds_for_mdlng(data_tb,
                                    dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                                    predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                    covar_var_nms_chr = covar_var_nms_chr)
  smry_of_sngl_predr_mdls_tb <- purrr::map_dfr(mdl_types_chr,
                                               ~ {
                                                 tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdl_types_lup,
                                                                                            match_var_nm_1L_chr = "short_name_chr",
                                                                                            match_value_xx = .x,
                                                                                            target_var_nm_1L_chr = "tfmn_chr",
                                                                                            evaluate_lgl = F)
                                                 write_mdl_type_sngl_outps(data_tb,
                                                                           n_folds_1L_int = n_folds_1L_int,
                                                                           dep_var_nm_1L_chr  = dep_var_nm_1L_chr,
                                                                           start_1L_chr = start_1L_chr,
                                                                           tfmn_1L_chr = tfmn_1L_chr,
                                                                           predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                                                           predr_var_desc_1L_chr = predr_var_desc_1L_chr,
                                                                           predr_vals_dbl = predr_vals_dbl,
                                                                           covar_var_nms_chr = covar_var_nms_chr,
                                                                           mdl_type_1L_chr = .x,
                                                                           path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                                                           mdl_types_lup = mdl_types_lup,
                                                                           mdl_fl_nm_1L_chr = paste0(fl_nm_pfx_1L_chr,
                                                                                                     predr_var_nm_1L_chr,
                                                                                                     "_",
                                                                                                     .x#,"_MDL"
                                                                           ),
                                                                           plt_idcs_int = plt_idcs_int)
                                               })
  if(!is.null(n_folds_1L_int))
    smry_of_sngl_predr_mdls_tb <- smry_of_sngl_predr_mdls_tb %>%
    dplyr::arrange(dplyr::desc(RsquaredP))
  return(smry_of_sngl_predr_mdls_tb)

}
write_ts_mdls <- function(data_tb,
                          dep_var_nm_1L_chr = "aqol6d_total_w",
                          predr_vars_nms_ls,
                          id_var_nm_1L_chr = "fkClientID",
                          round_var_nm_1L_chr = "round",
                          round_bl_val_1L_chr = "Baseline",
                          fn_ls,
                          mdl_nms_ls,
                          mdl_smry_dir_1L_chr,
                          iters_1L_int = 4000L,
                          seed_1L_int = 1000L){
  if(!dir.exists(mdl_smry_dir_1L_chr))
    dir.create(mdl_smry_dir_1L_chr)
  mdls_smry_tb <- purrr::map_dfr(1:length(mdl_nms_ls),
                                 ~ {
                                   idx_1L_int <- .x
                                   purrr::map2_dfr(fn_ls,
                                                   mdl_nms_ls[[idx_1L_int]],
                                                   ~ {
                                                     smry_ls <- make_smry_of_ts_mdl(data_tb = data_tb,
                                                                                    fn = .x,
                                                                                    predr_vars_nms_chr = predr_vars_nms_ls[[idx_1L_int]],
                                                                                    mdl_nm_1L_chr = .y,
                                                                                    path_to_write_to_1L_chr = mdl_smry_dir_1L_chr,
                                                                                    dep_var_nm_1L_chr = dep_var_nm_1L_chr,
                                                                                    id_var_nm_1L_chr = id_var_nm_1L_chr,
                                                                                    round_var_nm_1L_chr = round_var_nm_1L_chr,
                                                                                    round_bl_val_1L_chr = round_bl_val_1L_chr,
                                                                                    iters_1L_int = iters_1L_int,
                                                                                    seed_1L_int = seed_1L_int)
                                                     Sys.sleep(5)
                                                     smry_ls$smry_of_ts_mdl_tb
                                                   })
                                 })
  saveRDS(mdls_smry_tb, paste0(mdl_smry_dir_1L_chr,"/mdls_smry_tb.RDS"))
  return(mdls_smry_tb)
}
