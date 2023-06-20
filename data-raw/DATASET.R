library(ready4fun)
library(ready4use)
library(ready4show)
library(youthvars)
library(scorz)
library(specific)
ready4fun::write_fn_type_dirs()
# MANUAL STEP. Write all your functions to R files in the new "fns" directory.
fns_env_ls <- ready4fun::read_fns(c("data-raw/fns/","data-raw/mthds/"),
                                  fns_env = new.env(parent = globalenv()))
x <- ready4fun::make_pkg_desc_ls(pkg_title_1L_chr = "Implement Transfer to Utility Mapping Algorithms With Ready4",
                                                    pkg_desc_1L_chr = "Tools for developing, reporting and sharing utility mapping algorithms for use with the ready4 youth mental health systems model (https://ready4-dev.github.io/ready4/).
                            This development version of the TTU package has been made available as part of the process of testing and documenting the package.
                            If you have any questions, please contact the authors (matthew.hamilton@orygen.org.au).",
                                                    authors_prsn = c(utils::person(given = "Caroline",family = "Gao",email = "caroline.gao@orygen.org.au", role = c("aut"),comment = c(ORCID = "0000-0002-0987-2759")),
                                                                     utils::person(given = "Matthew",family = "Hamilton",email = "matthew.hamilton@orygen.org.au", role = c("aut", "cre"),comment = c(ORCID = "0000-0001-7407-9194")),
                                                                     utils::person("Orygen", role = c("cph", "fnd")),
                                                                     utils::person("Headspace", role = c( "fnd")),
                                                                     utils::person("National Health and Medical Research Council", role = c( "fnd"))),
                                                    urls_chr = c("https://ready4-dev.github.io/TTU/",
                                                                 "https://github.com/ready4-dev/TTU",
                                                                 "https://ready4-dev.github.io/ready4/")) %>%
  ready4fun::make_manifest(addl_pkgs_ls = ready4fun::make_addl_pkgs_ls(suggests_chr = c("betareg","caret","knitr","knitrBootstrap","rmarkdown"),
                                                                       #imports_chr = c(),
                                                                       depends_chr = "specific"),
                           build_ignore_ls = ready4fun::make_build_ignore_ls(file_nms_chr = c("initial_setup.R")),
                           check_type_1L_chr = "ready4",
                           copyright_holders_chr = "Orygen",
                           custom_dmt_ls = ready4fun::make_custom_dmt_ls(),##
                           dev_pkgs_chr = c("cmdstanr",
                                            "ready4",#"ready4fun",
                                            "ready4use","ready4show",
                                            #"youthvars","scorz",
                                            "specific"),
                           lifecycle_stage_1L_chr = "experimental",
                           path_to_pkg_logo_1L_chr = "../../../../../Documentation/Images/TTU-logo/default.png",
                           piggyback_to_1L_chr = "ready4-dev/ready4",
                           ready4_type_1L_chr = "modelling",
                           zenodo_badge_1L_chr = "[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5646593.svg)](https://doi.org/10.5281/zenodo.5646593)")
y <- ready4class::ready4class_constructor() %>%
  dplyr::bind_rows(ready4class::make_pt_ready4class_constructor(make_s3_lgl = FALSE,
                                                                name_stub_chr = "Synopsis",
                                                                slots_ls = list("a_Ready4showPaths",
                                                                                "b_SpecificResults",
                                                                                "c_SpecificParameters",
                                                                                "d_YouthvarsProfile",
                                                                                "e_Ready4useRepos") %>% list(),
                                                                pt_ls = list("Ready4showPaths",
                                                                             "SpecificResults",
                                                                             "SpecificParameters",
                                                                             "YouthvarsProfile",
                                                                             "Ready4useRepos") %>% list(),
                                                                class_desc_chr = "Input, Output and Authorship Data For Generating Utility Mapping Study Reports.",
                                                                parent_class_chr = "SpecificSynopsis"),
                   ready4class::make_pt_ready4class_constructor(make_s3_lgl = FALSE,
                                                                name_stub_chr = "Reports",
                                                                slots_ls = list("a_TTUSynopsis",
                                                                                "catalogue_tmpl_chr",
                                                                                "catalogue_fl_nms_ls",
                                                                                "manuscript_tmpl_chr",
                                                                                "manuscript_fl_nms_ls") %>% list(),
                                                                pt_ls = list("TTUSynopsis",
                                                                             "character",
                                                                             "list",
                                                                             "character",
                                                                             "list") %>% list(),
                                                                vals_ls = list(list(catalogue_tmpl_chr = "c(\"https://github.com/ready4-dev/ttu_mdl_ctlg\",\"0.0.9.7\")", #### UPDATE ###
                                                                                    catalogue_fl_nms_ls = "ready4show::make_rmd_fl_nms_ls(\"Lngl_Mdls_HTML\",
                                                                                                                           pdf_fl_nm_1L_chr = \"Lngl_Mdls_PDF\",
                                                                                                                           word_fl_nm_1L_chr = \"Lngl_Mdls_Word\")",
                                                                                    manuscript_tmpl_chr = "c(\"https://github.com/ready4-dev/ttu_lng_ss\",\"0.8.0.0\")",
                                                                                    manuscript_fl_nms_ls = "ready4show::make_rmd_fl_nms_ls(pdf_fl_nm_1L_chr = \"Main_PDF\",
                                                                                                                            word_fl_nm_1L_chr = \"Main_Word\")")),
                                                                class_desc_chr = "Metadata to produce utility mapping study reports.",
                                                                parent_class_chr = "Ready4Module",
                                                                inc_clss_ls = list("TTUSynopsis") %>% list()),
                     ready4class::make_pt_ready4class_constructor(make_s3_lgl = FALSE,
                                                                  name_stub_chr = "Project",
                                                                  slots_ls = list("a_ScorzProfile",
                                                                                  "b_SpecificParameters",
                                                                                  "c_SpecificProject",
                                                                                  "d_TTUReports") %>% list(),
                                                                  pt_ls = list("ScorzProfile",
                                                                               "SpecificParameters",
                                                                               "SpecificProject",
                                                                               "TTUReports") %>% list(),
                                                                  class_desc_chr = "Input And Output Data For Undertaking and Reporting Utility Mapping Studies.",
                                                                  parent_class_chr = "Ready4Module",
                                                                  inc_clss_ls = list("TTUReports") %>% list()))
datasets_ls <- list(tibble::tibble(short_name_chr = c("OLS_NTF",
                                                      "OLS_LOG",
                                                      "OLS_LOGIT",
                                                      "OLS_LOGLOG",
                                                      "OLS_CLL",
                                                      "GLM_GSN_LOG",
                                                      "GLM_BNL_LOG",
                                                      "GLM_BNL_LGT",
                                                      "GLM_BNL_CLL",
                                                      "BET_LOG",
                                                      "BET_LGT",
                                                      "BET_CLL"),
                                   long_name_chr = c("Ordinary Least Squares (no transformation)",
                                                     "Ordinary Least Squares (log transformation)",
                                                     "Ordinary Least Squares (logit transformation)",
                                                     "Ordinary Least Squares (log log transformation)",
                                                     "Ordinary Least Squares (complementary log log transformation)",
                                                     "Generalised Linear Model with Gaussian distribution and log link",
                                                     "Generalised Linear Model with Binomial distribution and log link",
                                                     "Generalised Linear Model with Binomial distribution and logit link",
                                                     "Generalised Linear Model with Binomial distribution and complementary log log link",
                                                     "Beta Regression Model with Binomial distribution and log link",
                                                     "Beta Regression Model with Binomial distribution and logit link",
                                                     "Beta Regression Model with Binomial distribution and complementary log log link")) %>%
                      dplyr::mutate(control_chr = c(rep(NA_character_,9),rep("betareg::betareg.control",3)),
                                    family_chr = c(rep(NA_character_,5),"gaussian(log)","quasibinomial(log)","quasibinomial(logit)","quasibinomial(cloglog)", rep(NA_character_,3)),
                                    fn_chr = c(rep("lm",5),rep("glm",4),rep("betareg::betareg",3)),
                                    start_chr = c(rep(NA_character_,5),
                                                  rep("-0.1,-0.1",4),
                                                  rep("-0.5,-0.1,3",3)),

                                    predn_type_chr = c(rep(NA_character_,5),
                                                       rep("response",7)),
                                    tfmn_chr = dplyr::case_when(startsWith(short_name_chr, "OLS_") ~ purrr::map_chr(short_name_chr, ~ {
                                      idx_1L_int <- 1 + stringi::stri_locate_last_fixed(.x,"_")[1,1] %>% as.vector()
                                      stringr::str_sub(.x, start = idx_1L_int)
                                    }),
                                    T ~ "NTF"),
                                    tfmn_for_bnml_lgl = c(rep(F,6),rep(T,3),rep(F,3)),
                                    fixed_acronym_chr = dplyr::case_when(startsWith(short_name_chr, "OLS_") ~ "OLS",
                                                                         T ~ "GLM"),
                                    mixed_acronym_chr = dplyr::case_when(startsWith(short_name_chr, "OLS_") ~ "LMM",
                                                                         T ~ "GLMM"),
                                    mixed_type_chr = dplyr::case_when(startsWith(short_name_chr, "OLS_") ~ "linear mixed model",
                                                                      T ~ "generalised linear mixed model"),
                                    with_chr = long_name_chr %>%
                                      purrr::map_chr(~ifelse(startsWith(.x,"Ordinary Least Squares"),
                                                             stringr::str_remove(.x,"Ordinary Least Squares ") %>%
                                                               stringr::str_sub(start = 2, end = -2),
                                                             stringr::str_remove(.x,"Generalised Linear Model with ") %>%
                                                               stringr::str_remove("Beta Regression Model with ")))) %>%
                      ready4fun::make_pkg_ds_ls(db_1L_chr = "mdl_types_lup",
                                                title_1L_chr = "Model types lookup table",
                                                desc_1L_chr = "A lookup table of abbreviations to describe the different model types supported by TTU functions"),
                    tibble::tibble(short_name_chr = c("coefs","hetg",
                                                      "dnst","sctr_plt",
                                                      "sim_dnst","sim_sctr",
                                                      "cnstrd_dnst","cnstrd_sctr_plt",
                                                      "cnstrd_sim_dnst","cnstrd_sim_sctr"),
                                   long_name_chr = c("population level effects",
                                                     "group level effects",
                                                     "comparative densities of observed data and predictions using mean model parameter values",
                                                     "comparative scatter plot of observed and predictions using mean model parameter values",
                                                     "comparative densities of observed data and predictions using sampled model parameter values",
                                                     "comparative scatter plot of observed and predictions using sampled model parameter values",
                                                     "comparative densities of observed data and predictions using mean model parameter values and transformation of out of range predictions to upper and lower bounds",
                                                     "comparative scatter plot of observed and predictions using mean model parameter values and transformation of out of range predictions to upper and lower bounds",
                                                     "comparative densities of observed data and predictions using sampled model parameter values and transformation of out of range predictions to upper and lower bounds",
                                                     "comparative scatter plot of observed and predictions using sampled model parameter values and transformation of out of range predictions to upper and lower bounds")) %>%
                      ready4fun::make_pkg_ds_ls(db_1L_chr = "plt_types_lup",
                                                title_1L_chr = "Model plot types lookup table",
                                                desc_1L_chr = "A lookup table of abbreviations to describe the different model plot types supported by TTU functions"),
                    # ADD TO PLOTS TABLE: # AUTOPLT # LNR_CMPRSN  # PRED_DNSTY # PRED_SCTR # SIM_DNSTY # BORUTA_VAR_IMP # RF_VAR_IMP
                    tibble::tibble(rprt_nms_chr = "AAA_TTU_MDL_CTG",
                                   title_chr = "Results supplement: longitudinal transfer to utility models.",
                                   paths_to_rmd_dir_1L_chr = NA_character_,
                                   pkg_dirs_chr = "Markdown",
                                   packages_chr = "TTU",
                                   nms_of_rmd_chr = "Report_TS_Mdls.RMD",
                                   rltv_paths_to_outp_yaml_chr = "_output.yml") %>%
                      tibble::add_case(rprt_nms_chr = "AAA_PMRY_ANLYS_MTH",
                                       title_chr = "Methods supplement: Main analysis algorithm",
                                       paths_to_rmd_dir_1L_chr = NA_character_,
                                       pkg_dirs_chr = "Markdown",
                                       packages_chr = "TTU",
                                       nms_of_rmd_chr = "Analyse.Rmd") %>%
                      tibble::add_case(rprt_nms_chr = "AAA_RPRT_WRTNG_MTH",
                                       title_chr = "Methods supplement: algorithm to auto-generate reports.",
                                       paths_to_rmd_dir_1L_chr = NA_character_,
                                       pkg_dirs_chr = "Markdown",
                                       packages_chr = "TTU",
                                       nms_of_rmd_chr = "Report.Rmd") %>%
                      ready4fun::make_pkg_ds_ls(db_1L_chr = "rprt_lup",
                                                title_1L_chr = "Report types lookup table",
                                                desc_1L_chr = "A lookup table of the different report types supported by TTU functions"))
z <- ready4pack::make_pt_ready4pack_manifest(x,
                                             constructor_r3 = y,
                                             pkg_ds_ls_ls = datasets_ls) %>%
  ready4pack::ready4pack_manifest()
z <- ready4::author(z)
# usethis::use_dev_package("youthvars",
#                          type = "Suggests",#D?
#                          remote = "ready4-dev/youthvars")
# usethis::use_dev_package("scorz",
#                          type = "Depends",
#                          remote = "ready4-dev/scorz")
usethis::use_dev_package("specific",
                         type = "Depends",
                         remote = "ready4-dev/specific")
devtools::build_vignettes()
# usethis::use_package("readr")
# MANUAL DELETION OF TRAILING INCLUDE
# usethis::use_dev_package("ready4",
#                          type = "Depends",
#                          remote = "ready4-dev/ready4")
# usethis::use_dev_package("specific",
#                          type = "Depends",
#                          remote = "ready4-dev/scorz")
# usethis::use_package("rgl")
# piggyback::pb_new_release("ready4-dev/TTU",
#                           tag = paste0("v",desc::desc_get_version()),
#                           body = "Version implemented following significant redevelopment of package dependencies.",
#                           prerelease = F)
# devtools::build_vignettes()
