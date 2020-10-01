## This script creates the data files embedded with this package.
## As it involves interaction, it must be run  in four separate parts.
##
## PART ONE
##
# 1. Load magrittr package to that the pipe operator ("%>%") can be used in this script.
library(magrittr)
# 2. Create "fns", "gnrcs" and "mthds" sub-directories.
ready4fun::write_fn_type_dirs()
# 3. MANUAL STEP. Write all your functions to R files in the new "fns" directory.
##
## PART TWO
##
# 4. Set-up package structure
options(usethis.description = list(
  Package = ready4fun::get_dev_pkg_nm(),
  Title =  "Map measures routinely collected in youth mental health services to AQOL 6D Health Utility.",
  Description = "Tools for mapping measures routinely collected in youth mental health services to AQOL 6D Health Utility. Part of the First Bounce model of primary youth mental health services.",
  `Authors@R` = c(utils::person(given = "Caroline",family = "Gao",email = "caroline.gao@orygen.org.au", role = c("aut"),comment = c(ORCID = "0000-0002-0987-2759")),
                  utils::person(given = "Matthew",family = "Hamilton",email = "matthew.hamilton@orygen.org.au", role = c("aut", "cre"),comment = c(ORCID = "0000-0001-7407-9194")),
                  utils::person("Orygen", role = c("cph", "fnd")),
                  utils::person("Headspace", role = c( "fnd")),
                  utils::person("National Health and Medical Research Council", role = c( "fnd"))
  ),
  License = usethis::use_gpl3_license("Orygen"),
  URL = c("https://github.com/orygen/FBaqol")
))
# Deletes contents of R directory and resets DESCRIPTION and NAMESPACE files.
ready4fun::write_pkg_setup_fls(incr_ver_1L_lgl = F,
                               delete_contents_of_R_dir = T,
                               copyright_holders_chr = "Orygen",
                               use_travis_1L_lgl = F,
                               path_to_pkg_logo_1L_chr = "man/figures/default.png",
                               github_repo_1L_chr = "orygen/FBaqol",
                               lifecycle_stage_1L_chr = "experimental")
# PAUSE FOR INTERACTIVE
##
## PART THREE
##
# 5. Create a lookup table of abbreviations used in this package and save it as a package dataset (data gets saved in the data directory, documentation script is created in R directory).
data("abbreviations_lup",package = "ready4use")
pkg_dss_tb <- ready4fun::write_abbr_lup(short_name_chr = c("adol","aqol","aqol6d","dim","disv","eq","scrg","unscrd"),
                            long_name_chr = c("adolescent",
                                                  "Assessment of Quality of Life health utility",
                                                  "Assessment of Quality of Life Six Dimension health utility",
                                                  "dimension",
                                                  "disvalue",
                                                  "equation",
                                                  "scoring",
                                                  "unscored"),
                            url_1L_chr = NA_character_,
                            seed_lup = abbreviations_lup)
data("abbreviations_lup")
# 5. Create function types and generics look-up tables
# 5.1 Create a lookup table of function types used in this package and save it as a package dataset (data gets saved in the data directory, documentation script is created in R directory).
data("fn_type_lup_tb",package = "ready4use")
# ready4fun::get_new_fn_types(abbreviations_lup = abbreviations_lup,
#                             fn_type_lup_tb = fn_type_lup_tb)
pkg_dss_tb <- fn_type_lup_tb %>%
  ready4fun::add_rows_to_fn_type_lup(fn_type_nm_chr = c("Calculate","Extract","Impute","Randomise","Reorder","Reset","Scramble"),
                                  fn_type_desc_chr = c("Calculates a numeric value.",
                                                       "Extracts data from an object.",
                                                       "Imputes data.",
                                                       "Randomly samples from data.",
                                                       "Reorders an object to conform to a pre-specified schema.",
                                                       "Resets the value of an object to a default.",
                                                       "Randomly reorders an object."),
                                  is_generic_lgl = F,
                                  is_method_lgl = F) %>% # Add to ready4fun template.
  dplyr::arrange(fn_type_nm_chr) %>%
  ready4fun::write_dmtd_fn_type_lup(url_1L_chr = NA_character_,
                                    abbreviations_lup = abbreviations_lup,
                                    pkg_dss_tb = pkg_dss_tb)
data("fn_type_lup_tb")

#
# 6. Create a table of all functions to document
fns_dmt_tb <- ready4fun::make_dmt_for_all_fns(paths_ls = ready4fun::make_fn_nms()[1],
                                              undocumented_fns_dir_chr = ready4fun::make_undmtd_fns_dir_chr()[1],
                                              custom_dmt_ls = list(details_ls = NULL,
                                                                   inc_for_main_user_lgl_ls = list(force_true_chr = c("add_aqol_items_tbs_ls","add_aqol_scores_tbs_ls",
                                                                                                                      "add_aqol6dU_to_aqol6d_items_tb_tb","add_aqol6dU_to_tbs_ls",
                                                                                                                      "add_corrs_and_uts_to_tbs_ls_ls","add_dmn_disu_to_aqol6d_items_tb_tb",
                                                                                                                      "add_dmn_scores_to_aqol6d_items_tb_tb", "add_domain_unwtd_tots_tb",
                                                                                                                      "add_itm_disu_to_aqol6d_itms_tb_tb","add_labels_to_aqol6d_tb",
                                                                                                                      "add_uids_to_tbs_ls","calculate_aqol6d_d1_disu_dbl",
                                                                                                                      "calculate_aqol6d_d2_disu_dbl","calculate_aqol6d_d3_disu_dbl",
                                                                                                                      "calculate_aqol6d_d4_disu_dbl","calculate_aqol6d_d5_disu_dbl",
                                                                                                                      "calculate_aqol6d_d6_disu_dbl","calculate_aqol6dU_dbl",
                                                                                                                      "extract_g_legend_1L_chr","force_min_max_and_int_cnstrs_tb",
                                                                                                                      "force_vec_to_sum_to_int","impute_miss_itms_in_aqol6d_items_tb_tb",
                                                                                                                      #"make_aqol_items_props_tbs_ls",
                                                                                                                      "make_aqol6d_fns_ls",
                                                                                                                      "make_aqol6d_items_tb", "make_correlated_data_tb",
                                                                                                                      "make_corstars_tbl_xx","make_dim_sclg_cons_dbl",
                                                                                                                      "make_domain_items_ls","make_item_wrst_wghts_ls_ls",
                                                                                                                      "make_pdef_corr_mat_mat", "make_synth_series_tbs_ls",
                                                                                                                      "make_vec_with_sum_of_int", "randomise_ptl_fup_fct",
                                                                                                                      "reorder_tbs_for_target_cors","replace_var_vals_with_missing_tbl",
                                                                                                                      "scramble_xx","transform_raw_aqol_tb_to_aqol6d_tb",
                                                                                                                      "write_results_csv"),
                                                                                                   force_false_chr = NA_character_),
                                                                   args_ls_ls = NULL),
                                              fn_type_lup_tb = fn_type_lup_tb,
                                              abbreviations_lup = abbreviations_lup)
##
pkg_dss_tb <- tibble::tribble(
  ~var_name_chr, ~coeff_dbl,
  "vD1", 0.0719264,
  "vD2", 0.1027818,
  "vD3", 0.2519563,
  "vD4", 0.3201172,
  "vD5", 0.1288289,
  "vD6", 0.2052164,
  "Constant", - 0.0444493
) %>% ready4fun::write_and_doc_ds(db_1L_chr = "aqol6d_from_8d_coeffs_lup_tb",
                                    title_1L_chr = "Model 2A Coefficients To Weight AQoL6D",
                                    desc_1L_chr = "Coefficients for model to predict AQoL-6D utility score from AQoL-8D. The optimal model is Model 2A (see Richardson et al (2011, 18-19)*/",
                                    url_1L_chr = "https://www.aqol.com.au/index.php/scoring-algorithms",
                              abbreviations_lup = abbreviations_lup,
                              pkg_dss_tb = pkg_dss_tb)
pkg_dss_tb <- tibble::tribble(
  ~Question_chr, ~Answer_1_dbl, ~Answer_2_dbl, ~Answer_3_dbl, ~Answer_4_dbl, ~Answer_5_dbl, ~Answer_6_dbl,
  "Q1", 0, 0.073, 0.435, 0.820, 1, NA_real_,
  "Q2", 0, 0.033, 0.240, 0.471, 0.840,1,
  "Q3", 0, 0.041, 0.251, 0.570, 0.830, 1,
  "Q4", 0, 0.040, 0.297, 0.797, 1, NA_real_,
  "Q5", 0, 0.074, 0.461, 0.841, 1, NA_real_,
  "Q6", 0, 0.193, 0.759, 1, NA_real_,NA_real_,
  "Q7", 0, 0.197, 0.648, 1, NA_real_, NA_real_,
  "Q8", 0, 0.133, 0.392, 0.838, 1, NA_real_,
  "Q9", 0, 0.142, 0.392, 0.824, 1, NA_real_,
  "Q10", 0, 0.097, 0.330, 0.784, 1, NA_real_,
  "Q11", 0, 0.064, 0.368, 0.837, 1, NA_real_,
  "Q12", 0, 0.056, 0.338, 0.722, 1, NA_real_,
  "Q13", 0, 0.055, 0.382, 0.774, 1, NA_real_,
  "Q14", 0, 0.057, 0.423, 0.826, 1, NA_real_,
  "Q15", 0, 0.133, 0.642, 1, NA_real_,NA_real_,
  "Q16", 0, 0.200, 0.758, 1, NA_real_, NA_real_,
  "Q17", 0, 0.072, 0.338, 0.752, 1, NA_real_,
  "Q18", 0, 0.033, 0.223, 0.621, 0.843, 1,
  "Q19", 0, 0.024, 0.205, 0.586, 0.826, 1,
  "Q20", 0, 0.187, 0.695, 1, NA_real_,NA_real_
) %>% ready4fun::write_and_doc_ds(db_1L_chr = "disvalues_lup_tb",
                                    title_1L_chr = "AQoL6D item disvalues lookup table",
                                    desc_1L_chr = "Disutility weights for individual AQoL6D items.",
                                    url_1L_chr = "https://www.aqol.com.au/index.php/scoring-algorithms",
                                    abbreviations_lup = abbreviations_lup,
                                  pkg_dss_tb = pkg_dss_tb)
pkg_dss_tb <- tibble::tribble(
  ~Dimension_chr, ~Constant_dbl,
  "IL",-0.978,
  "RL", -0.923,
  "MH", -0.983,
  "COP", -0.930,
  "P", -0.96,
  "SEN", -0.851
) %>% ready4fun::write_and_doc_ds(db_1L_chr = "dim_sclg_constant_lup_tb",
                                    title_1L_chr = "AQoL6D dimension scaling constants lookup table",
                                    desc_1L_chr = "Scaling constants for each dimension of AQoL6D.",
                                    url_1L_chr = "https://www.aqol.com.au/index.php/scoring-algorithms",
                                    abbreviations_lup = abbreviations_lup,
                                  pkg_dss_tb = pkg_dss_tb)
pkg_dss_tb <- tibble::tribble(
  ~Question_chr, ~Worst_Weight_dbl,
  "Q1", 0.385412,
  "Q2", 0.593819,
  "Q3", 0.630323,
  "Q4", 0.794888,
  "Q5", 0.64303,
  "Q6", 0.697742,
  "Q7", 0.508658,
  "Q8", 0.640377,
  "Q9", 0.588422,
  "Q10", 0.648748,
  "Q11", 0.71122,
  "Q12", 0.415694,
  "Q13", 0.636994,
  "Q14", 0.773296,
  "Q15", 0.631833,
  "Q16", 0.767573,
  "Q17", 0.652241,
  "Q18", 0.580696,
  "Q19", 0.463022,
  "Q20", 0.604613
) %>% ready4fun::write_and_doc_ds(db_1L_chr = "itm_wrst_wghts_lup_tb",
                                    title_1L_chr = "AQoL6D item worst weightings lookup table",
                                    desc_1L_chr = "Worst weightings for individual items in AQoL6D.",
                                    url_1L_chr = "https://www.aqol.com.au/index.php/scoring-algorithms",
                                    abbreviations_lup = abbreviations_lup,
                                  pkg_dss_tb = pkg_dss_tb)
pkg_dss_tb <- read.csv("data-raw/AQoL_6D_Dim_Scaling.csv", stringsAsFactors = F, fileEncoding="UTF-8-BOM") %>%
  ready4fun::write_and_doc_ds(db_1L_chr = "adol_dim_scalg_eqs_lup",
                              title_1L_chr = "AQoL6D item worst weightings lookup table",
                              desc_1L_chr = "Dimension scaling equations for adolescent version of AQoL6D scoring algorithm.",
                              url_1L_chr = "https://www.aqol.com.au/index.php/scoring-algorithms",
                              abbreviations_lup = abbreviations_lup,
                              pkg_dss_tb = pkg_dss_tb)
pkg_dss_tb <- read.csv("vignettes/Data/aqol_valid_stata.csv") %>%
  ready4fun::write_and_doc_ds(db_1L_chr = "syn_pop_with_STATA_adults_scoring_tb",
                              title_1L_chr = "STATA comparison validation synthetic population",
                              desc_1L_chr = "Synthetic population following application of STATA adult scoring algorithm.",
                              url_1L_chr = "https://www.aqol.com.au/index.php/scoring-algorithms",
                              abbreviations_lup = abbreviations_lup,
                              pkg_dss_tb = pkg_dss_tb)
## 7 Document.
# Write documented methods to R directory.
## Note files to be rewritten cannot be open in RStudio.
ready4fun::write_and_doc_fn_fls(fns_dmt_tb,
                                r_dir_1L_chr = "R",
                                dev_pkgs_chr = c("ready4fun","ready4class","ready4use"),
                                update_pkgdown_1L_lgl = T)
##
## PART FOUR
##
# Update Description file with imported packages.
pkgdown::build_site()
##
## Add, Commit and Push
