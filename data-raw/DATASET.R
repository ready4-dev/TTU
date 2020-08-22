##
## This script creates the data files embedded with this package.
# 1. Load magrittr package to that the pipe operator ("%>%") can be used in this script.
library(magrittr)
#
# 2. Specify package name
pkg_nm_chr <- "FBaqol"
#
# 3. Create "fns", "gnrcs" and "mthds" sub-directories.
undocumented_fns_dir_chr <- paste0("data-raw/",c("fns","gnrcs","mthds"))
paths_ls <- undocumented_fns_dir_chr %>% purrr::map(~{
  if(!dir.exists(.x))
    dir.create(.x)
  ready4fun::read_fns(.x)
}) %>% stats::setNames(c("fns","gnrcs","mthds"))
# 4. Create a lookup table of abbreviations used in this package and save it as a package dataset (data gets saved in the data directory, documentation script is created in R directory).
data("abbreviations_lup",package = "ready4class")
# ready4fun::make_abbr_lup_tb(short_name_chr_vec = NA_character_,
#                             long_name_chr_vec = NA_character_,
#                             url_chr = NA_character_,
#                             pkg_nm_chr = pkg_nm_chr,
#                             seed_lup = abbreviations_lup)
# data("abbreviations_lup")
# 5. Create function types and generics look-up tables
# 5.1 Create a lookup table of function types used in this package and save it as a package dataset (data gets saved in the data directory, documentation script is created in R directory).
data("fn_type_lup_tb",package = "ready4class")
# fn_type_lup_tb %>%
#   dplyr::bind_rows(tibble::tibble(fn_type_nm_chr = NA_character_,
#                                   fn_type_desc_chr = NA_character_,
#                                   first_arg_desc_chr = NA_character_,
#                                   second_arg_desc_chr = NA_character_,
#                                   is_generic_lgl = F)) %>% # Add to ready4fun template.
#   dplyr::arrange(fn_type_nm_chr) %>%
#   ready4fun::make_and_doc_fn_type_R(pkg_nm_chr = pkg_nm_chr,
#                                     url_chr = NA_character_,
#                                     abbreviations_lup = abbreviations_lup)
data("fn_type_lup_tb")
# 5.2 Create a look-up table of the generics used in this package.
generics_lup_tb <- NULL
# ready4fun::make_and_doc_generics_tb_R(generic_nm_chr = NA_character_,
#                                       description_chr = NA_character_,
#                                       pkg_nm_chr = pkg_nm_chr,
#                                       url_chr = NA_character_,
#                                       abbreviations_lup = abbreviations_lup)
# data("generics_lup_tb")
#
# 6. Create a table of all functions to document
all_fns_dmt_tb <- make_all_fns_dmt_tb(paths_ls = paths_ls,
                                      undocumented_fns_dir_chr = undocumented_fns_dir_chr,
                                      custom_dmt_ls = list(details_ls = NULL,
                                                           export_ls = list(force_true_chr_vec = NA_character_,
                                                                            force_false_chr_vec = NA_character_),
                                                           args_ls_ls = NULL),
                                      fn_type_lup_tb = fn_type_lup_tb,
                                      generics_lup_tb = generics_lup_tb,
                                      abbreviations_lup = abbreviations_lup)
##
aqol6d_from_8d_coeffs_lup_tb <- tibble::tribble(
  ~var_name_chr, ~coeff_dbl,
  "vD1", 0.0719264,
  "vD2", 0.1027818,
  "vD3", 0.2519563,
  "vD4", 0.3201172,
  "vD5", 0.1288289,
  "vD6", 0.2052164,
  "Constant", - 0.0444493
)
disutilities_lup_tb <- tibble::tribble(
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
)
dim_sclg_constant_lup_tb <- tibble::tribble(
  ~Dimension_chr, ~Constant_dbl,
  "IL",-0.978,
  "RL", -0.923,
  "MH", -0.983,
  "COP", -0.930,
  "P", -0.96,
  "SEN", -0.851
)
itm_wrst_wghts_lup_tb <- tibble::tribble(
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
)
## 7 Document.
# Write documented methods to R directory.
## Note files to be rewritten cannot be open in RStudio.
ready4fun::write_and_doc_fn_fls_R(all_fns_dmt_tb,
                                  r_dir_chr = "R")
#
# Update Description file with imported packages.
ready4fun::write_ns_imps_to_desc()
##
## 8. Run script to make package classes.
# source("data-raw/MAKE_CLASSES.R")
# prototype_lup <- prototype_lup %>%
#   ready4_class_pt_lup()
# usethis::use_data(prototype_lup,overwrite = T, internal = T)
# ## 9. Remake the classes we previously created, this time using the new, preferred make_and_update method, which appends the metadata on the new classes to our instance of the ready4_class_pt_lup class.
# prototype_lup <- make_and_update(classes_to_make_tb,
#                                  dev_pckg_namespace = dev_pckg_namespace,
#                                  name_prefix = name_prefix,
#                                  output_dir = "R",
#                                  file_exists_logic = "overwrite")
# ## 10. Update the internal system data.
# usethis::use_data(prototype_lup,overwrite = T, internal = T)
##
# 11. Create vignettes
usethis::use_vignette("ready4class")
devtools::document()
ready4fun::write_ns_imps_to_desc()
