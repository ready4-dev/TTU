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
