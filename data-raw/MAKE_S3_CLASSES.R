library(magrittr)
s3_classes_to_make_tb <- ready4class::ready4_constructor_tbl() %>%
  dplyr::bind_rows(tibble::tribble(
    ~ make_s3_lgl, ~ name_stub_chr, ~ pt_ls, ~ pt_chkr_pfx_ls, ~ pt_ns_ls, ~ vals_ls, ~ allowed_vals_ls, ~ min_max_vals_ls, ~ start_end_vals_ls, ~ class_desc_chr, ~ parent_class_chr, ~ slots_ls, ~ meaningful_nms_ls, ~inc_clss_ls,
    TRUE, "phq9", list("numeric"), list("is."),list("base"),NULL, NULL,
    list(c(0,27)),
    NULL, "Readyforwhatsnext S3 class for PHQ-9", NA_character_, NULL, NULL, NULL))

data("prototype_lup",package = "ready4class")
prototype_lup <- tibble::add_case(prototype_lup, type_chr = "integer", val_chr = "NA_integer_", pt_ns_chr = "base",fn_to_call_chr="",default_val_chr = "NA_integer_", old_class_lgl = F)
data("abbreviations_lup")
classes_to_make_tb <- s3_classes_to_make_tb
name_pfx_1L_chr <- "firstbounce_"
pkg_dss_tb <- ready4fun::write_abbr_lup(short_name_chr = paste0(name_pfx_1L_chr,classes_to_make_tb$name_stub_chr),
                                        long_name_chr = c(classes_to_make_tb$class_desc_chr),
                                        custom_plural_ls = NULL,
                                        no_plural_chr = NA_character_,
                                        url_1L_chr = "https://readyforwhatsnext.github.io/readyforwhatsnext/",
                                        seed_lup = abbreviations_lup)
data("abbreviations_lup")

classes_to_make_tb %>%
  ready4class::write_classes_and_make_lup(dev_pkg_ns_1L_chr = ready4fun::get_dev_pkg_nm(),
                                          name_pfx_1L_chr = name_pfx_1L_chr,
                                          output_dir_1L_chr = "R",
                                          file_exists_cdn_1L_chr = "overwrite",
                                          abbreviations_lup = abbreviations_lup,
                                          init_class_pt_lup = prototype_lup)
