classes_to_make_tb <- ready4class::ready4_constructor_tbl() %>%
  dplyr::bind_rows(tibble::tribble(
    ~ make_s3_lgl, ~ name_stub_chr, ~ pt_ls, ~ pt_chkr_pfx_ls, ~ pt_ns_ls, ~ vals_ls, ~ allowed_vals_ls, ~ min_max_vals_ls, ~ start_end_vals_ls, ~ class_desc_chr, ~ parent_class_chr, ~ slots_ls, ~ meaningful_nms_ls, ~ inc_clss_ls, ~ asserts_ls,
    TRUE, "predictors_lup", list("tibble"), list("is_"),list("tibble"),list(short_name_chr = "character(0)",
                                                                            long_name_chr = "character(0)",
                                                                            min_val_dbl = "numeric(0)",
                                                                            max_val_dbl = "numeric(0)",
                                                                            class_chr = "character(0)",
                                                                            increment_dbl = "numeric(0)",
                                                                            class_fn_chr = "character(0)",
                                                                            mdl_scaling_dbl = "numeric(0)",
                                                                            covariate_lgl = "logical(0)"), NULL,NULL, NULL, "TTU S3 class for candidate predictors lookup table", NA_character_, NULL, NULL, NULL, NULL)
  )
name_pfx_1L_chr <- "TTU_"
