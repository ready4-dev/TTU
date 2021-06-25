# get model ingredients
ingredients_ls <- youthu::get_mdl_from_dv("mdl_ingredients",
                                          dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/612HDC")
# Inspect which predictors are used in dataset models
ingredients_ls$predictors_tb
# Predict
data_tb <- make_fake_eq5d_ds(prop_with_fup_data_1L_dbl = 0.6, seed_1L_int = 1, force_attach_1L_lgl = F)
data_tb <- data_tb%>%
  dplyr::filter(uid %in% sample(data_tb$uid, 500))
data_tb %>% add_utl_predn_to_new_ds(ingredients_ls = ingredients_ls,
                        mdl_nm_1L_chr = mdl_nms_chr[3],
                        new_data_is_1L_chr = "Simulated",
                        predr_vars_nms_chr = c(k10 = "k10_int"),
                        round_var_nm_1L_chr = "Timepoint",
                        round_bl_val_1L_chr = "BL")
# FOR YOUTHU
# Get names of models using predictors present in new data
# Review summary information on each model
## Retrieve and review model catalogue
# outp_smry_ls <- readRDS(normalizePath("Fake/Output/I_ALL_OUTPUT_.RDS"))
# rmarkdown::render(normalizePath("Fake/Markdown/Report_TS_Mdls.Rmd"),
#                   params = list(outp_smry_ls = outp_smry_ls,
#                                 output_type_1L_chr = "PDF"),
#                   output_format = "bookdown::pdf_book",
#                   output_file = "test.pdf")

