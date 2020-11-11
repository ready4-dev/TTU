devtools::load_all()
path_to_rprt_dir_1L_chr = normalizePath("../../../../../Data/Project/Utility_Models/Reports")
write_rndrd_rprt(path_to_RMD_1L_chr = 'vignettes/test.Rmd', #'vignettes/Child_RMDs/_TEST_EXPANSION.Rmd',# 'vignettes/Child_RMDs/_AQoL_6D_Replication_Mixed.Rmd',#
                 path_to_rprt_dir_1L_chr = path_to_rprt_dir_1L_chr,
                 file_nm_1L_chr = "SI",#"main_models"
                 params_ls = list(output_type_1L_chr = "HTML",
                                  mdl_smry_dir_1L_chr = path_to_rprt_dir_1L_chr))





