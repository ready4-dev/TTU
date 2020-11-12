devtools::load_all()
write_rndrd_rprt("vignettes/Child_RMDs",
                 nm_of_RMD_1L_chr = "_Mdls_Report.RMD",
                 params_ls = list(mdl_smry_dir_1L_chr = "../Output",
                                  output_type_1L_chr = "Word",
                                  section_type_1L_chr = "#"),
                 rltv_path_to_outpt_yaml_1L_chr = "output_yml",
                 paths_to_fls_to_copy_chr = list.files(path_to_RMD_dir_1L_chr, full.names = T)[c(15:17,22:23)],
                 path_to_write_fls_to_1L_chr = normalizePath("../../../../../Data/Project/Utility_Models"),
                 nm_of_rprt_dir_1L_chr = "Markdown",
                 path_to_outpt_rtrp_1L_chr = normalizePath("../../../../../Data/Project/Utility_Models/Reports"),
                 file_nm_1L_chr = "SI")


