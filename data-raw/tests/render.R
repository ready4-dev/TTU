devtools::load_all()
# write_rndrd_rprt("vignettes/Child_RMDs",
#                  nm_of_RMD_1L_chr = "_AQoL_6D_Replication_Mixed.RMD",
#                  params_ls = list(mdl_smry_dir_1L_chr = "../Output",
#                                   output_type_1L_chr = "Word",
#                                   path_to_data_1L_chr = normalizePath("vignettes/Data/fake_pop_tb.rds")),
#                  rltv_path_to_outp_yaml_1L_chr = "_output.yml",
#                  paths_to_fls_to_copy_chr = list.files("vignettes/Child_RMDs", full.names = T)[c(1,5,9,10)],
#                  path_to_write_fls_to_1L_chr = normalizePath("../../../../../Data/Project/Utility_Models"),
#                  nm_of_rprt_dir_1L_chr = "Markdown",
#                  path_to_outp_rtrp_1L_chr = normalizePath("../../../../../Data/Project/Utility_Models/Reports"),
#                  file_nm_1L_chr = "Mixed_Model_Code")
ready4show::write_rndrd_rprt("vignettes/Child_RMDs",
                 nm_of_RMD_1L_chr = "_Mdls_Report.RMD",
                 params_ls = list(mdl_smry_dir_1L_chr = "../Output",
                                  output_type_1L_chr = "PDF",
                                  section_type_1L_chr = "#"),
                 rltv_path_to_outp_yaml_1L_chr = "_output.yml",
                 paths_to_fls_to_copy_chr = list.files("vignettes/Child_RMDs", full.names = T)[c(15:17,22:23)],
                 path_to_write_fls_to_1L_chr = normalizePath("../../../../../Data/Project/Utility_Models"),
                 nm_of_rprt_dir_1L_chr = "Markdown",
                 path_to_outp_rtrp_1L_chr = normalizePath("../../../../../Data/Project/Utility_Models/Reports"),
                 file_nm_1L_chr = "SI_Mdls")


