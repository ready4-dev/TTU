## -----------------------------------------------------------------------------
library(ready4)
library(ready4use)
library(youthvars)
library(scorz)
library(TTU)

## -----------------------------------------------------------------------------
A <- Ready4useDyad(ds_tb = Ready4useRepos(dv_nm_1L_chr = "fakes",
                                          dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/HJXYKQ",
                                          dv_server_1L_chr = "dataverse.harvard.edu") %>%
                     ingest(fls_to_ingest_chr = c("ymh_clinical_tb"),
                            metadata_1L_lgl = F) %>%
                     youthvars::transform_raw_ds_for_analysis(),
                   dictionary_r3 = Ready4useRepos(dv_nm_1L_chr = "TTU", 
                                                  dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0", 
                                                  dv_server_1L_chr = "dataverse.harvard.edu") %>%
                     ingest(fls_to_ingest_chr = c("dictionary_r3"),
                            metadata_1L_lgl = F)) %>%
  renew(type_1L_chr = "label")

## -----------------------------------------------------------------------------
A <- YouthvarsSeries(a_Ready4useDyad = A,
                     id_var_nm_1L_chr = "fkClientID",
                     timepoint_var_nm_1L_chr = "round",
                     timepoint_vals_chr = levels(procureSlot(A,
                                                             "ds_tb")$round))

## -----------------------------------------------------------------------------
A <- TTUProject(a_ScorzProfile = ScorzAqol6Adol(a_YouthvarsProfile = A))

## -----------------------------------------------------------------------------
# Make into renew method
A <- renewSlot(A, "a_ScorzProfile")
A <- renewSlot(A, "b_SpecificParameters", SpecificConverter(a_ScorzProfile = A@a_ScorzProfile) %>%
                 metamorphose() %>%
                 procureSlot("b_SpecificParameters"))

## -----------------------------------------------------------------------------
A <- renewSlot(A, "b_SpecificParameters@predictors_lup", Ready4useRepos(dv_nm_1L_chr = "TTU", 
                                                                        dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0", 
                                                                        dv_server_1L_chr = "dataverse.harvard.edu") %>%
                 ingest(fls_to_ingest_chr = c("predictors_r3"),
                        metadata_1L_lgl = F)) 

## -----------------------------------------------------------------------------
exhibitSlot(A, "b_SpecificParameters@predictors_lup")

## -----------------------------------------------------------------------------
A <- renewSlot(A, "b_SpecificParameters@depnt_var_min_max_dbl", c(0.03,1)) %>% # Inherit From TTUAqolAdol
  renewSlot("b_SpecificParameters@candidate_predrs_chr", c("BADS","GAD7", "K6", "OASIS", "PHQ9", "SCARED")) %>%
  renewSlot("b_SpecificParameters@candidate_covars_chr", c("d_sex_birth_s", "d_age",  "d_sexual_ori_s", 
                                                           "d_studying_working", "c_p_diag_s", "c_clinical_staging_s",
                                                           "SOFAS")) %>%
  renewSlot("b_SpecificParameters@descv_var_nms_chr", c("d_age","Gender","d_relation_s", "d_sexual_ori_s", 
                                                        "Region", "d_studying_working", "c_p_diag_s", 
                                                        "c_clinical_staging_s","SOFAS")) %>%
  renewSlot("b_SpecificParameters@msrmnt_date_var_nm_1L_chr", "d_interview_date") 

## -----------------------------------------------------------------------------
A <-  renewSlot(A, "b_SpecificParameters@fake_1L_lgl", T) 

## -----------------------------------------------------------------------------
A <- renewSlot(A, "c_SpecificProject", SpecificModels(a_YouthvarsProfile = A@a_ScorzProfile@a_YouthvarsProfile,
                                                      b_SpecificParameters = A@b_SpecificParameters,
                                                      paths_chr = tempdir())) 

## -----------------------------------------------------------------------------
A <- ratifySlot(A, "c_SpecificProject")

## -----------------------------------------------------------------------------
A <- renewSlot(A, "c_SpecificProject", 
               authorSlot(A, "c_SpecificProject", what_1L_chr = "workspace"))

## -----------------------------------------------------------------------------
A <- renewSlot(A, "c_SpecificProject",
               authorSlot(A, "c_SpecificProject", what_1L_chr = "descriptives",
                          digits_1L_int = 3L))

## -----------------------------------------------------------------------------
A <- renewSlot(A, "c_SpecificProject",
               investigateSlot(A, "c_SpecificProject",
                               depnt_var_max_val_1L_dbl = 0.99,# Parent method argument
                               session_ls = sessionInfo()))

## -----------------------------------------------------------------------------
A <- renewSlot(A, "c_SpecificProject",
               renew(procureSlot(A, "c_SpecificProject"),
                     new_val_xx = c("GLM_GSN_LOG", "OLS_CLL"),
                     type_1L_chr = "results",
                     what_1L_chr = "prefd_mdls"))

## -----------------------------------------------------------------------------
A <- renewSlot(A, "c_SpecificProject",
               investigateSlot(A,"c_SpecificProject"))

## -----------------------------------------------------------------------------
A <- renewSlot(A, "c_SpecificProject",
               renew(procureSlot(A, "c_SpecificProject"),
                     new_val_xx = "SOFAS",
                     type_1L_chr = "results",
                     what_1L_chr = "prefd_covars"))

## -----------------------------------------------------------------------------
A <- renewSlot(A, "c_SpecificProject",
               investigateSlot(A, "c_SpecificProject"))

