## ----eval=F-------------------------------------------------------------------
#  @article {Hamilton2021.07.07.21260129,
#  	author = {Hamilton, Matthew P and Gao, Caroline X and Filia, Kate M and Menssink, Jana M and Sharmin, Sonia and Telford, Nic and Herrman, Helen and Hickie, Ian B and Mihalopoulos, Cathrine and Rickwood, Debra J and McGorry, Patrick D and Cotton, Sue M},
#  	title = {Predicting Quality Adjusted Life Years in young people attending primary mental health services},
#  	elocation-id = {2021.07.07.21260129},
#  	year = {2021},
#  	doi = {10.1101/2021.07.07.21260129},
#  	publisher = {Cold Spring Harbor Laboratory Press},
#  	URL = {https://www.medrxiv.org/content/early/2021/07/12/2021.07.07.21260129},
#  	eprint = {https://www.medrxiv.org/content/early/2021/07/12/2021.07.07.21260129.full.pdf},
#  	journal = {medRxiv}
#  }

## ----eval=FALSE---------------------------------------------------------------
#  @software{hamilton_matthew_2022_6212704,
#    author       = {Hamilton, Matthew and
#                    Gao, Caroline},
#    title        = {{Complete study program to reproduce all steps from
#                     data ingest through to results dissemination for a
#                     study to map mental health measures to AQoL-6D
#                     health utility}},
#    month        = feb,
#    year         = 2022,
#    note         = {{Matthew Hamilton and Caroline Gao  (2022).
#                     Complete study program to reproduce all steps from
#                     data ingest through to results dissemination for a
#                     study to map mental health measures to AQoL-6D
#                     health utility. Zenodo.
#                     https://doi.org/10.5281/zenodo.6116077. Version
#                     0.0.9.3}},
#    publisher    = {Zenodo},
#    version      = {0.0.9.3},
#    doi          = {10.5281/zenodo.6212704},
#    url          = {https://doi.org/10.5281/zenodo.6212704}
#  }

## ----message=FALSE, warning=FALSE---------------------------------------------
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
A <- renewSlot(A, "a_ScorzProfile")

## ----echo=FALSE---------------------------------------------------------------
# Make next bit into renew method

## -----------------------------------------------------------------------------
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
exhibitSlot(A, "b_SpecificParameters@predictors_lup",
         scroll_box_args_ls = list(width = "100%"))

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

