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
library(ready4show)
library(ready4use)
library(youthvars)
library(scorz)
library(TTU)

## ----eval = FALSE-------------------------------------------------------------
#  consent_1L_chr <- "" # Default value - asks for consent prior to writing each file.

## ----echo = FALSE-------------------------------------------------------------
consent_1L_chr <- "Y" # Gives consent to write files without additional requests.

## -----------------------------------------------------------------------------
A <- Ready4useDyad(ds_tb = Ready4useRepos(dv_nm_1L_chr = "fakes", dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/HJXYKQ", dv_server_1L_chr = "dataverse.harvard.edu") %>%
                     ingest(fls_to_ingest_chr = c("ymh_clinical_tb"), metadata_1L_lgl = F) %>% youthvars::transform_raw_ds_for_analysis(),
                   dictionary_r3 = Ready4useRepos(dv_nm_1L_chr = "TTU", dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0", dv_server_1L_chr = "dataverse.harvard.edu") %>%
                     ingest(fls_to_ingest_chr = c("dictionary_r3"), metadata_1L_lgl = F)) %>%
  renew(type_1L_chr = "label")

## -----------------------------------------------------------------------------
A <- YouthvarsSeries(a_Ready4useDyad = A, id_var_nm_1L_chr = "fkClientID", timepoint_var_nm_1L_chr = "round",
                     timepoint_vals_chr = levels(procureSlot(A, "ds_tb")$round))

## -----------------------------------------------------------------------------
A <- TTUProject(a_ScorzProfile = ScorzAqol6Adol(a_YouthvarsProfile = A))
A <- renew(A, what_1L_chr = "utility") 

## -----------------------------------------------------------------------------
A <- renew(A, what_1L_chr = "parameters")

## -----------------------------------------------------------------------------
A <- renew(A, "use_renew_mthd", fl_nm_1L_chr = "predictors_r3", type_1L_chr = "predictors_lup", 
           y_Ready4useRepos = Ready4useRepos(dv_nm_1L_chr = "TTU", dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0", 
                                             dv_server_1L_chr = "dataverse.harvard.edu"),
           what_1L_chr = "parameters")

## -----------------------------------------------------------------------------
exhibit(A, scroll_box_args_ls = list(width = "100%"))

## -----------------------------------------------------------------------------
A <- renew(A, c(0.03,1), type_1L_chr = "range", what_1L_chr = "parameters") %>%
  renew(c("BADS","GAD7", "K6", "OASIS", "PHQ9", "SCARED"),
        type_1L_chr = "predictors_vars", what_1L_chr = "parameters") %>%
  renew(c("d_sex_birth_s", "d_age",  "d_sexual_ori_s", "d_studying_working", "c_p_diag_s", "c_clinical_staging_s", "SOFAS"),     
        type_1L_chr = "covariates", what_1L_chr = "parameters") %>%
  renew(c("d_age","Gender","d_relation_s", "d_sexual_ori_s" ,"Region", "d_studying_working", "c_p_diag_s", "c_clinical_staging_s","SOFAS"), 
        type_1L_chr = "descriptives", what_1L_chr = "parameters") %>%
  renew("d_interview_date", type_1L_chr = "temporal", what_1L_chr = "parameters")

## -----------------------------------------------------------------------------
A <- renew(A, T, type_1L_chr = "is_fake", what_1L_chr = "parameters")

## -----------------------------------------------------------------------------
A <- renew(A, consent_1L_chr = consent_1L_chr, paths_chr = tempdir(), what_1L_chr = "project")

## ----message=FALSE, results='hide', warning=FALSE-----------------------------
A <- author(A, consent_1L_chr = consent_1L_chr, digits_1L_int = 3L, what_1L_chr = "descriptives")

