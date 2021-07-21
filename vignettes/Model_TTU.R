## ----message=FALSE, warning=FALSE---------------------------------------------
library(TTU)

## -----------------------------------------------------------------------------
seed_1L_int = 12345
set.seed(seed_1L_int)

## ----warning = F--------------------------------------------------------------
path_to_write_to_1L_chr <- "Output" 
dir.create(path_to_write_to_1L_chr)

## ----message=FALSE------------------------------------------------------------
data_tb <- readRDS(url("https://dataverse.harvard.edu/api/access/datafile/4750597"))

## ----inputds, echo=F, eval = knitr::is_html_output(), results='asis'----------
data_tb %>%
    head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Input dataset",
                          mkdn_tbl_ref_1L_chr = "tab:inputds")

## -----------------------------------------------------------------------------
predictors_lup <- TTU_predictors_lup(make_pt_TTU_predictors_lup(short_name_chr = c("k10_int","psych_well_int"),
                                              long_name_chr = c("Kessler Psychological Distress - 10 Item Total Score",
                                                   "Overall Wellbeing Measure (Winefield et al. 2012)"),
                                              min_val_dbl = c(10,18),
                                              max_val_dbl = c(50,90),
                                              class_chr = "integer",
                                              increment_dbl = 1,
                                              class_fn_chr = "as.integer",
                                              mdl_scaling_dbl = 0.01,
                                              covariate_lgl = F))

