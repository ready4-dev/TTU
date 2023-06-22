renew_TTUProject <- function(x,
                             new_val_xx = NULL,
                             consent_1L_chr = "",
                             fl_nm_1L_chr = character(0),
                             paths_chr = character(0),
                             predictors_lup = NULL,
                             type_1L_chr = "default",
                             y_Ready4useRepos = ready4use::Ready4useRepos(),
                             what_1L_chr = "utility",
                             ...){
  if(what_1L_chr == "parameters"){
    if(type_1L_chr=="default"){
      x <- renewSlot(x, "b_SpecificParameters", SpecificConverter(a_ScorzProfile = x@a_ScorzProfile) %>%
                       metamorphose() %>%
                       procureSlot("b_SpecificParameters"))
    }
    if(type_1L_chr == "range"){
      x <- renewSlot(x, "b_SpecificParameters@depnt_var_min_max_dbl", new_val_xx)
    }
    if(type_1L_chr=="predictors_lup"){
      if(is.null(predictors_lup)){
        predictors_lup <- Ready4useRepos(dv_nm_1L_chr = "TTU",
                                         dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0",
                                         dv_server_1L_chr = "dataverse.harvard.edu") %>%
          ingest(fls_to_ingest_chr = c(fl_nm_1L_chr),
                 metadata_1L_lgl = F)
      }
      x <- renewSlot(x, "b_SpecificParameters@predictors_lup", predictors_lup)
    }
    if(type_1L_chr == "predictors_vars"){
      x <- renewSlot(x, "b_SpecificParameters@candidate_predrs_chr", new_val_xx)
    }
    if(type_1L_chr == "covariates"){
      x <- renewSlot(x, "b_SpecificParameters@candidate_covars_chr", new_val_xx)
    }
    if(type_1L_chr == "descriptives"){
      x <- renewSlot(x, "b_SpecificParameters@descv_var_nms_chr", new_val_xx)
    }
    if(type_1L_chr == "is_fake"){
      x <- renewSlot(x, "b_SpecificParameters@fake_1L_lgl", new_val_xx)
    }
    if(type_1L_chr == "temporal"){
      x <- renewSlot(x, "b_SpecificParameters@msrmnt_date_var_nm_1L_chr", new_val_xx)
    }
  }
  if(what_1L_chr == "project"){
    x <- renewSlot(x, "c_SpecificProject", SpecificModels(a_YouthvarsProfile = x@a_ScorzProfile@a_YouthvarsProfile,
                                                          b_SpecificParameters = x@b_SpecificParameters,
                                                          paths_chr = paths_chr))
    x <- ratifySlot(x, "c_SpecificProject")
    x <- renewSlot(x, "c_SpecificProject",
                   authorSlot(x, "c_SpecificProject", consent_1L_chr = consent_1L_chr, what_1L_chr = "workspace"))
  }
  if(what_1L_chr == "reporting"){
    if(type_1L_chr=="default"){
      x <- renewSlot(x, "d_TTUReports",
                     {
                       Y <- metamorphoseSlot(x, "c_SpecificProject")
                       Y <- TTUSynopsis(a_Ready4showPaths = Y@a_Ready4showPaths,
                                        b_SpecificResults = Y@b_SpecificResults,
                                        c_SpecificParameters = Y@c_SpecificParameters,
                                        d_YouthvarsProfile = Y@d_YouthvarsProfile,
                                        rmd_fl_nms_ls = Y@rmd_fl_nms_ls)
                       Y <- TTUReports(a_TTUSynopsis = Y)
                       Y
                     })
    }
    if(type_1L_chr == "authors"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis",
                     procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>%
                       renewSlot("authors_r3", new_val_xx))
    }
    if(type_1L_chr == "changes"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis",
                     procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>%
                       renewSlot("correspondences_r3", new_val_xx))
    }
    if(type_1L_chr == "digits"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis",
                     procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>%
                       renewSlot("digits_int", new_val_xx))
    }
    if(type_1L_chr == "formats"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis",
                     procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>%
                       renewSlot("outp_formats_chr", new_val_xx))
    }
    if(type_1L_chr == "institutes"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis",
                     procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>%
                       renewSlot("institutes_r3", new_val_xx))
    }
    if(type_1L_chr == "repos"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis",
                     procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>%
                       renewSlot("e_Ready4useRepos", new_val_xx))
    }
    if(type_1L_chr == "title"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis",
                     procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>%
                       renewSlot("title_1L_chr", new_val_xx))
    }
  }
  if(what_1L_chr == "results"){
    if(type_1L_chr =="covariates"){
      x <- renewSlot(x, "c_SpecificProject",
                     renew(procureSlot(x, "c_SpecificProject"),
                           new_val_xx = new_val_xx,
                           type_1L_chr = "results",
                           what_1L_chr = "prefd_covars"))
    }
    if(type_1L_chr == "models"){
      x <- renewSlot(x, "c_SpecificProject",
                     renew(procureSlot(x, "c_SpecificProject"),
                           new_val_xx = new_val_xx,
                           type_1L_chr = "results",
                           what_1L_chr = "prefd_mdls"))
    }

  }
  if(what_1L_chr == "utility"){
    x <- renewSlot(x, "a_ScorzProfile")
  }
  return(x)
}
