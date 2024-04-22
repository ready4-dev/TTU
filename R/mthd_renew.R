#' 
#' Renew (update) values
#' @name renew-TTUProject
#' @description renew method applied to TTUProject
#' @param x An object of class TTUProject
#' @param new_val_xx New value (an output object of multiple potential types), Default: NULL
#' @param consent_1L_chr Consent (a character vector of length one), Default: ''
#' @param depnt_var_min_val_1L_dbl Dependent variable minimum value (a double vector of length one), Default: numeric(0)
#' @param fl_nm_1L_chr File name (a character vector of length one), Default: character(0)
#' @param paths_chr Paths (a character vector), Default: character(0)
#' @param type_1L_chr Type (a character vector of length one), Default: 'default'
#' @param y_Ready4useRepos PARAM_DESCRIPTION, Default: ready4use::Ready4useRepos()
#' @param what_1L_chr What (a character vector of length one), Default: 'utility'
#' @param ... Additional arguments
#' @return x (An object of class TTUProject)
#' @rdname renew-methods
#' @aliases renew,TTUProject-method
#' @export 
#' @importFrom ready4use Ready4useRepos
#' @importFrom ready4 renew
methods::setMethod("renew", "TTUProject", function (x, new_val_xx = NULL, consent_1L_chr = "", depnt_var_min_val_1L_dbl = numeric(0), 
    fl_nm_1L_chr = character(0), paths_chr = character(0), type_1L_chr = "default", 
    y_Ready4useRepos = ready4use::Ready4useRepos(), what_1L_chr = "utility", 
    ...) 
{
    if (what_1L_chr == "parameters") {
        if (type_1L_chr == "default") {
            x <- renewSlot(x, "b_SpecificParameters", SpecificConverter(a_ScorzProfile = x@a_ScorzProfile) %>% 
                metamorphose() %>% procureSlot("b_SpecificParameters"))
        }
        if (type_1L_chr == "range") {
            x <- renewSlot(x, "b_SpecificParameters@depnt_var_min_max_dbl", 
                new_val_xx)
        }
        if (type_1L_chr == "predictors_lup") {
            if (identical(new_val_xx, "use_renew_mthd")) {
                predictors_lup <- y_Ready4useRepos %>% ingest(fls_to_ingest_chr = c(fl_nm_1L_chr), 
                  metadata_1L_lgl = F)
            }
            else {
                predictors_lup <- new_val_xx
            }
            x <- renewSlot(x, "b_SpecificParameters@predictors_lup", 
                predictors_lup)
        }
        if (type_1L_chr == "predictors_vars") {
            x <- renewSlot(x, "b_SpecificParameters@candidate_predrs_chr", 
                new_val_xx)
        }
        if (type_1L_chr == "covariates") {
            x <- renewSlot(x, "b_SpecificParameters@candidate_covars_chr", 
                new_val_xx)
        }
        if (type_1L_chr == "descriptives") {
            x <- renewSlot(x, "b_SpecificParameters@descv_var_nms_chr", 
                new_val_xx)
        }
        if (type_1L_chr == "is_fake") {
            x <- renewSlot(x, "b_SpecificParameters@fake_1L_lgl", 
                new_val_xx)
        }
        if (type_1L_chr == "temporal") {
            x <- renewSlot(x, "b_SpecificParameters@msrmnt_date_var_nm_1L_chr", 
                new_val_xx)
        }
    }
    if (what_1L_chr == "project") {
        if (type_1L_chr == "default") {
            x <- renewSlot(x, "c_SpecificProject", SpecificModels(a_YouthvarsProfile = x@a_ScorzProfile@a_YouthvarsProfile, 
                b_SpecificParameters = x@b_SpecificParameters, 
                paths_chr = paths_chr))
            x <- ratifySlot(x, "c_SpecificProject")
            x <- renewSlot(x, "c_SpecificProject", authorSlot(x, 
                "c_SpecificProject", consent_1L_chr = consent_1L_chr, 
                what_1L_chr = "workspace"))
        }
        if (type_1L_chr == "dummys") {
            x <- renewSlot(x, "c_SpecificProject", renew(x@c_SpecificProject, 
                new_val_xx, what_1L_chr = type_1L_chr))
        }
    }
    if (what_1L_chr == "records") {
        if (type_1L_chr == "ds") {
            x <- renewSlot(x, "a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad@ds_tb", 
                new_val_xx)
        }
        if (type_1L_chr == "dict") {
            x <- renewSlot(x, "a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad@dictionary_r3", 
                new_val_xx)
        }
    }
    if (what_1L_chr == "reporting") {
        if (type_1L_chr == "default") {
            Y <- metamorphoseSlot(x, "c_SpecificProject")
            Y <- TTUSynopsis(a_Ready4showPaths = Y@a_Ready4showPaths, 
                b_SpecificResults = Y@b_SpecificResults, c_SpecificParameters = Y@c_SpecificParameters, 
                d_YouthvarsProfile = Y@d_YouthvarsProfile, rmd_fl_nms_ls = Y@rmd_fl_nms_ls)
            Y <- TTUReports(a_TTUSynopsis = Y)
            x <- renewSlot(x, "d_TTUReports", Y)
        }
        if (type_1L_chr == "abstract") {
            if (identical(new_val_xx, "use_renew_mthd")) {
                descs_ls <- x@d_TTUReports@a_TTUSynopsis@b_SpecificResults@a_SpecificShareable@shareable_outp_ls$results_ls$study_descs_ls
                x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis@abstract_args_ls", 
                  manufactureSlot(x, "d_TTUReports@a_TTUSynopsis", 
                    what_1L_chr = "abstract_args_ls", depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                    depnt_var_nms_chr = c(descs_ls$health_utl_nm_1L_chr, 
                      descs_ls$health_utl_long_nm_1L_chr)))
            }
            else {
                x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", 
                  procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% 
                    renewSlot("abstract_args_ls", new_val_xx))
            }
        }
        if (type_1L_chr == "authors") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("authors_r3", 
                new_val_xx))
        }
        if (type_1L_chr == "background") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("background_1L_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "changes") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("correspondences_r3", 
                new_val_xx))
        }
        if (type_1L_chr == "conflicts") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("coi_1L_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "conclusion") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("conclusion_1L_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "digits") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("digits_int", 
                new_val_xx))
        }
        if (type_1L_chr == "ethics") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("ethics_1L_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "formats") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("outp_formats_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "figures-body") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("figures_in_body_lgl", 
                new_val_xx))
        }
        if (type_1L_chr == "funding") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("funding_1L_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "institutes") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("institutes_r3", 
                new_val_xx))
        }
        if (type_1L_chr == "interval") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("interval_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "keywords") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("keywords_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "naming") {
            x <- enhanceSlot(x, "d_TTUReports@a_TTUSynopsis", 
                depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                depnt_var_nms_chr = new_val_xx, with_1L_chr = "results_ls")
        }
        if (type_1L_chr == "repos") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("e_Ready4useRepos", 
                new_val_xx))
        }
        if (type_1L_chr == "sample") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("sample_desc_1L_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "tables-body") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("tables_in_body_lgl", 
                new_val_xx))
        }
        if (type_1L_chr == "template-catalaogue") {
            x <- renewSlot(x, "d_TTUReports", procureSlot(x, 
                "d_TTUReports") %>% renewSlot("catalogue_tmpl_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "template-manuscript") {
            x <- renewSlot(x, "d_TTUReports", procureSlot(x, 
                "d_TTUReports") %>% renewSlot("manuscript_tmpl_chr", 
                new_val_xx))
        }
        if (type_1L_chr == "title") {
            x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, 
                "d_TTUReports@a_TTUSynopsis") %>% renewSlot("title_1L_chr", 
                new_val_xx))
        }
    }
    if (what_1L_chr == "results") {
        if (type_1L_chr == "covariates") {
            x <- renewSlot(x, "c_SpecificProject", renew(procureSlot(x, 
                "c_SpecificProject"), new_val_xx = new_val_xx, 
                type_1L_chr = "results", what_1L_chr = "prefd_covars"))
        }
        if (type_1L_chr == "models") {
            x <- renewSlot(x, "c_SpecificProject", renew(procureSlot(x, 
                "c_SpecificProject"), new_val_xx = new_val_xx, 
                type_1L_chr = "results", what_1L_chr = "prefd_mdls"))
        }
    }
    if (what_1L_chr == "utility") {
        x <- renewSlot(x, "a_ScorzProfile")
    }
    return(x)
})
