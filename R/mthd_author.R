#' 
#' Author and save files
#' @name author-TTUReports
#' @description author method applied to TTUReports
#' @param x An object of class TTUReports
#' @param args_ls Arguments (a list), Default: NULL
#' @param consent_1L_chr Consent (a character vector of length one), Default: ''
#' @param depnt_var_desc_1L_chr Dependent variable description (a character vector of length one), Default: 'NA'
#' @param depnt_var_min_val_1L_dbl Dependent variable minimum value (a double vector of length one), Default: numeric(0)
#' @param download_tmpl_1L_lgl Download template (a logical vector of length one), Default: T
#' @param fl_type_1L_chr File type (a character vector of length one), Default: '.eps'
#' @param timepoint_new_nms_chr Timepoint new names (a character vector), Default: 'NA'
#' @param type_1L_chr Type (a character vector of length one), Default: 'Report'
#' @param what_1L_chr What (a character vector of length one), Default: 'NA'
#' @param ... Additional arguments
#' @return NULL
#' @rdname author-methods
#' @aliases author,TTUReports-method
#' @export 
#' @importFrom purrr map flatten_chr discard reduce map_chr pluck
#' @importFrom utils packageDescription
#' @importFrom stringr str_replace_all str_locate str_sub
#' @importFrom dplyr mutate
#' @importFrom ready4 write_with_consent make_list_phrase author
#' @importFrom ggplot2 ggsave
methods::setMethod("author", "TTUReports", function (x, args_ls = NULL, consent_1L_chr = "", depnt_var_desc_1L_chr = NA_character_, 
    depnt_var_min_val_1L_dbl = numeric(0), download_tmpl_1L_lgl = T, 
    fl_type_1L_chr = ".eps", timepoint_new_nms_chr = NA_character_, 
    type_1L_chr = "Report", what_1L_chr = NA_character_, ...) 
{
    if (type_1L_chr == "Report") {
        if (download_tmpl_1L_lgl) {
            authorData(x@a_TTUSynopsis, consent_1L_chr = consent_1L_chr, 
                tmpl_url_1L_chr = ifelse(what_1L_chr == "Catalogue", 
                  x@catalogue_tmpl_chr[1], x@manuscript_tmpl_chr[1]), 
                tmpl_version_1L_chr = ifelse(what_1L_chr == "Catalogue", 
                  x@catalogue_tmpl_chr[2], x@manuscript_tmpl_chr[2]), 
                what_1L_chr = what_1L_chr)
        }
        if (what_1L_chr == "Catalogue") {
            x@a_TTUSynopsis@rmd_fl_nms_ls <- x@catalogue_fl_nms_ls
        }
        else {
            x@a_TTUSynopsis@rmd_fl_nms_ls <- x@manuscript_fl_nms_ls
        }
        if (what_1L_chr == "Catalogue") {
            author(x@a_TTUSynopsis, args_ls = args_ls, consent_1L_chr = consent_1L_chr, 
                type_1L_chr = type_1L_chr, what_1L_chr = what_1L_chr)
        }
        else {
            authorReport(x@a_TTUSynopsis, args_ls = args_ls, 
                consent_1L_chr = consent_1L_chr, type_1L_chr = type_1L_chr, 
                what_1L_chr = what_1L_chr)
        }
    }
    else {
        dir_1L_chr <- paste0(x@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr, 
            "/", x@a_TTUSynopsis@a_Ready4showPaths@mkdn_data_dir_1L_chr, 
            "/", what_1L_chr)
        if (type_1L_chr == "Dependencies") {
            df <- data.frame(Package = c("youthvars", "scorz", 
                "specific", "TTU") %>% purrr::map(~{
                desc_ls <- utils::packageDescription(.x)
                desc_ls[c("Depends", "Imports")] %>% purrr::map(~{
                  if (is.null(.x)) {
                    character(0)
                  }
                  else {
                    .x %>% strsplit(",\\n") %>% purrr::flatten_chr() %>% 
                      purrr::map(~strsplit(.x, ", ") %>% purrr::flatten_chr()) %>% 
                      purrr::flatten_chr() %>% sort() %>% purrr::discard(~startsWith(.x, 
                      "R "))
                  }
                }) %>% purrr::flatten_chr() %>% unique() %>% 
                  sort()
            }) %>% purrr::reduce(~c(.x, .y)) %>% purrr::map_chr(~{
                updated_1L_chr <- stringr::str_replace_all(.x, 
                  "\\n", " ")
                problem_idx_1L_chr <- stringr::str_locate(updated_1L_chr, 
                  " ")[1, 1] %>% unname()
                if (!is.na(problem_idx_1L_chr)) 
                  updated_1L_chr <- updated_1L_chr %>% stringr::str_sub(end = problem_idx_1L_chr - 
                    1)
                updated_1L_chr %>% trimws(which = "left")
            }) %>% unique() %>% sort())
            df <- df %>% dplyr::mutate(Version = Package %>% 
                purrr::map_chr(~utils::packageDescription(.x) %>% 
                  purrr::pluck("Version")), Citation = Package %>% 
                purrr::map_chr(~get_pkg_citation(.x)))
            ready4::write_with_consent(consented_fn = saveRDS, 
                prompt_1L_chr = paste0("Do you confirm that you want to write the file ", 
                  "packages.RDS", " to ", dir_1L_chr, "?"), consent_1L_chr = consent_1L_chr, 
                consented_args_ls = list(object = df, file = paste0(dir_1L_chr, 
                  "/packages.RDS")), consented_msg_1L_chr = paste0("File ", 
                  "packages.RDS", " has been written to ", dir_1L_chr, 
                  "."), declined_msg_1L_chr = "Write request cancelled - no new files have been written.")
        }
        if (type_1L_chr == "Plots") {
            composite_1_plt <- depict(x@a_TTUSynopsis, consent_1L_chr = consent_1L_chr, 
                depnt_var_desc_1L_chr = depnt_var_desc_1L_chr, 
                depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                timepoint_old_nms_chr = procureSlot(x, "a_TTUSynopsis@d_YouthvarsProfile@timepoint_vals_chr"), 
                timepoint_new_nms_chr = timepoint_new_nms_chr, 
                what_1L_chr = "composite_mdl", write_1L_lgl = T)
            composite_2_plt <- depict(x@a_TTUSynopsis, consent_1L_chr = consent_1L_chr, 
                what_1L_chr = "composite_utl", write_1L_lgl = T)
            if (!is.na(what_1L_chr)) {
                consented_fn <- function(composite_1_plt, composite_2_plt, 
                  dir_1L_chr, fl_type_1L_chr) {
                  ggplot2::ggsave(file = paste0(dir_1L_chr, "/fig1", 
                    fl_type_1L_chr), composite_2_plt)
                  ggplot2::ggsave(file = paste0(dir_1L_chr, "/fig2", 
                    fl_type_1L_chr), composite_1_plt)
                }
                ready4::write_with_consent(consented_fn = consented_fn, 
                  prompt_1L_chr = paste0("Do you confirm that you want to write the files ", 
                    ready4::make_list_phrase(paste0("fig", 1:2, 
                      fl_type_1L_chr)), " to ", dir_1L_chr, "?"), 
                  consent_1L_chr = consent_1L_chr, consented_args_ls = list(composite_1_plt = composite_1_plt, 
                    composite_2_plt = composite_2_plt, dir_1L_chr = dir_1L_chr, 
                    fl_type_1L_chr = fl_type_1L_chr), consented_msg_1L_chr = paste0("Files ", 
                    ready4::make_list_phrase(paste0("fig", 1:2, 
                      fl_type_1L_chr)), " have been written to ", 
                    dir_1L_chr, "."), declined_msg_1L_chr = "Write request cancelled - no new files have been written.")
            }
        }
    }
})
#' 
#' Author and save files
#' @name author-TTUProject
#' @description author method applied to TTUProject
#' @param x An object of class TTUProject
#' @param consent_1L_chr Consent (a character vector of length one), Default: ''
#' @param custom_args_ls Custom arguments (a list), Default: NULL
#' @param custom_fn Custom (a function), Default: NULL
#' @param depnt_var_min_val_1L_dbl Dependent variable minimum value (a double vector of length one), Default: numeric(0)
#' @param digits_1L_int Digits (an integer vector of length one), Default: 2
#' @param download_tmpl_1L_lgl Download template (a logical vector of length one), Default: T
#' @param fl_nm_1L_chr File name (a character vector of length one), Default: 'TTUProject'
#' @param items_as_domains_1L_lgl Items as domains (a logical vector of length one), Default: F
#' @param supplement_fl_nm_1L_chr Supplement file name (a character vector of length one), Default: 'TA_PDF'
#' @param timepoint_new_nms_chr Timepoint new names (a character vector), Default: 'NA'
#' @param type_1L_chr Type (a character vector of length one), Default: 'auto'
#' @param what_1L_chr What (a character vector of length one), Default: 'default'
#' @param ... Additional arguments
#' @return x (An object of class TTUProject)
#' @rdname author-methods
#' @aliases author,TTUProject-method
#' @export 
#' @importFrom ready4 write_with_consent author
#' @importFrom Hmisc capitalize
#' @importFrom R.utils copyDirectory
#' @importFrom ready4show make_rmd_fl_nms_ls
#' @importFrom methods callNextMethod
methods::setMethod("author", "TTUProject", function (x, consent_1L_chr = "", custom_args_ls = NULL, custom_fn = NULL, 
    depnt_var_min_val_1L_dbl = numeric(0), digits_1L_int = 2L, 
    download_tmpl_1L_lgl = T, fl_nm_1L_chr = "TTUProject", items_as_domains_1L_lgl = F, 
    supplement_fl_nm_1L_chr = "TA_PDF", timepoint_new_nms_chr = NA_character_, 
    type_1L_chr = "auto", what_1L_chr = "default", ...) 
{
    if (what_1L_chr %in% c("catalogue", "Catalogue", "dependencies", 
        "Dependencies", "descriptives", "Descriptives", "manuscript", 
        "Manuscript", "models", "Models", "plots", "Plots", "purge", 
        "Purge", "self", "Self", "supplement", "Supplement")) {
        if (what_1L_chr %in% c("self", "Self")) {
            to_1L_chr <- paste0(x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr, 
                "/", fl_nm_1L_chr, ".RDS")
            ready4::write_with_consent(consented_fn = saveRDS, 
                prompt_1L_chr = paste0("Do you confirm that you want to write a copy of this TTUProject module to ", 
                  to_1L_chr, "?"), consent_1L_chr = consent_1L_chr, 
                consented_args_ls = list(object = x, file = to_1L_chr), 
                consented_msg_1L_chr = paste0("A copy of this TTUProject module has been written to ", 
                  to_1L_chr, "."), declined_msg_1L_chr = "Write request cancelled - no new file has been written.")
        }
        if (what_1L_chr %in% c("catalogue", "Catalogue")) {
            authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, 
                download_tmpl_1L_lgl = download_tmpl_1L_lgl, 
                what_1L_chr = Hmisc::capitalize(what_1L_chr))
        }
        if (what_1L_chr %in% c("descriptives", "Descriptives")) {
            if (items_as_domains_1L_lgl == T) {
                x_labels_chr <- manufacture(x@a_ScorzProfile, 
                  what_1L_chr = "domains", custom_args_ls = list(string = x@b_SpecificParameters@itm_labels_chr), 
                  custom_fn = Hmisc::capitalize)
            }
            else {
                x_labels_chr <- manufacture(x@a_ScorzProfile, 
                  what_1L_chr = "domains")
            }
            x <- renewSlot(x, "c_SpecificProject", authorSlot(x, 
                "c_SpecificProject", consent_1L_chr = consent_1L_chr, 
                digits_1L_int = digits_1L_int, what_1L_chr = tolower(what_1L_chr), 
                x_labels_chr = x_labels_chr))
        }
        if (what_1L_chr %in% c("manuscript", "Manuscript")) {
            if (type_1L_chr == "auto") {
                x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr <- paste0(Hmisc::capitalize(what_1L_chr), 
                  "_", Hmisc::capitalize(type_1L_chr))
                authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, 
                  download_tmpl_1L_lgl = download_tmpl_1L_lgl, 
                  type_1L_chr = "Report", what_1L_chr = x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr)
            }
            if (type_1L_chr == "copy") {
                from_1L_chr <- paste0(x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr, 
                  "/", x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@mkdn_data_dir_1L_chr, 
                  "/Manuscript_Auto")
                to_1L_chr <- paste0(x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr, 
                  "/", x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@mkdn_data_dir_1L_chr, 
                  "/Manuscript_Submission")
                ready4::write_with_consent(consented_fn = R.utils::copyDirectory, 
                  prompt_1L_chr = paste0("Do you confirm that you want to copy the directory ", 
                    from_1L_chr, " (and all its contents) to ", 
                    to_1L_chr, "?"), consent_1L_chr = consent_1L_chr, 
                  consented_args_ls = list(from = from_1L_chr, 
                    to = to_1L_chr), consented_msg_1L_chr = paste0("The directory ", 
                    from_1L_chr, " has been copied to ", to_1L_chr, 
                    "."), declined_msg_1L_chr = "Write request cancelled - no new directory copy has been written.")
            }
            if (type_1L_chr %in% c("dependencies", "Dependencies")) {
                author(x@d_TTUReports, consent_1L_chr = consent_1L_chr, 
                  type_1L_chr = Hmisc::capitalize(type_1L_chr), 
                  what_1L_chr = "Manuscript_Submission")
            }
            if (type_1L_chr %in% c("plots", "Plots")) {
                authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, 
                  depnt_var_desc_1L_chr = x@d_TTUReports@a_TTUSynopsis@b_SpecificResults@a_SpecificShareable@shareable_outp_ls$results_ls$study_descs_ls$health_utl_nm_1L_chr, 
                  depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                  timepoint_new_nms_chr = timepoint_new_nms_chr, 
                  type_1L_chr = Hmisc::capitalize(type_1L_chr), 
                  what_1L_chr = "Manuscript_Submission")
            }
            if (type_1L_chr == "submission") {
                x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr <- paste0(Hmisc::capitalize(what_1L_chr), 
                  "_", Hmisc::capitalize(type_1L_chr))
                authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, 
                  download_tmpl_1L_lgl = download_tmpl_1L_lgl, 
                  type_1L_chr = "Report", what_1L_chr = x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr)
            }
        }
        if (what_1L_chr %in% c("models", "Models")) {
            x <- renewSlot(x, "c_SpecificProject", authorData(procureSlot(x, 
                "c_SpecificProject"), consent_1L_chr = consent_1L_chr))
        }
        if (what_1L_chr %in% c("purge")) {
            authorSlot(x, "c_SpecificProject", type_1L_chr = "purge_write", 
                consent_1L_chr = consent_1L_chr)
        }
        if (what_1L_chr %in% c("plots", "Plots")) {
            authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, 
                depnt_var_desc_1L_chr = x@d_TTUReports@a_TTUSynopsis@b_SpecificResults@a_SpecificShareable@shareable_outp_ls$results_ls$study_descs_ls$health_utl_nm_1L_chr, 
                type_1L_chr = Hmisc::capitalize(what_1L_chr))
        }
        if (what_1L_chr %in% c("supplement", "Supplement")) {
            x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr <- paste0("Manuscript_", 
                Hmisc::capitalize(type_1L_chr))
            authorReport(procureSlot(x, "d_TTUReports") %>% renewSlot("a_TTUSynopsis@rmd_fl_nms_ls", 
                ready4show::make_rmd_fl_nms_ls(pdf_fl_nm_1L_chr = supplement_fl_nm_1L_chr)) %>% 
                renewSlot("a_TTUSynopsis@outp_formats_chr", rep(x@d_TTUReports@a_TTUSynopsis@outp_formats_chr[2], 
                  2)) %>% procureSlot("a_TTUSynopsis"), consent_1L_chr = consent_1L_chr, 
                fl_nm_1L_chr = "Supplement", what_1L_chr = x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr)
        }
    }
    else {
        x <- methods::callNextMethod()
    }
    return(x)
})
