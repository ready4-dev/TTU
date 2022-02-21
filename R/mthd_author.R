#' 
#' Author and save files
#' @name author-TTUReports
#' @description author method applied to TTUReports
#' @param x An object of class TTUReports
#' @param depnt_var_desc_1L_chr Dependent variable description (a character vector of length one), Default: 'NA'
#' @param download_tmpl_1L_lgl Download template (a logical vector of length one), Default: T
#' @param fl_type_1L_chr File type (a character vector of length one), Default: '.eps'
#' @param timepoint_new_nms_chr Timepoint new names (a character vector), Default: 'NA'
#' @param type_1L_chr Type (a character vector of length one), Default: 'Report'
#' @param what_1L_chr What (a character vector of length one), Default: 'NA'
#' @return NULL
#' @rdname author-methods
#' @aliases author,TTUReports-method
#' @export 
#' @importFrom purrr map flatten_chr discard reduce map_chr pluck
#' @importFrom utils packageDescription
#' @importFrom stringr str_replace_all str_locate str_sub
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggsave
#' @importFrom ready4 author
methods::setMethod("author", "TTUReports", function (x, depnt_var_desc_1L_chr = NA_character_, download_tmpl_1L_lgl = T, 
    fl_type_1L_chr = ".eps", timepoint_new_nms_chr = NA_character_, 
    type_1L_chr = "Report", what_1L_chr = NA_character_) 
{
    if (type_1L_chr == "Report") {
        if (download_tmpl_1L_lgl) {
            authorData(x@a_TTUSynopsis, tmpl_url_1L_chr = ifelse(what_1L_chr == 
                "Catalogue", x@catalogue_tmpl_chr[1], x@manuscript_tmpl_chr[1]), 
                tmpl_version_1_L_chr = ifelse(what_1L_chr == 
                  "Catalogue", x@catalogue_tmpl_chr[2], x@manuscript_tmpl_chr[2]), 
                what_1L_chr = what_1L_chr)
        }
        if (what_1L_chr == "Catalogue") {
            x@a_TTUSynopsis@rmd_fl_nms_ls <- x@catalogue_fl_nms_ls
        }
        else {
            x@a_TTUSynopsis@rmd_fl_nms_ls <- x@manuscript_fl_nms_ls
        }
        if (what_1L_chr == "Catalogue") {
            author(x@a_TTUSynopsis, type_1L_chr = "Report", what_1L_chr = what_1L_chr)
        }
        else {
            authorReport(x@a_TTUSynopsis, what_1L_chr = what_1L_chr)
        }
    }
    else {
        dir_1L_chr <- paste0(x@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr, 
            "/", x@a_TTUSynopsis@a_Ready4showPaths@mkdn_data_dir_1L_chr, 
            "/", what_1L_chr)
        if (type_1L_chr == "Dependencies") {
            df <- data.frame(Package = c("youthvars", "scorz", 
                "specific", "TTU", ) %>% purrr::map(~{
                utils::packageDescription(.x) %>% c("Depends", 
                  "Imports")[] %>% purrr::map(~{
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
            saveRDS(df, paste0(dir_1L_chr, "/packages.RDS"))
        }
        if (type_1L_chr == "Plots") {
            composite_1_plt <- depictSlot(x, "a_TTUSynopsis", 
                depnt_var_desc_1L_chr = depnt_var_desc_1L_chr, 
                timepoint_old_nms_chr = procureSlot(x, "a_TTUSynopsis@d_YouthvarsProfile@timepoint_vals_chr"), 
                timepoint_new_nms_chr = timepoint_new_nms_chr, 
                what_1L_chr = "composite_mdl", write_1L_lgl = T)
            composite_2_plt <- depictSlot(x, slot_nm_1L_chr = "a_TTUSynopsis", 
                what_1L_chr = "composite_utl", write_1L_lgl = T)
            if (!is.na(what_1L_chr)) {
                ggplot2::ggsave(file = paste0(dir_1L_chr, "/fig1", 
                  fl_type_1L_chr), composite_2_plt)
                ggplot2::ggsave(file = paste0(dir_1L_chr, "/fig2", 
                  fl_type_1L_chr), composite_1_plt)
            }
        }
    }
})
