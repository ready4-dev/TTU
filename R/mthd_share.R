#' 
#' Share data via an online repository
#' @name share-TTUProject
#' @description share method applied to TTUProject
#' @param x An object of class TTUProject
#' @param formats_chr Formats (a character vector), Default: c(".docx", ".pdf", ".tex")
#' @param types_chr Types (a character vector), Default: 'auto'
#' @param what_chr What (a character vector), Default: c("catalogue", "instrument", "models")
#' @return x (An object of class TTUProject)
#' @rdname share-methods
#' @aliases share,TTUProject-method
#' @export 
#' @importFrom purrr walk map_lgl map_chr
#' @importFrom Hmisc capitalize
#' @importFrom ready4 write_fls_to_dv share
#' @importFrom dplyr filter
methods::setMethod("share", "TTUProject", function (x, formats_chr = c(".docx", ".pdf", ".tex"), types_chr = "auto", 
    what_chr = c("catalogue", "instrument", "models")) 
{
    if ("manuscript" %in% what_chr | "supplement" %in% what_chr) {
        purrr::walk(paste0("Manuscript_", Hmisc::capitalize(types_chr)), 
            ~{
                x <- x@d_TTUReports@a_TTUSynopsis
                what_1L_chr <- .x
                files_chr <- list.files(paste0(x@a_Ready4showPaths@outp_data_dir_1L_chr, 
                  "/", x@a_Ready4showPaths@reports_dir_1L_chr, 
                  "/", what_1L_chr))
                files_chr <- files_chr[files_chr %>% purrr::map_lgl(~{
                  string_1L_chr <- .x
                  any(formats_chr %>% purrr::map_lgl(~endsWith(string_1L_chr, 
                    .x)))
                })]
                files_chr <- files_chr[files_chr %>% purrr::map_lgl(~{
                  string_1L_chr <- .x
                  any(Hmisc::capitalize(what_chr) %>% purrr::map_lgl(~startsWith(string_1L_chr, 
                    .x)))
                })]
                ms_nm_1L_chr <- "Manuscript"
                idx_1L_int <- ready4::write_fls_to_dv(paste0(paste0(x@a_Ready4showPaths@outp_data_dir_1L_chr, 
                  "/", x@a_Ready4showPaths@reports_dir_1L_chr, 
                  "/", what_1L_chr, "/"), files_chr), descriptions_chr = files_chr %>% 
                  purrr::map_chr(~paste0("Scientific summary of utility mapping study", 
                    " ", ifelse(startsWith(.x, "Supplement"), 
                      "(Supplement) ", ""), ifelse(endsWith(.x, 
                      ".tex"), "(LaTeX) ", ""), ifelse(what_1L_chr == 
                      "Manuscript_Auto", " (algorithm generated)", 
                      ""))), ds_url_1L_chr = x@e_Ready4useRepos@dv_ds_nm_1L_chr)
                Sys.sleep(5L)
            })
    }
    if ("catalogue" %in% what_chr) {
        shareSlot(x, "d_TTUReports@a_TTUSynopsis", type_1L_chr = "Report", 
            what_1L_chr = "Catalogue")
    }
    if ("instrument" %in% what_chr) {
        descs_ls <- x@d_TTUReports@a_TTUSynopsis@b_SpecificResults@a_SpecificShareable@shareable_outp_ls$results_ls$study_descs_ls
        Y <- x
        Y@a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad@ds_tb <- Y@a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad@ds_tb %>% 
            dplyr::filter(F)
        instrument_ls <- list(X_ScorzProfile = Y@a_ScorzProfile, 
            depnt_var_nms_chr = c(descs_ls$health_utl_nm_1L_chr, 
                descs_ls$health_utl_long_nm_1L_chr))
        Y <- share(x@d_TTUReports@a_TTUSynopsis@e_Ready4useRepos, 
            description_1L_chr = "R list object with details of the utility instrument used in this study.", 
            obj_to_share_xx = instrument_ls, fl_nm_1L_chr = "instrument")
    }
    if ("models" %in% what_chr) {
        shareSlot(A, "d_TTUReports@a_TTUSynopsis", type_1L_chr = "Models", 
            what_1L_chr = "ingredients")
    }
    return(x)
})
