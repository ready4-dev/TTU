share_TTUProject <- function(x,
                             formats_chr = c(".docx",".pdf", ".tex"),
                             types_chr = "auto",
                             what_chr = c("catalogue", "models")){
  if("manuscript" %in% what_chr | "supplement" %in% what_chr){
    purrr::walk(paste0("Manuscript_", Hmisc::capitalize(types_chr)),
                 ~ {
                   x <- x@d_TTUReports@a_TTUSynopsis
                   what_1L_chr <- .x

                   files_chr <- list.files(paste0(x@a_Ready4showPaths@outp_data_dir_1L_chr,
                                       "/",
                                       x@a_Ready4showPaths@reports_dir_1L_chr,
                                       "/",
                                       what_1L_chr))
                   files_chr <- files_chr[files_chr %>% purrr::map_lgl(~{
                     string_1L_chr <- .x
                     any(formats_chr %>% purrr::map_lgl(~endsWith(string_1L_chr,.x)))
                     })]
                   files_chr <- files_chr[files_chr %>% purrr::map_lgl(~{
                     string_1L_chr <- .x
                     any(Hmisc::capitalize(what_chr) %>% purrr::map_lgl(~startsWith(string_1L_chr,.x)))
                   })]
                     ms_nm_1L_chr <- "Manuscript"
                     idx_1L_int <- ready4::write_fls_to_dv(paste0(paste0(x@a_Ready4showPaths@outp_data_dir_1L_chr,
                                                                  "/",
                                                                  x@a_Ready4showPaths@reports_dir_1L_chr,
                                                                  "/",
                                                                  what_1L_chr,
                                                                  "/"),
                                                                  files_chr),
                                                           descriptions_chr = files_chr %>% purrr::map_chr(~paste0("Scientific summary of utility mapping study",
                                                                                                                   " ",
                                                                                                                   ifelse(startsWith(.x,"Supplement"),"(Supplement) ",""),
                                                                                                                   ifelse(endsWith(.x,".tex"),"(LaTeX) ",""),
                                                                                                                   ifelse(what_1L_chr == "Manuscript_Auto", " (algorithm generated)",""))),
                                                           ds_url_1L_chr = x@e_Ready4useRepos@dv_ds_nm_1L_chr)
                   Sys.sleep(5L)
                 }
    )
    }
    if("catalogue" %in% what_chr){
      shareSlot(x, "d_TTUReports@a_TTUSynopsis", type_1L_chr = "Report", what_1L_chr = "Catalogue")
    }
    if("models" %in% what_chr){
      shareSlot(A, "d_TTUReports@a_TTUSynopsis", type_1L_chr = "Models", what_1L_chr = "ingredients")
    }

  return(x)



}
