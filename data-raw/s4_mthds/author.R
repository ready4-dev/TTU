author_TTUReports <- function(x,
                              download_tmpl_1L_lgl = T,
                              what_1L_chr = "Catalogue"){
  if(download_tmpl_1L_lgl){
    authorData(x@a_SpecificSynopsis,
               tmpl_url_1L_chr = ifelse(what_1L_chr == "Catalogue",
                                        x@catalogue_tmpl_chr[1],
                                        x@manuscript_tmpl_chr[1]),
               tmpl_version_1_L_chr = ifelse(what_1L_chr == "Catalogue",
                                             x@catalogue_tmpl_chr[2],
                                             x@manuscript_tmpl_chr[2]),
               what_1L_chr = what_1L_chr)
  }
  if(what_1L_chr == "Catalogue"){
    x@a_SpecificSynopsis@rmd_fl_nms_ls <- x@catalogue_fl_nms_ls
  }else{
    x@a_SpecificSynopsis@rmd_fl_nms_ls <- x@manuscript_fl_nms_ls
  }
  if(what_1L_chr == "Catalogue"){
    author(x@a_SpecificSynopsis,
           type_1L_chr = "Report",
           what_1L_chr = what_1L_chr)
  }else{
    authorReport(x@a_SpecificSynopsis,
                 what_1L_chr = what_1L_chr)
  }
}
