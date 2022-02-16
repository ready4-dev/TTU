make_report_fl_nms_ls <- function(){
  report_fl_nms_ls <- ready4show::make_rmd_fl_nms_ls("Lngl_Mdls_HTML",
                                                     pdf_fl_nm_1L_chr = "Lngl_Mdls_PDF",
                                                     word_fl_nm_1L_chr = "Lngl_Mdls_Word")
  return(report_fl_nms_ls)
}
