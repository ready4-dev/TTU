#' Make report file names list
#' @description make_report_fl_nms_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make report file names list. The function returns Report file names (a list).

#' @return Report file names (a list)
#' @rdname make_report_fl_nms_ls
#' @export 
#' @importFrom ready4show make_rmd_fl_nms_ls
#' @keywords internal
make_report_fl_nms_ls <- function () 
{
    report_fl_nms_ls <- ready4show::make_rmd_fl_nms_ls("Lngl_Mdls_HTML", 
        pdf_fl_nm_1L_chr = "Lngl_Mdls_PDF", word_fl_nm_1L_chr = "Lngl_Mdls_Word")
    return(report_fl_nms_ls)
}
