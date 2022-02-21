#' Get package citation
#' @description get_pkg_citation() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get package citation. Function argument pkg_nm_1L_chr specifies the where to look for the required object. The function returns Citation (a character vector of length one).
#' @param pkg_nm_1L_chr Package name (a character vector of length one)
#' @return Citation (a character vector of length one)
#' @rdname get_pkg_citation
#' @export 
#' @keywords internal
get_pkg_citation <- function (pkg_nm_1L_chr) 
{
    citation_chr <- suppressWarnings(citation(pkg_nm_1L_chr)) %>% 
        capture.output()
    start_idx_1L_int <- 4
    end_idx_1L_int <- which(citation_chr == "")[which(which(citation_chr == 
        "") > start_idx_1L_int)[1]] - 1
    citation_1L_chr <- citation_chr[start_idx_1L_int:end_idx_1L_int] %>% 
        paste0(collapse = "")
    return(citation_1L_chr)
}
