get_pkg_citation <- function(pkg_nm_1L_chr){
  citation_chr <- suppressWarnings(citation(pkg_nm_1L_chr)) %>% capture.output()
  start_idx_1L_int <- 4
  end_idx_1L_int <- which(citation_chr== "")[which(which(citation_chr== "")>start_idx_1L_int)[1]]-1
  citation_1L_chr<- citation_chr[start_idx_1L_int:end_idx_1L_int] %>% paste0(collapse = "")
  return(citation_1L_chr)
}
