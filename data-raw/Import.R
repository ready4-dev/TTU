author_TTUReports <- function(x,
                              depnt_var_desc_1L_chr = NA_character_,
                              depnt_var_min_val_1L_dbl = numeric(0),
                              download_tmpl_1L_lgl = T,
                              fl_type_1L_chr = ".eps",
                              timepoint_new_nms_chr = NA_character_,
                              type_1L_chr = "Report",
                              what_1L_chr = NA_character_,
                              ...){
  if(type_1L_chr == "Report"){
    if(download_tmpl_1L_lgl){
      authorData(x@a_TTUSynopsis,
                 tmpl_url_1L_chr = ifelse(what_1L_chr == "Catalogue",
                                          x@catalogue_tmpl_chr[1],
                                          x@manuscript_tmpl_chr[1]),
                 tmpl_version_1L_chr = ifelse(what_1L_chr == "Catalogue",
                                               x@catalogue_tmpl_chr[2],
                                               x@manuscript_tmpl_chr[2]),
                 what_1L_chr = what_1L_chr)
    }
    if(what_1L_chr == "Catalogue"){
      x@a_TTUSynopsis@rmd_fl_nms_ls <- x@catalogue_fl_nms_ls
    }else{
      x@a_TTUSynopsis@rmd_fl_nms_ls <- x@manuscript_fl_nms_ls
    }
    if(what_1L_chr == "Catalogue"){
      author(x@a_TTUSynopsis,
             type_1L_chr = "Report",
             what_1L_chr = what_1L_chr)
    }else{
      author_SpecificSynopsis(x@a_TTUSynopsis, # authorReport when Exported
                              what_1L_chr = what_1L_chr,
                              ...)
    }
  }else{
    dir_1L_chr <- paste0(x@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr,
                         "/",
                         x@a_TTUSynopsis@a_Ready4showPaths@mkdn_data_dir_1L_chr,
                         "/",
                         what_1L_chr)
    if(type_1L_chr == "Dependencies"){
      df <- data.frame(Package = c("youthvars","scorz","specific","TTU") %>%
                         purrr::map(~ {
                           utils::packageDescription(.x) %>%
                             `[`(c("Depends", "Imports")) %>%
                             purrr::map(~{
                               if(is.null(.x)){
                                 character(0)
                               }else{
                                 .x %>%
                                   strsplit(",\\n") %>%
                                   purrr::flatten_chr() %>%
                                   purrr::map(~strsplit(.x,", ") %>%
                                                purrr::flatten_chr()) %>%
                                   purrr::flatten_chr() %>% sort() %>%
                                   purrr::discard(~startsWith(.x,"R "))
                               }
                             }) %>%
                             purrr::flatten_chr() %>%
                             unique() %>%
                             sort()
                         }) %>%
                         purrr::reduce(~c(.x,.y)) %>%
                         purrr::map_chr(~{
                           updated_1L_chr <- stringr::str_replace_all(.x,"\\n"," ")
                           problem_idx_1L_chr <- stringr::str_locate(updated_1L_chr," ")[1,1] %>%
                             unname()
                           if(!is.na(problem_idx_1L_chr))
                             updated_1L_chr <- updated_1L_chr %>%
                             stringr::str_sub(end = problem_idx_1L_chr-1)
                           updated_1L_chr %>% trimws(which = "left")
                         }) %>% unique() %>% sort())
      df <- df %>%
        dplyr::mutate(Version = Package %>%
                        purrr::map_chr(~utils::packageDescription(.x) %>%
                                         purrr::pluck("Version")),
                      Citation = Package %>%
                        purrr::map_chr(~get_pkg_citation(.x)))
      saveRDS(df,
              paste0(dir_1L_chr,
                     "/packages.RDS"))
      
    }
    if(type_1L_chr == "Plots"){
      composite_1_plt <- depict_SpecificSynopsis(x@a_TTUSynopsis,#depictSlot(x,"a_TTUSynopsis", when exported
                                                 depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                                                 depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                 timepoint_old_nms_chr = procureSlot(x,
                                                                                     "a_TTUSynopsis@d_YouthvarsProfile@timepoint_vals_chr"),
                                                 timepoint_new_nms_chr = timepoint_new_nms_chr,
                                                 what_1L_chr = "composite_mdl",
                                                 write_1L_lgl = T)
      composite_2_plt <- depict_SpecificSynopsis(x@a_TTUSynopsis,#depictSlot(x,"a_TTUSynopsis", when exported
                                                 what_1L_chr = "composite_utl",
                                                 write_1L_lgl = T)
      if(!is.na(what_1L_chr)){
        ggplot2::ggsave(file = paste0(dir_1L_chr,
                                      "/fig1",
                                      fl_type_1L_chr),
                        composite_2_plt)
        ggplot2::ggsave(file = paste0(dir_1L_chr,
                                      "/fig2",
                                      fl_type_1L_chr),
                        composite_1_plt)
        
      }
      
    }
  }
}
