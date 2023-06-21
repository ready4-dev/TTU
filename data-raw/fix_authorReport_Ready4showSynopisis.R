authorReport_Ready4showSynopsis <- function (x,
                                             args_ls = NULL,
                                             consent_1L_chr = "",
                                             consent_indcs_int = 1L,
                                             fl_nm_1L_chr = NA_character_,
                                             options_chr = c("Y", "N"),
                                             rmd_fl_nms_ls = NULL,
                                             what_1L_chr = NA_character_,
                                             ...) {
  if(!is.na(what_1L_chr)){
    x@a_Ready4showPaths@ms_mkdn_dir_1L_chr <- what_1L_chr
    x@a_Ready4showPaths@ms_dir_1L_chr <- what_1L_chr
  }
  if(is.na(fl_nm_1L_chr)){
    if(is.na(x@fl_nm_1L_chr)){
      fl_nm_1L_chr <- "Manuscript"
    }else{
      fl_nm_1L_chr <- x@fl_nm_1L_chr
    }
  }
  x@fl_nm_1L_chr <- fl_nm_1L_chr
  if(!is.null(rmd_fl_nms_ls)){
    x@rmd_fl_nms_ls <- rmd_fl_nms_ls
  }
  paths_ls <- manufacture(x,
                          what_1L_chr = "paths_ls")
  write_new_dirs(paths_ls %>%
                   purrr::flatten_chr(),
                 consent_1L_chr = consent_1L_chr,
                 consent_indcs_int = consent_indcs_int,
                 options_chr = options_chr)
  header_yaml_args_ls <- make_header_yaml_args_ls(x@authors_r3,
                                                  institutes_tb = x@institutes_r3,
                                                  keywords_chr = x@keywords_chr,
                                                  title_1L_chr = x@title_1L_chr) # Fake data / fl_nm ?
  if (identical(x@abstract_args_ls,list())) {
    x@abstract_args_ls <- NULL
  }
  write_header_fls(path_to_header_dir_1L_chr = paste0(paths_ls$path_to_ms_mkdn_dir_1L_chr,
                                                      "/Header"),
                   header_yaml_args_ls = header_yaml_args_ls,
                   abstract_args_ls = x@abstract_args_ls,
                   consent_1L_chr = consent_1L_chr,
                   consent_indcs_int = consent_indcs_int,
                   options_chr = options_chr)
  write_custom_authors(paths_ls,
                       rmd_fl_nms_ls = x@rmd_fl_nms_ls,
                       consent_1L_chr = consent_1L_chr,
                       consent_indcs_int = consent_indcs_int,
                       options_chr = options_chr)
  params_ls <- append(list(X=x), args_ls) #list(X = x,...)
  output_fl_1L_chr <- paste0(x@fl_nm_1L_chr,
                             ifelse(x@outp_formats_chr[1] == "Word",
                                    ".docx",
                                    paste0(".",
                                           tolower(x@outp_formats_chr[1]))))
  file_to_render_1L_chr <- paste0(paths_ls$path_to_ms_mkdn_dir_1L_chr,
                                  "/Parent_",
                                  x@outp_formats_chr[1],
                                  "/",
                                  x@rmd_fl_nms_ls[names(x@rmd_fl_nms_ls)==x@outp_formats_chr[1]][[1]],
                                  ".Rmd")
  ready4::write_with_consent(consented_fn = rmarkdown::render,
                             prompt_1L_chr = paste0("Do you confirm that you want to render the file ",
                                                    file_to_render_1L_chr,
                                                    "?"),
                             consent_1L_chr = consent_1L_chr,
                             consent_indcs_int = consent_indcs_int,
                             consented_args_ls = list(input = file_to_render_1L_chr,
                                                      output_format = NULL,
                                                      params = params_ls,
                                                      output_file = output_fl_1L_chr,
                                                      output_dir = paths_ls$path_to_ms_outp_dir_1L_chr),
                             consented_msg_1L_chr = paste0("File ",
                                                           file_to_render_1L_chr,
                                                           " has been rendered."),
                             declined_msg_1L_chr = "Render request cancelled.",
                             options_chr = options_chr)

}
