author_SpecificSynopsis <- function(x,
                                    args_ls = NULL,
                                    consent_1L_chr = "",
                                    reference_1L_int = NA_integer_,
                                    type_1L_chr = "Report",
                                    what_1L_chr = "Catalogue",
                                    ...){
  if(what_1L_chr == "Catalogue"){
    outp_smry_ls_ls <- manufacture(x@b_SpecificResults,
                                   what_1L_chr = "indexed_shareable")
    refs_int <- 1:length(outp_smry_ls_ls)
    if(!is.na(reference_1L_int)){
      outp_smry_ls_ls <- outp_smry_ls_ls[reference_1L_int]
      refs_int <- reference_1L_int
    }
    ctlg_nms_chr <- purrr::map2_chr(outp_smry_ls_ls,
                                    refs_int,
                                    ~ {
                                      fl_nm_1L_chr <- paste0("AAA_TTU_MDL_CTG",
                                                             ifelse(.y==1,
                                                                    "",
                                                                    paste0("-",(.y-1))))
                                      authorReport_Ready4showSynopsis(x %>% #
                                                                        renewSlot("b_SpecificResults@a_SpecificShareable@shareable_outp_ls",
                                                                                  .x),
                                                                      args_ls = args_ls,
                                                                      consent_1L_chr = consent_1L_chr,
                                                                      fl_nm_1L_chr = fl_nm_1L_chr,
                                                                      what_1L_chr = what_1L_chr#,...
                                                                      )
                                      fl_nm_1L_chr
                                    }
    )
  }
  return(x)
}
