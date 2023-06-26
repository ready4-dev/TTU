exhibit_TTUProject <- function(x,
                               captions_chr = NULL,
                               display_1L_chr = "all",
                               header_1L_chr = "",
                               header_col_nms_chr = " ",
                               mkdn_tbl_refs_chr = NULL,
                               profile_idx_int = NA_integer_,
                               output_type_1L_chr = "HTML",
                               type_1L_chr = "default",
                               use_lbls_as_col_nms_1L_lgl = T,
                               use_rdocx_1L_lgl = F,
                               variables_chr = character(0),
                               what_1L_chr = "predictors",
                               ...){
  if(what_1L_chr == "predictors"){
    exhibitSlot(x, "b_SpecificParameters@predictors_lup",
                ... = ...)
  }
  if(what_1L_chr == "profile"){
    exhibit(procure(x, variables_chr = variables_chr, type_1L_chr == "default", what_1L_chr = what_1L_chr, ...=...),
            captions_chr = captions_chr,
            header_1L_chr = header_1L_chr,
            header_col_nms_chr = header_col_nms_chr,
            mkdn_tbl_refs_chr = mkdn_tbl_refs_chr,
            profile_idx_int = profile_idx_int,
            output_type_1L_chr = output_type_1L_chr,
            what_1L_chr = type_1L_chr) #descriptives
  }
  if(what_1L_chr == "records"){
    if(type_1L_chr %in% c("ds","dict")){
      exhibit(x@a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad,
              caption_1L_chr = {if(is.null(captions_chr)){NA_character_}else{captions_chr[1]}},
              display_1L_chr = display_1L_chr,
              mkdn_tbl_ref_1L_chr = {if(is.null(mkdn_tbl_refs_chr)){""}else{mkdn_tbl_refs_chr[1]}},
              output_type_1L_chr = "HTML",
              type_1L_chr = type_1L_chr,
              use_lbls_as_col_nms_1L_lgl = T,
              use_rdocx_1L_lgl = F,
              ... = ...)
    }
  }
}
