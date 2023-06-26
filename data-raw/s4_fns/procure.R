procure_TTUProject <- function(x,
                               type_1L_chr = "default",
                               variables_chr = character(0),
                               what_1L_chr = "records",
                               ...){
  if(what_1L_chr == "parameters"){
    if(type_1L_chr == "models_lup"){
      object_xx <- x@b_SpecificParameters@candidate_mdls_lup
    }
  }
  if(what_1L_chr %in% c("project")){
    if(type_1L_chr == "models"){
      object_xx <- procureSlot(x,"c_SpecificProject", use_procure_mthd_1L_lgl = T, what_1L_chr = "prefd_mdls")
    }
  }
  if(what_1L_chr %in% c("profile", "records")){
    object_xx <- x@a_ScorzProfile@a_YouthvarsProfile
    if(!identical(variables_chr, character(0))){
      object_xx@a_Ready4useDyad@ds_tb <- object_xx@a_Ready4useDyad@ds_tb %>% dplyr::select(variables_chr)
    }
    if(what_1L_chr == "records"){
      if(type_1L_chr %in% c("default")){
        object_xx <- object_xx@a_Ready4useDyad@ds_tb
      }
    }
  }
  return(object_xx)
}
