manufacture_TTUProject <- function(x,
                                   type_1L_chr = "dummys",
                                   what_1L_chr = "factors",
                                   ...){
  object_xx <- manufacture(x@c_SpecificProject@a_YouthvarsProfile@a_Ready4useDyad, type_1L_chr = type_1L_chr, what_1L_chr = what_1L_chr,
                           restrict_to_chr = x@c_SpecificProject@b_SpecificParameters@candidate_covars_chr)
  return(object_xx)
}
