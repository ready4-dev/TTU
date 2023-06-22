exhibit_TTUProject <- function(x,
                               what_1L_chr = "predictors",
                               ...){
  if(what_1L_chr == "predictors"){
    exhibitSlot(x, "b_SpecificParameters@predictors_lup",
                ... = ...)
  }

}
