#' Model types lookup table
#' 
#' A lookup table of abbreviations to describe the different model types supported by TTU functions
#' 
#' A tibble
#' 
#' \describe{
#'   \item{short_name_chr}{Short name (a character vector)}
#'   \item{long_name_chr}{Long name (a character vector)}
#'   \item{control_chr}{Control (a character vector)}
#'   \item{family_chr}{Family (a character vector)}
#'   \item{fn_chr}{Function (a character vector)}
#'   \item{start_chr}{Start (a character vector)}
#'   \item{predn_type_chr}{Prediction type (a character vector)}
#'   \item{tfmn_chr}{Transformation (a character vector)}
#'   \item{tfmn_for_bnml_lgl}{Transformation for binomial (a logical vector)}
#'   \item{fixed_acronym_chr}{Fixed acronym (a character vector)}
#'   \item{mixed_acronym_chr}{Mixed acronym (a character vector)}
#'   \item{mixed_type_chr}{Mixed type (a character vector)}
#'   \item{with_chr}{With (a character vector)}
#' }
"mdl_types_lup"
