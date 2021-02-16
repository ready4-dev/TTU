
setOldClass(c("TTU_predictors_lup","tbl_df", "tbl", "data.frame"))
#' TTU S3 class for candidate predictors lookup table
#' @description Create a new valid instance of the TTU S3 class for candidate predictors lookup table
#' @param x A prototype for the TTU S3 class for candidate predictors lookup table, Default: make_pt_TTU_predictors_lup()
#' @return A validated instance of the TTU S3 class for candidate predictors lookup table
#' @details TTU S3 class for candidate predictors lookup table
#' @rdname TTU_predictors_lup
#' @export 

TTU_predictors_lup <- function(x = make_pt_TTU_predictors_lup()){ 
validate_TTU_predictors_lup(make_new_TTU_predictors_lup(x))
}
#' Make new TTU S3 class for candidate predictors lookup table
#' @description Create a new unvalidated instance of the TTU S3 class for candidate predictors lookup table
#' @param x A prototype for the TTU S3 class for candidate predictors lookup table
#' @return An unvalidated instance of the TTU S3 class for candidate predictors lookup table
#' @details TTU S3 class for candidate predictors lookup table
#' @rdname make_new_TTU_predictors_lup
#' @export 
#' @importFrom tibble is_tibble
make_new_TTU_predictors_lup <- function(x){ 
stopifnot(tibble::is_tibble(x))
class(x) <- append(c("TTU_predictors_lup",setdiff(make_pt_TTU_predictors_lup() %>% class(),class(x))),
class(x))
x
}
#' Make prototype TTU S3 class for candidate predictors lookup table
#' @description Create a new prototype for the TTU S3 class for candidate predictors lookup table
#' @param short_name_chr Short name (a character vector), Default: character(0)
#' @param long_name_chr Long name (a character vector), Default: character(0)
#' @param min_val_dbl Min value (a double vector), Default: numeric(0)
#' @param max_val_dbl Max value (a double vector), Default: numeric(0)
#' @param class_chr Class (a character vector), Default: character(0)
#' @param increment_dbl Increment (a double vector), Default: numeric(0)
#' @param class_fn_chr Class function (a character vector), Default: character(0)
#' @param mdl_scaling_dbl Mdl scaling (a double vector), Default: numeric(0)
#' @param covariate_lgl Covariate (a logical vector), Default: logical(0)
#' @return A prototype for TTU S3 class for candidate predictors lookup table
#' @details TTU S3 class for candidate predictors lookup table
#' @rdname make_pt_TTU_predictors_lup
#' @export 
#' @importFrom ready4class update_pt_fn_args_ls
#' @importFrom rlang exec
#' @importFrom tibble tibble
make_pt_TTU_predictors_lup <- function(short_name_chr = character(0),
long_name_chr = character(0),
min_val_dbl = numeric(0),
max_val_dbl = numeric(0),
class_chr = character(0),
increment_dbl = numeric(0),
class_fn_chr = character(0),
mdl_scaling_dbl = numeric(0),
covariate_lgl = logical(0)){ 
args_ls <- list(short_name_chr = short_name_chr,
long_name_chr = long_name_chr,
min_val_dbl = min_val_dbl,
max_val_dbl = max_val_dbl,
class_chr = class_chr,
increment_dbl = increment_dbl,
class_fn_chr = class_fn_chr,
mdl_scaling_dbl = mdl_scaling_dbl,
covariate_lgl = covariate_lgl) %>% ready4class::update_pt_fn_args_ls()
rlang::exec(tibble::tibble,!!!args_ls)
}
#' Validate TTU S3 class for candidate predictors lookup table
#' @description Validate an instance of the TTU S3 class for candidate predictors lookup table
#' @param x An unvalidated instance of the TTU S3 class for candidate predictors lookup table
#' @return A prototpe for TTU S3 class for candidate predictors lookup table
#' @details TTU S3 class for candidate predictors lookup table
#' @rdname validate_TTU_predictors_lup
#' @export 
#' @importFrom stringr str_detect str_c
#' @importFrom dplyr summarise_all arrange filter pull
#' @importFrom tidyr gather
#' @importFrom purrr map2_chr
validate_TTU_predictors_lup <- function(x){
if(sum(stringr::str_detect(names(x)[names(x) %in% names(make_pt_TTU_predictors_lup())],
names(make_pt_TTU_predictors_lup())))!=length(names(make_pt_TTU_predictors_lup()))){
stop(paste0("TIBBLE must include columns named: ",
names(make_pt_TTU_predictors_lup()) %>% stringr::str_c(sep="", collapse = ", ")),
call. = FALSE)
}
 if(!identical(make_pt_TTU_predictors_lup() %>% 
dplyr::summarise_all(class) %>% 
 tidyr::gather(variable,class) %>% 
dplyr::arrange(variable),
x %>% 
dplyr::summarise_all(class) %>% 
 tidyr::gather(variable,class) %>% 
dplyr::filter(variable %in% names(make_pt_TTU_predictors_lup())) %>% dplyr::arrange(variable))){
stop(paste0("TIBBLE columns should be of the following classes: ",
purrr::map2_chr(make_pt_TTU_predictors_lup() %>% 
dplyr::summarise_all(class) %>% 
 tidyr::gather(variable,class) %>% 
dplyr::pull(1),
 make_pt_TTU_predictors_lup() %>% 
dplyr::summarise_all(class) %>% 
 tidyr::gather(variable,class) %>% 
dplyr::pull(2),
 ~ paste0(.x,": ",.y)) %>% 
stringr::str_c(sep="", collapse = ", ")),
call. = FALSE)
}
x}
#' Is TTU S3 class for candidate predictors lookup table
#' @description Check whether an object is a valid instance of the TTU S3 class for candidate predictors lookup table
#' @param x An object of any type
#' @return A logical value, TRUE if a valid instance of the TTU S3 class for candidate predictors lookup table
#' @details TTU S3 class for candidate predictors lookup table
#' @rdname is_TTU_predictors_lup
#' @export 

is_TTU_predictors_lup <- function(x) inherits(validate_TTU_predictors_lup(x), "TTU_predictors_lup")
