#' Predict utility
#' @description predict_utility() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict utility. The function returns Predd utl (a double vector).
#' @param data_tb Data (a tibble)
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param model_mdl PARAM_DESCRIPTION
#' @param force_min_max_1L_lgl Force min max (a logical vector of length one), Default: T
#' @param utl_min_val_1L_dbl Utl min value (a double vector of length one), Default: 0.03
#' @param impute_1L_lgl Impute (a logical vector of length one), Default: T
#' @param utl_cls_fn Utl class (a function), Default: NULL
#' @return Predd utl (a double vector)
#' @rdname predict_utility
#' @export 
#' @importFrom rlang exec
#' @keywords internal
predict_utility <- function (data_tb, tfmn_1L_chr = "NTF", model_mdl, force_min_max_1L_lgl = T, 
    utl_min_val_1L_dbl = 0.03, impute_1L_lgl = T, utl_cls_fn = NULL) 
{
    predd_utl_dbl <- predict(model_mdl, newdata = data_tb) %>% 
        calculate_dep_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, tfmn_is_outp_1L_lgl = T) %>% 
        as.vector()
    if (impute_1L_lgl) 
        predd_utl_dbl[which(is.na(predd_utl_dbl))] <- predd_utl_dbl %>% 
            na.omit() %>% mean()
    if (force_min_max_1L_lgl) {
        predd_utl_dbl[which(predd_utl_dbl > 1)] <- 1
        predd_utl_dbl[which(predd_utl_dbl < utl_min_val_1L_dbl)] <- utl_min_val_1L_dbl
    }
    if (!is.null(utl_cls_fn)) {
        predd_utl_dbl <- predd_utl_dbl %>% rlang::exec(.fn = utl_cls_fn)
    }
    return(predd_utl_dbl)
}
