#' Randomise ptl fup
#' @description randomise_ptl_fup_fct() is a Randomise function that randomly samples from data. Specifically, this function implements an algorithm to randomise ptl fup. The function is called for its side effects and does not return a value.
#' @param vector_fct Vector (a factor)
#' @param prob_unchanged_dbl Prob unchanged (a double vector of length 1)
#' @return NULL
#' @rdname randomise_ptl_fup_fct
#' @export 
#' @importFrom purrr map2_dbl
#' @keywords internal
randomise_ptl_fup_fct <- function (vector_fct, prob_unchanged_dbl) 
{
    labels_chr <- levels(vector_fct)
    levels_dbl <- 1:length(labels_chr)
    unchanged_lgl <- runif(length(vector_fct)) > 0.5
    vector_fct %>% purrr::map2_dbl(unchanged_lgl, ~{
        ifelse(.y, .x, sample(levels_dbl[levels_dbl != .x], 1))
    }) %>% factor(levels = levels_dbl, labels = labels_chr)
}
