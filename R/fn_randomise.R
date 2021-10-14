#' Randomise changes in factor vector levels
#' @description randomise_changes_in_fct_lvls() is a Randomise function that randomly samples from data. Specifically, this function implements an algorithm to randomise changes in factor vector levels. The function is called for its side effects and does not return a value.
#' @param vector_fct Vector (a factor vector)
#' @param prob_unchanged_dbl Probability unchanged (a double vector)
#' @return NULL
#' @rdname randomise_changes_in_fct_lvls
#' @export 
#' @importFrom stats runif
#' @importFrom purrr map2_dbl
#' @keywords internal
randomise_changes_in_fct_lvls <- function (vector_fct, prob_unchanged_dbl) 
{
    labels_chr <- levels(vector_fct)
    levels_dbl <- 1:length(labels_chr)
    unchanged_lgl <- stats::runif(length(vector_fct)) > prob_unchanged_dbl
    vector_fct %>% purrr::map2_dbl(unchanged_lgl, ~{
        ifelse(.y, .x, sample(levels_dbl[levels_dbl != .x], 1))
    }) %>% factor(levels = levels_dbl, labels = labels_chr)
}
