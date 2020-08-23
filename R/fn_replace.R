#' Replace var vals with missing table
#' @description replace_var_vals_with_missing_tbl() is a Replace function that edits an object, replacing a specified element with another specified element. Specifically, this function implements an algorithm to replace var vals with missing a table. Function argument tbl specifies the object to be updated. Argument synth_data_spine_ls provides the object to be updated. The function is called for its side effects and does not return a value.
#' @param tbl PARAM_DESCRIPTION
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param idx_int Idx (an integer vector of length 1)
#' @return NULL
#' @rdname replace_var_vals_with_missing_tbl
#' @export 
#' @importFrom purrr reduce
#' @importFrom simstudy defMiss genMiss genObs
#' @keywords internal
replace_var_vals_with_missing_tbl <- function (tbl, synth_data_spine_ls, idx_int) 
{
    missing_def_tbl <- purrr::reduce(1:length(synth_data_spine_ls$var_names_chr), 
        .init = NULL, ~simstudy::defMiss(.x, varname = synth_data_spine_ls$var_names_chr[.y], 
            formula = synth_data_spine_ls$missing_ls[[idx_int]][.y]/synth_data_spine_ls$nbr_obs_dbl[idx_int], 
            logit.link = FALSE))
    missing_mat <- simstudy::genMiss(tbl, missing_def_tbl, idvars = "id")
    synth_tbl <- simstudy::genObs(tbl, missing_mat, idvars = "id")
    return(synth_tbl)
}
