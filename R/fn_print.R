#' Print all plots for model set
#' @description print_all_plts_for_mdl_set() is a Print function that prints output to console Specifically, this function implements an algorithm to print all plots for model set. The function is called for its side effects and does not return a value.
#' @param output_ls Output (a list)
#' @param start_from_1L_int Start from (an integer vector of length one), Default: 0
#' @return NULL
#' @rdname print_all_plts_for_mdl_set
#' @export 
#' @importFrom purrr flatten walk
#' @keywords internal
print_all_plts_for_mdl_set <- function (output_ls, start_from_1L_int = 0L) 
{
    label_refs_mat <- paste0("lab", 1:(output_ls %>% purrr::flatten() %>% 
        length()) + start_from_1L_int) %>% matrix(nrow = length(output_ls), 
        byrow = T)
    purrr::walk(1:length(output_ls), ~print_ts_mdl_plts(output_ls[[.x]], 
        title_1L_chr = names(output_ls)[.x], label_refs_chr = label_refs_mat[.x, 
            ]))
}
#' Print time series model plots
#' @description print_ts_mdl_plts() is a Print function that prints output to console Specifically, this function implements an algorithm to print time series model plots. The function is called for its side effects and does not return a value.
#' @param paths_to_plts_chr Paths to plots (a character vector)
#' @param title_1L_chr Title (a character vector of length one)
#' @param label_refs_chr Label references (a character vector)
#' @param mdl_smry_ls Model summary (a list)
#' @return NULL
#' @rdname print_ts_mdl_plts
#' @export 
#' @importFrom purrr pwalk
#' @keywords internal
print_ts_mdl_plts <- function (paths_to_plts_chr, title_1L_chr, label_refs_chr, mdl_smry_ls) 
{
    cat("\n")
    cat("## ", title_1L_chr, "\n")
    cat("\n")
    purrr::pwalk(list(paths_to_plts_chr, paste0("Caption", 1:length(paths_to_plts_chr)), 
        label_refs_chr), ~{
        cat(paste0("\n![", "](", ..1, ")\n\n"))
        cat(fig_nums(..3, ..2))
        cat("\n\n")
    })
    cat("\n\n")
}
