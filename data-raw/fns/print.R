print_all_plts_for_mdl_set <- function (output_ls, start_from_1L_int = 0L)
{
    label_refs_mat <- paste0("lab", 1:(output_ls %>% purrr::flatten() %>%
        length()) + start_from_1L_int) %>% matrix(nrow = length(output_ls),
        byrow = T)
    purrr::walk(1:length(output_ls), ~print_ts_mdl_plts(output_ls[[.x]],
        title_1L_chr = names(output_ls)[.x], label_refs_chr = label_refs_mat[.x,
            ]))
}
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
