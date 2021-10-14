#' Knit from template
#' @description knit_from_tmpl() is a Knit function that knits a rmarkdown file Specifically, this function implements an algorithm to knit from template. The function is called for its side effects and does not return a value.
#' @param params_to_expand_ls Parameters to expand (a list)
#' @param path_to_tmpl_1L_chr Path to template (a character vector of length one)
#' @return NULL
#' @rdname knit_from_tmpl
#' @export 
#' @importFrom purrr map
#' @importFrom rlang exec
#' @importFrom knitr knit_expand knit_child
#' @keywords internal
knit_from_tmpl <- function (params_to_expand_ls, path_to_tmpl_1L_chr) 
{
    src <- purrr::map(params_to_expand_ls, ~{
        args_ls <- append(list(file = path_to_tmpl_1L_chr), .x)
        rlang::exec(knitr::knit_expand, !!!args_ls)
    })
    res <- knitr::knit_child(text = unlist(src), quiet = TRUE)
    cat(res, sep = "\n")
}
#' Knit model report
#' @description knit_mdl_rprt() is a Knit function that knits a rmarkdown file Specifically, this function implements an algorithm to knit model report. The function is called for its side effects and does not return a value.
#' @param knit_pars_ls Knit pars (a list)
#' @param path_to_mdl_rprt_tmpl_1L_chr Path to model report template (a character vector of length one)
#' @return NULL
#' @rdname knit_mdl_rprt
#' @export 
#' @importFrom purrr pmap
#' @importFrom knitr knit_expand knit_child
#' @keywords internal
knit_mdl_rprt <- function (knit_pars_ls, path_to_mdl_rprt_tmpl_1L_chr) 
{
    src <- purrr::pmap(append(knit_pars_ls, list(names(knit_pars_ls$plt_nms_ls))), 
        ~knitr::knit_expand(path_to_mdl_rprt_tmpl_1L_chr, plt_nms_chr = ..1 %>% 
            deparse(), path_to_mdl_1L_chr = ..2 %>% deparse(), 
            caption_1L_chr = ..3 %>% deparse(), label_stub_1L_chr = ..4 %>% 
                deparse(), output_type_1L_chr = ..5 %>% deparse(), 
            section_ttl_1L_chr = ..6, paths_to_mdl_plts_chr = ..7 %>% 
                deparse(), mdl_nms_chr = ..8 %>% deparse()))
    res <- knitr::knit_child(text = unlist(src), quiet = TRUE)
    cat(res, sep = "\n")
}
