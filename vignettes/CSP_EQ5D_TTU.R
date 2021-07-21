## ----echo = TRUE, message=FALSE-----------------------------------------------
library(magrittr)
library(TTU)

## -----------------------------------------------------------------------------
authors_tb <- ready4show::authors_tb

## ----authorsds, eval = knitr::is_html_output(), echo=F, results='asis'--------
authors_tb %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Authors table",
                          mkdn_tbl_ref_1L_chr = "tab:authorsds")

## -----------------------------------------------------------------------------
institutes_tb <- ready4show::institutes_tb

## ----institutesds, eval = knitr::is_html_output(), echo = F,results='asis'----
institutes_tb %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Author institutions table",
                          mkdn_tbl_ref_1L_chr = "tab:institutesds")

