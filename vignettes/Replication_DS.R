## ----message=F----------------------------------------------------------------
library(magrittr)
library(TTU)

## ----message=FALSE------------------------------------------------------------
dictionary_tb <- repln_ds_dict_r3

## ----ddic---------------------------------------------------------------------
dictionary_tb %>%
  ready4show::print_table(output_type_1L_chr = "HTML", 
                          caption_1L_chr = "Data dictionary",  
                          use_lbls_as_col_nms_1L_lgl = T, 
                          mkdn_tbl_ref_1L_chr = "tab:ddic")

## -----------------------------------------------------------------------------
raw_data_tb <- replication_popl_tb 

## -----------------------------------------------------------------------------
raw_data_tb <- raw_data_tb %>%
  ready4use::add_labels_from_dictionary(dictionary_tb)

## ----repds--------------------------------------------------------------------
raw_data_tb %>%
  head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Replication dataset", 
                          use_lbls_as_col_nms_1L_lgl = T, 
                          mkdn_tbl_ref_1L_chr = "tab:repds") 

