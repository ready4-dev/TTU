add_utility_predn_to_ds <- function (data_tb, mdl, tfmn_1L_chr) 
{
    dep_vars_chr <- c(outp_smry_ls$dep_var_nm_1L_chr, transform_dep_var_nm(dep_var_nm_1L_chr = outp_smry_ls$dep_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)) %>% unique()
    data_tb <- purrr::reduce(dep_vars_chr, .init = data_tb, ~dplyr::mutate(.x, 
        `:=`(!!rlang::sym(.y), NA_real_)))
    data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(outp_smry_ls$dep_var_nm_1L_chr), 
        predict_utility(data_tb = data_tb, tfmn_1L_chr = tfmn_1L_chr, 
            mdl = mdl)))
    return(data_tb)
}
