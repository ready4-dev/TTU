#' Transform dep var name for cll
#' @description transform_dep_var_nm_for_cll() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dep var name for cll. Function argument dep_var_nm_1L_chr specifies the object to be updated. The function returns Transformed dep var name (a character vector of length one).
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one)
#' @return Transformed dep var name (a character vector of length one)
#' @rdname transform_dep_var_nm_for_cll
#' @export 

transform_dep_var_nm_for_cll <- function (dep_var_nm_1L_chr) 
{
    tfd_dep_var_nm_1L_chr <- paste0(dep_var_nm_1L_chr, "_cloglog")
    return(tfd_dep_var_nm_1L_chr)
}
#' Transform raw Assessment of Quality of Life tibble to Assessment of Quality of Life Six Dimension
#' @description transform_raw_aqol_tb_to_aqol6d_tb() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform raw assessment of quality of life tibble to assessment of quality of life six dimension tibble. Function argument raw_aqol_tb specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension (a tibble).
#' @param raw_aqol_tb Raw Assessment of Quality of Life (a tibble)
#' @return Assessment of Quality of Life Six Dimension (a tibble)
#' @rdname transform_raw_aqol_tb_to_aqol6d_tb
#' @export 
#' @importFrom dplyr mutate filter select contains rename
#' @importFrom Hmisc label
transform_raw_aqol_tb_to_aqol6d_tb <- function (raw_aqol_tb) 
{
    aqol6d_tb <- raw_aqol_tb %>% dplyr::mutate(d_agegroup = cut(d_age, 
        breaks = c(11, 17, 30), labels = c("Age 12-17", "Age 18-26"))) %>% 
        dplyr::filter(!is.na(aqol6d_total_w)) %>% dplyr::mutate(round = factor(round, 
        labels = c("Baseline", "Follow-up"))) %>% dplyr::select(fkClientID, 
        c_p_diag_s, s_centre, c_clinical_staging_s, d_age, d_agegroup, 
        d_gender, d_sex_birth_s, d_sexual_ori_s, d_country_bir_s, 
        d_ATSI, d_english_home, d_english_native, d_relation_s, 
        d_studying_working, k6_total, phq9_total, bads_total, 
        gad7_total, oasis_total, scared_total, dplyr::contains("aqol6d"), 
        c_sofas, round) %>% dplyr::mutate(Gender = factor(ifelse(d_gender == 
        "Genderqueer/gender nonconforming/agender" | d_gender == 
        "Transgender", "Other", as.character(d_gender)))) %>% 
        dplyr::mutate(Region = as.factor(ifelse(s_centre == "Canberra" | 
            s_centre == "Southport" | s_centre == "Knox", "Metro", 
            "Regional"))) %>% dplyr::mutate(CALD = factor(ifelse(d_country_bir_s == 
        "Other" | d_english_home == "No" | d_english_native == 
        "No" | d_ATSI == "Yes", "Yes", "No"))) %>% dplyr::rename(PHQ9 = phq9_total, 
        BADS = bads_total, GAD7 = gad7_total, OASIS = oasis_total, 
        SCARED = scared_total, K6 = k6_total, SOFAS = c_sofas)
    Hmisc::label(aqol6d_tb$CALD) = "Culturally and linguistically diverse (CALD) background"
    Hmisc::label(aqol6d_tb$d_agegroup) = "Age group"
    return(aqol6d_tb)
}
#' Transform tibble to mdl input
#' @description transform_tb_to_mdl_inp() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform tibble to mdl input. Function argument data_tb specifies the object to be updated. Argument dep_var_nm_1L_chr provides the object to be updated. The function returns Transformed for gsn log mdl (a tibble).
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @return Transformed for gsn log mdl (a tibble)
#' @rdname transform_tb_to_mdl_inp
#' @export 
#' @importFrom dplyr select all_of group_by arrange mutate across first lag
#' @importFrom rlang sym
transform_tb_to_mdl_inp <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", predr_vars_nms_chr, 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline") 
{
    data_tb <- data.frame(data_tb)
    tfd_for_gsn_log_mdl_tb <- data_tb %>% dplyr::select(dplyr::all_of(dep_var_nm_1L_chr), 
        dplyr::all_of(id_var_nm_1L_chr), dplyr::all_of(round_var_nm_1L_chr), 
        dplyr::all_of(predr_vars_nms_chr)) %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>% 
        dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr), !!rlang::sym(round_var_nm_1L_chr))
    tfd_for_gsn_log_mdl_tb <- tfd_for_gsn_log_mdl_tb %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr), 
        .fns = list(baseline = ~dplyr::first(.)/100, change = ~ifelse(round == 
            "Baseline", 0, (. - dplyr::lag(.))/100)))) %>% dplyr::mutate(`:=`(!!rlang::sym(transform_dep_var_nm_for_cll(dep_var_nm_1L_chr)), 
        log(-log(1 - ifelse(!!rlang::sym(dep_var_nm_1L_chr) == 
            1, 0.999, !!rlang::sym(dep_var_nm_1L_chr)))))) %>% 
        na.omit()
    return(tfd_for_gsn_log_mdl_tb)
}
#' Transform ts mdl data
#' @description transform_ts_mdl_data() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform ts mdl data. Function argument mdl_ls specifies the object to be updated. Argument data_tb provides the object to be updated. The function returns Cnfdl mdl (a list).
#' @param mdl_ls Mdl (a list)
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param mdl_nm_1L_chr Mdl name (a character vector of length one)
#' @return Cnfdl mdl (a list)
#' @rdname transform_ts_mdl_data
#' @export 
#' @importFrom dplyr select all_of summarise across everything
#' @importFrom purrr map flatten_chr
transform_ts_mdl_data <- function (mdl_ls, data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", 
    predr_vars_nms_chr, id_var_nm_1L_chr = "fkClientID", mdl_nm_1L_chr) 
{
    old_data_tb <- data_tb %>% dplyr::select(c(dplyr::all_of(id_var_nm_1L_chr), 
        dplyr::all_of(dep_var_nm_1L_chr), predr_vars_nms_chr %>% 
            purrr::map(~paste0(.x, c("", "_baseline", "_change"))) %>% 
            purrr::flatten_chr()))
    cnfdl_mdl_ls <- mdl_ls
    cnfdl_mdl_ls$data <- old_data_tb %>% as.data.frame() %>% 
        dplyr::summarise(dplyr::across(dplyr::everything(), ~sample(.x, 
            1)))
    return(cnfdl_mdl_ls)
}
