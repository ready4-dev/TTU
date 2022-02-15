#' Get background text
#' @description get_background_text() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get background text. Function argument results_ls specifies the where to look for the required object. The function returns Text (a character vector of length one).
#' @param results_ls Results (a list)
#' @return Text (a character vector of length one)
#' @rdname get_background_text
#' @export 
#' @keywords internal
get_background_text <- function (results_ls) 
{
    text_1L_chr <- results_ls$study_descs_ls$background_1L_chr
    return(text_1L_chr)
}
#' Get candidates for mixed models
#' @description get_cndts_for_mxd_mdls() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get candidates for mixed models. Function argument mdl_types_lup specifies the where to look for the required object. The function returns Candidates for mixed models (a lookup table).
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @return Candidates for mixed models (a lookup table)
#' @rdname get_cndts_for_mxd_mdls
#' @export 
#' @importFrom utils data
#' @importFrom dplyr filter
#' @keywords internal
get_cndts_for_mxd_mdls <- function (mdl_types_lup = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", package = "TTU", envir = environment())
    cndts_for_mxd_mdls_lup <- mdl_types_lup %>% dplyr::filter(!tfmn_for_bnml_lgl, 
        short_name_chr != "BET_LOG")
    return(cndts_for_mxd_mdls_lup)
}
#' Get conclusion text
#' @description get_conclusion_text() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get conclusion text. Function argument results_ls specifies the where to look for the required object. The function returns Text (a character vector of length one).
#' @param results_ls Results (a list)
#' @return Text (a character vector of length one)
#' @rdname get_conclusion_text
#' @export 
#' @keywords internal
get_conclusion_text <- function (results_ls) 
{
    text_1L_chr <- results_ls$study_descs_ls$conclusion_1L_chr
    return(text_1L_chr)
}
#' Get covariate category categoriess
#' @description get_covar_ctgs() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get covariate category categoriess. Function argument results_ls specifies the where to look for the required object. The function returns Covariate category categoriess (a character vector).
#' @param results_ls Results (a list)
#' @param collapse_1L_lgl Collapse (a logical vector of length one), Default: T
#' @return Covariate category categoriess (a character vector)
#' @rdname get_covar_ctgs
#' @export 
#' @importFrom stringi stri_replace_last_fixed
#' @keywords internal
get_covar_ctgs <- function (results_ls, collapse_1L_lgl = T) 
{
    covar_ctgs_chr <- names(results_ls$candidate_covars_ls) %>% 
        tolower()
    if (collapse_1L_lgl) {
        covar_ctgs_chr <- covar_ctgs_chr %>% paste0(collapse = ", ") %>% 
            stringi::stri_replace_last_fixed(",", " and")
    }
    return(covar_ctgs_chr)
}
#' Get covariates by category categories
#' @description get_covars_by_ctg() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get covariates by category categories. Function argument results_ls specifies the where to look for the required object. The function returns Covariates by category categories (a list).
#' @param results_ls Results (a list)
#' @param collapse_1L_lgl Collapse (a logical vector of length one), Default: F
#' @return Covariates by category categories (a list)
#' @rdname get_covars_by_ctg
#' @export 
#' @importFrom purrr map map2
#' @importFrom stats setNames
#' @importFrom stringi stri_replace_last_fixed
#' @importFrom Hmisc capitalize
#' @keywords internal
get_covars_by_ctg <- function (results_ls, collapse_1L_lgl = F) 
{
    covars_by_ctg_ls <- results_ls$candidate_covars_ls %>% purrr::map(~.x) %>% 
        stats::setNames(get_covar_ctgs(results_ls, collapse_1L_lgl = F))
    if (collapse_1L_lgl) {
        covars_by_ctg_ls <- covars_by_ctg_ls %>% purrr::map2(names(covars_by_ctg_ls), 
            ~{
                covars_1L_chr <- .x %>% paste0(collapse = ", ") %>% 
                  stringi::stri_replace_last_fixed(",", " and")
                paste0(ifelse(length(.x) > 1, .y %>% Hmisc::capitalize(), 
                  paste0("The ", .y)), " covariate", ifelse(length(.x) > 
                  1, "s were ", " was "), covars_1L_chr, ".")
            })
    }
    return(covars_by_ctg_ls)
}
#' Get health utility name
#' @description get_hlth_utl_nm() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get health utility name. Function argument results_ls specifies the where to look for the required object. The function returns Health utility name (a character vector of length one).
#' @param results_ls Results (a list)
#' @param short_nm_1L_lgl Short name (a logical vector of length one), Default: T
#' @return Health utility name (a character vector of length one)
#' @rdname get_hlth_utl_nm
#' @export 
#' @keywords internal
get_hlth_utl_nm <- function (results_ls, short_nm_1L_lgl = T) 
{
    health_utl_nm_1L_chr <- ifelse(short_nm_1L_lgl, results_ls$study_descs_ls$health_utl_nm_1L_chr, 
        results_ls$study_descs_ls$health_utl_long_nm_1L_chr)
    return(health_utl_nm_1L_chr)
}
#' Get health utility statistic
#' @description get_hlth_utl_stat() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get health utility statistic. Function argument results_ls specifies the where to look for the required object. The function returns Health utility statistic (a character vector of length one).
#' @param results_ls Results (a list)
#' @param stat_1L_chr Statistic (a character vector of length one), Default: 'bl_mean'
#' @return Health utility statistic (a character vector of length one)
#' @rdname get_hlth_utl_stat
#' @export 
#' @keywords internal
get_hlth_utl_stat <- function (results_ls, stat_1L_chr = "bl_mean") 
{
    hlth_utl_stat_1L_chr <- switch(stat_1L_chr, bl_mean = results_ls$hlth_utl_and_predrs_ls$bl_hu_mean_1L_dbl, 
        bl_sd = results_ls$hlth_utl_and_predrs_ls$bl_hu_sd_1L_dbl, 
        fup_mean = results_ls$hlth_utl_and_predrs_ls$fup_hu_mean_1L_dbl, 
        fup_sd = results_ls$hlth_utl_and_predrs_ls$fup_hu_sd_1L_dbl)
    return(hlth_utl_stat_1L_chr)
}
#' Get link from transformation
#' @description get_link_from_tfmn() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get link from transformation. Function argument tfmn_1L_chr specifies the where to look for the required object. The function returns Link (a character vector of length one).
#' @param tfmn_1L_chr Transformation (a character vector of length one)
#' @param is_OLS_1L_lgl Is ordinary least squares (a logical vector of length one), Default: F
#' @return Link (a character vector of length one)
#' @rdname get_link_from_tfmn
#' @export 
#' @keywords internal
get_link_from_tfmn <- function (tfmn_1L_chr, is_OLS_1L_lgl = F) 
{
    link_1L_chr <- ifelse(is_OLS_1L_lgl, "identity", ifelse(tfmn_1L_chr == 
        "LOG", "log", ifelse(tfmn_1L_chr == "LGT", "logit", ifelse(tfmn_1L_chr == 
        "CLL", "cloglog", ifelse(tfmn_1L_chr == "LOGLOG", "loglog", 
        ifelse(tfmn_1L_chr == "NTF", "identity", "ERROR"))))))
    if (link_1L_chr == "ERROR") 
        stop("Link cannot be identified - incorrect transformation argument tfmn_1L_chr")
    return(link_1L_chr)
}
#' Get longitudinal transfer to utility algorithm types
#' @description get_lngl_ttu_types() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get longitudinal transfer to utility algorithm types. Function argument results_ls specifies the where to look for the required object. The function returns Model types (a character vector).
#' @param results_ls Results (a list)
#' @param collapse_1L_lgl Collapse (a logical vector of length one), Default: T
#' @return Model types (a character vector)
#' @rdname get_lngl_ttu_types
#' @export 
#' @importFrom stringi stri_replace_last
#' @keywords internal
get_lngl_ttu_types <- function (results_ls, collapse_1L_lgl = T) 
{
    mdl_types_chr <- results_ls$ttu_lngl_ls$best_mdls_tb$model_type
    if (collapse_1L_lgl) {
        mdl_types_chr <- mdl_types_chr %>% paste0(collapse = ", ") %>% 
            stringi::stri_replace_last(fixed = ",", " and")
    }
    return(mdl_types_chr)
}
#' Get model comparisons
#' @description get_mdl_cmprsns() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get model comparisons. Function argument results_ls specifies the where to look for the required object. The function returns Model comparisons (an output object of multiple potential types).
#' @param results_ls Results (a list)
#' @param describe_1L_lgl Describe (a logical vector of length one), Default: T
#' @param mixed_1L_lgl Mixed (a logical vector of length one), Default: F
#' @param as_list_1L_lgl As list (a logical vector of length one), Default: F
#' @return Model comparisons (an output object of multiple potential types)
#' @rdname get_mdl_cmprsns
#' @export 
#' @importFrom purrr map map_chr
#' @importFrom stats setNames
#' @importFrom stringi stri_replace_last_fixed
#' @keywords internal
get_mdl_cmprsns <- function (results_ls, describe_1L_lgl = T, mixed_1L_lgl = F, 
    as_list_1L_lgl = F) 
{
    if (as_list_1L_lgl) {
        mdl_types_chr <- c("OLS", "GLM")[c("OLS", "GLM") %in% 
            results_ls$tables_ls$tenf_sngl_predr_tb$Model]
        mdl_cmprsns_ls <- mdl_types_chr %>% purrr::map(~{
            if (.x == "OLS") {
                mdls_chr <- results_ls$tables_ls$tenf_sngl_predr_tb$Model[(which(results_ls$tables_ls$tenf_sngl_predr_tb$Model == 
                  "OLS") + 1):(which(results_ls$tables_ls$tenf_sngl_predr_tb$Model == 
                  "GLM") - 1)] %>% unique()
            }
            if (.x == "GLM") {
                mdls_chr <- results_ls$tables_ls$tenf_sngl_predr_tb$Model[(which(results_ls$tables_ls$tenf_sngl_predr_tb$Model == 
                  "GLM") + 1):length(results_ls$tables_ls$tenf_sngl_predr_tb$Model)] %>% 
                  unique()
            }
            mdls_chr
        }) %>% stats::setNames(mdl_types_chr)
        mdl_cmprsns_xx <- mdl_cmprsns_ls
    }
    else {
        mdl_cmprsns_1L_chr <- paste0(ifelse(!"OLS" %in% results_ls$tables_ls$tenf_sngl_predr_tb$Model, 
            "", ifelse(describe_1L_lgl, paste0("OLS regression models used ", 
                results_ls$tables_ls$tenf_sngl_predr_tb$Model[(which(results_ls$tables_ls$tenf_sngl_predr_tb$Model == 
                  "OLS") + 1):(which(results_ls$tables_ls$tenf_sngl_predr_tb$Model == 
                  "GLM") - 1)] %>% unique() %>% purrr::map_chr(~.x %>% 
                  stringi::stri_replace_last_fixed("(", "(measured on a scale of ")) %>% 
                  paste0(collapse = ", ") %>% stringi::stri_replace_last_fixed(",", 
                  " and") %>% tolower(), "."), ifelse(mixed_1L_lgl, 
                "linear mixed effect models (LMMs)", "ordinary least squares (OLS) regression models"))), 
            ifelse(!describe_1L_lgl & length(intersect(c("OLS", 
                "GLM"), results_ls$tables_ls$tenf_sngl_predr_tb$Model)) == 
                2, " and ", ifelse(describe_1L_lgl, " ", "")), 
            ifelse(!"GLM" %in% results_ls$tables_ls$tenf_sngl_predr_tb$Model, 
                "", ifelse(describe_1L_lgl, paste0("GLMs used ", 
                  results_ls$tables_ls$tenf_sngl_predr_tb$Model[(which(results_ls$tables_ls$tenf_sngl_predr_tb$Model == 
                    "GLM") + 1):length(results_ls$tables_ls$tenf_sngl_predr_tb$Model)] %>% 
                    unique() %>% purrr::map_chr(~.x %>% stringi::stri_replace_last_fixed("(", 
                    "(measured on a scale of ")) %>% paste0(collapse = ", ") %>% 
                    stringi::stri_replace_last_fixed(",", " and") %>% 
                    tolower()), ifelse(mixed_1L_lgl, "generalised linear mixed effect models (GLMMs)", 
                  "generalised linear models (GLMs)"))))
        mdl_cmprsns_xx <- mdl_cmprsns_1L_chr
    }
    return(mdl_cmprsns_xx)
}
#' Get model type from name
#' @description get_mdl_type_from_nm() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get model type from name. Function argument mdl_nm_1L_chr specifies the where to look for the required object. The function returns Model type (a character vector of length one).
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @return Model type (a character vector of length one)
#' @rdname get_mdl_type_from_nm
#' @export 
#' @importFrom utils data
#' @importFrom dplyr pull
#' @importFrom purrr map_lgl
#' @keywords internal
get_mdl_type_from_nm <- function (mdl_nm_1L_chr, mdl_types_lup = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", package = "TTU", envir = environment())
    mdl_type_1L_chr <- (mdl_types_lup %>% dplyr::pull(short_name_chr))[mdl_types_lup %>% 
        dplyr::pull(short_name_chr) %>% purrr::map_lgl(~endsWith(mdl_nm_1L_chr, 
        .x))]
    return(mdl_type_1L_chr)
}
#' Get models with significant covariates
#' @description get_mdls_with_signft_covars() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get models with significant covariates. Function argument outp_smry_ls specifies the where to look for the required object. The function returns Models with significant covariates (a list).
#' @param outp_smry_ls Output summary (a list)
#' @param params_ls_ls Parameters (a list of lists)
#' @return Models with significant covariates (a list)
#' @rdname get_mdls_with_signft_covars
#' @export 
#' @importFrom purrr map flatten map_lgl
#' @importFrom dplyr filter pull
#' @importFrom stats setNames
#' @keywords internal
get_mdls_with_signft_covars <- function (outp_smry_ls, params_ls_ls) 
{
    signft_covars_chr <- outp_smry_ls$mdls_with_covars_smry_tb %>% 
        get_signft_covars(covar_var_nms_chr = params_ls_ls$params_ls$candidate_covar_nms_chr)
    signft_vars_ls <- outp_smry_ls[["mdls_with_covars_smry_tb"]]$Significant %>% 
        purrr::map(~strsplit(.x, " ")) %>% purrr::flatten()
    mdls_with_signft_covars_ls <- signft_covars_chr %>% purrr::map(~{
        covar_nm_1L_chr <- .x
        mdls_chr <- outp_smry_ls$mdls_with_covars_smry_tb %>% 
            dplyr::filter(purrr::map_lgl(signft_vars_ls, ~any(.x == 
                covar_nm_1L_chr))) %>% dplyr::pull(variable)
        mdls_chr
    }) %>% stats::setNames(signft_covars_chr)
    return(mdls_with_signft_covars_ls)
}
#' Get number of predictors
#' @description get_nbr_of_predrs() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get number of predictors. Function argument results_ls specifies the where to look for the required object. The function returns Number of predictors (an output object of multiple potential types).
#' @param results_ls Results (a list)
#' @param as_words_1L_lgl As words (a logical vector of length one), Default: T
#' @return Number of predictors (an output object of multiple potential types)
#' @rdname get_nbr_of_predrs
#' @export 
#' @importFrom purrr map_int
#' @importFrom xfun numbers_to_words
#' @keywords internal
get_nbr_of_predrs <- function (results_ls, as_words_1L_lgl = T) 
{
    nbr_of_predrs_xx <- results_ls$study_descs_ls$predr_ctgs_ls %>% 
        purrr::map_int(~length(.x[.x %in% results_ls$candidate_predrs_chr])) %>% 
        sum()
    if (as_words_1L_lgl) 
        nbr_of_predrs_xx <- nbr_of_predrs_xx %>% xfun::numbers_to_words()
    return(nbr_of_predrs_xx)
}
#' Get number of predictors by category categories
#' @description get_nbr_of_predrs_by_ctg() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get number of predictors by category categories. Function argument results_ls specifies the where to look for the required object. The function returns Predictors by category categories (a character vector of length one).
#' @param results_ls Results (a list)
#' @return Predictors by category categories (a character vector of length one)
#' @rdname get_nbr_of_predrs_by_ctg
#' @export 
#' @importFrom purrr map_lgl map2_chr
#' @importFrom xfun numbers_to_words
#' @importFrom stringi stri_replace_last_fixed
#' @keywords internal
get_nbr_of_predrs_by_ctg <- function (results_ls) 
{
    multiple_1L_lgl <- length(get_predr_ctgs(results_ls, collapse_1L_lgl = F) > 
        1)
    predrs_by_ctg_1L_chr <- results_ls$study_descs_ls$predr_ctgs_ls[names(results_ls$study_descs_ls$predr_ctgs_ls) %>% 
        tolower() %>% purrr::map_lgl(~.x %in% get_predr_ctgs(results_ls, 
        collapse_1L_lgl = F))] %>% purrr::map2_chr(get_predr_ctgs(results_ls, 
        collapse_1L_lgl = F), ~paste0(.y, ifelse(multiple_1L_lgl, 
        paste0(" (", length(.x[.x %in% results_ls$candidate_predrs_chr]) %>% 
            xfun::numbers_to_words(), " measure", ifelse(length(.x[.x %in% 
            results_ls$candidate_predrs_chr]) > 1, "s", ""), 
            ")"), ""))) %>% paste0(collapse = ", ") %>% stringi::stri_replace_last_fixed(",", 
        " and") %>% tolower()
    return(predrs_by_ctg_1L_chr)
}
#' Get number of secondary analyses
#' @description get_nbr_of_scndry_analyses() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get number of secondary analyses. Function argument results_ls specifies the where to look for the required object. The function returns Number of secondary analyses length one (an output object of multiple potential types).
#' @param results_ls Results (a list)
#' @param as_words_1L_lgl As words (a logical vector of length one), Default: T
#' @param capitalise_1L_lgl Capitalise (a logical vector of length one), Default: T
#' @return Number of secondary analyses length one (an output object of multiple potential types)
#' @rdname get_nbr_of_scndry_analyses
#' @export 
#' @importFrom xfun numbers_to_words
#' @importFrom Hmisc capitalize
#' @keywords internal
get_nbr_of_scndry_analyses <- function (results_ls, as_words_1L_lgl = T, capitalise_1L_lgl = T) 
{
    nbr_of_scndry_analyses_1L_xx <- names(results_ls$mdl_ingredients_ls) %>% 
        startsWith("secondary") %>% sum()
    if (as_words_1L_lgl) {
        nbr_of_scndry_analyses_1L_xx <- nbr_of_scndry_analyses_1L_xx %>% 
            xfun::numbers_to_words()
        if (capitalise_1L_lgl) {
            nbr_of_scndry_analyses_1L_xx <- nbr_of_scndry_analyses_1L_xx %>% 
                Hmisc::capitalize()
        }
    }
    return(nbr_of_scndry_analyses_1L_xx)
}
#' Get ordered single cross-sectional models
#' @description get_ordered_sngl_csnl_mdls() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get ordered single cross-sectional models. Function argument results_ls specifies the where to look for the required object. The function returns Ordered single cross-sectional models (a character vector).
#' @param results_ls Results (a list)
#' @param select_int Select (an integer vector), Default: NULL
#' @param collapse_1L_lgl Collapse (a logical vector of length one), Default: F
#' @return Ordered single cross-sectional models (a character vector)
#' @rdname get_ordered_sngl_csnl_mdls
#' @export 
#' @importFrom stringi stri_replace_last
#' @keywords internal
get_ordered_sngl_csnl_mdls <- function (results_ls, select_int = NULL, collapse_1L_lgl = F) 
{
    ordered_sngl_csnl_mdls_chr <- results_ls$ttu_cs_ls$cs_mdls_predrs_seq_dscdng_chr
    if (!is.null(select_int)) {
        ordered_sngl_csnl_mdls_chr <- ordered_sngl_csnl_mdls_chr[select_int]
    }
    if (collapse_1L_lgl) {
        ordered_sngl_csnl_mdls_chr <- ordered_sngl_csnl_mdls_chr %>% 
            paste0(collapse = ", ") %>% stringi::stri_replace_last(fixed = ",", 
            " and")
    }
    return(ordered_sngl_csnl_mdls_chr)
}
#' Get population descriptives
#' @description get_popl_descvs() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get population descriptives. Function argument results_ls specifies the where to look for the required object. The function returns Population descriptives (a character vector of length one).
#' @param results_ls Results (a list)
#' @return Population descriptives (a character vector of length one)
#' @rdname get_popl_descvs
#' @export 
#' @importFrom stringi stri_replace_last_fixed
#' @keywords internal
get_popl_descvs <- function (results_ls) 
{
    popl_descvs_1L_chr <- results_ls$tables_ls$participant_descs$variable %>% 
        unique() %>% paste0(collapse = ", ") %>% stringi::stri_replace_last_fixed(",", 
        " and")
    return(popl_descvs_1L_chr)
}
#' Get predictor category categoriess
#' @description get_predr_ctgs() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get predictor category categoriess. Function argument results_ls specifies the where to look for the required object. The function returns Predictor category categoriess (a character vector).
#' @param results_ls Results (a list)
#' @param collapse_1L_lgl Collapse (a logical vector of length one), Default: T
#' @return Predictor category categoriess (a character vector)
#' @rdname get_predr_ctgs
#' @export 
#' @importFrom purrr map_int
#' @importFrom stringi stri_replace_last_fixed
#' @keywords internal
get_predr_ctgs <- function (results_ls, collapse_1L_lgl = T) 
{
    predr_ctgs_chr <- (results_ls$study_descs_ls$predr_ctgs_ls %>% 
        names())[(results_ls$study_descs_ls$predr_ctgs_ls %>% 
        purrr::map_int(~length(.x[.x %in% results_ls$candidate_predrs_chr]))) > 
        0] %>% tolower()
    if (collapse_1L_lgl) {
        predr_ctgs_chr <- predr_ctgs_chr %>% paste0(collapse = ", ") %>% 
            stringi::stri_replace_last_fixed(",", " and")
    }
    return(predr_ctgs_chr)
}
#' Get predictors by category categories
#' @description get_predrs_by_ctg() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get predictors by category categories. Function argument results_ls specifies the where to look for the required object. The function returns Predictors by category categories (a list).
#' @param results_ls Results (a list)
#' @param long_desc_1L_lgl Long description (a logical vector of length one), Default: F
#' @param transform_1L_lgl Transform (a logical vector of length one), Default: F
#' @param collapse_1L_lgl Collapse (a logical vector of length one), Default: F
#' @return Predictors by category categories (a list)
#' @rdname get_predrs_by_ctg
#' @export 
#' @importFrom purrr map_lgl map map2 map_chr flatten_chr
#' @importFrom stats setNames
#' @importFrom Hmisc capitalize
#' @importFrom ready4 get_from_lup_obj
#' @importFrom stringi stri_replace_last_fixed stri_replace_last
#' @keywords internal
get_predrs_by_ctg <- function (results_ls, long_desc_1L_lgl = F, transform_1L_lgl = F, 
    collapse_1L_lgl = F) 
{
    predrs_by_ctg_ls <- results_ls$study_descs_ls$predr_ctgs_ls[names(results_ls$study_descs_ls$predr_ctgs_ls) %>% 
        tolower() %>% purrr::map_lgl(~.x %in% get_predr_ctgs(results_ls, 
        collapse_1L_lgl = F))] %>% purrr::map(~.x[.x %in% results_ls$candidate_predrs_chr]) %>% 
        stats::setNames(get_predr_ctgs(results_ls, collapse_1L_lgl = F))
    if (long_desc_1L_lgl) {
        predrs_by_ctg_ls <- predrs_by_ctg_ls %>% purrr::map2(names(predrs_by_ctg_ls) %>% 
            Hmisc::capitalize(), ~{
            predr_descs_1L_chr <- .x %>% purrr::map_chr(~paste0(ready4::get_from_lup_obj(results_ls$mdl_ingredients_ls$dictionary_tb, 
                match_value_xx = .x, match_var_nm_1L_chr = "var_nm_chr", 
                target_var_nm_1L_chr = "var_desc_chr", evaluate_1L_lgl = F), 
                " (", .x %>% transform_names(rename_lup = results_ls$var_nm_change_lup), 
                " - measured on a scale of ", ready4::get_from_lup_obj(results_ls$mdl_ingredients_ls$predictors_lup, 
                  match_value_xx = .x, match_var_nm_1L_chr = "short_name_chr", 
                  target_var_nm_1L_chr = "min_val_dbl", evaluate_1L_lgl = F), 
                "-", ready4::get_from_lup_obj(results_ls$mdl_ingredients_ls$predictors_lup, 
                  match_value_xx = .x, match_var_nm_1L_chr = "short_name_chr", 
                  target_var_nm_1L_chr = "max_val_dbl", evaluate_1L_lgl = F), 
                ")")) %>% paste0(collapse = ", ") %>% stringi::stri_replace_last_fixed(",", 
                " and")
            paste0(.y, " was measured by ", predr_descs_1L_chr, 
                ".")
        })
    }
    else {
        if (transform_1L_lgl) {
            predrs_by_ctg_ls <- predrs_by_ctg_ls %>% purrr::map(~{
                purrr::map_chr(.x, ~transform_names(.x, rename_lup = results_ls$var_nm_change_lup))
            }) %>% purrr::flatten_chr()
            if (collapse_1L_lgl) {
                predrs_by_ctg_ls <- predrs_by_ctg_ls %>% paste0(collapse = ", ") %>% 
                  stringi::stri_replace_last(fixed = ",", " and")
            }
        }
    }
    return(predrs_by_ctg_ls)
}
#' Get preferred model predictors
#' @description get_prefd_mdl_predrs() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get preferred model predictors. Function argument results_ls specifies the where to look for the required object. The function returns Predictors (a character vector of length one).
#' @param results_ls Results (a list)
#' @return Predictors (a character vector of length one)
#' @rdname get_prefd_mdl_predrs
#' @export 
#' @importFrom stringi stri_replace_last
#' @keywords internal
get_prefd_mdl_predrs <- function (results_ls) 
{
    predrs_1L_chr <- results_ls$predr_var_nms_chr %>% paste0(collapse = ", ") %>% 
        stringi::stri_replace_last(fixed = ",", " and")
    return(predrs_1L_chr)
}
#' Get random intercept
#' @description get_random_intercept() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get random intercept. Function argument mdls_smry_tb specifies the where to look for the required object. The function returns Standard deviation (a double vector).
#' @param mdls_smry_tb Models summary (a tibble)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param deterministic_1L_lgl Deterministic (a logical vector of length one), Default: T
#' @return Standard deviation (a double vector)
#' @rdname get_random_intercept
#' @export 
#' @importFrom dplyr filter
#' @importFrom ready4 get_from_lup_obj
#' @keywords internal
get_random_intercept <- function (mdls_smry_tb, mdl_nm_1L_chr, deterministic_1L_lgl = T) 
{
    mdl_smry_tb <- mdls_smry_tb %>% dplyr::filter(Model == mdl_nm_1L_chr)
    sd_dbl <- c(mdl_smry_tb %>% ready4::get_from_lup_obj(match_value_xx = "SD (Intercept)", 
        match_var_nm_1L_chr = "Parameter", target_var_nm_1L_chr = "Estimate", 
        evaluate_1L_lgl = F), ifelse(deterministic_1L_lgl, 0, 
        mdl_smry_tb %>% ready4::get_from_lup_obj(match_value_xx = "SD (Intercept)", 
            match_var_nm_1L_chr = "Parameter", target_var_nm_1L_chr = "SE", 
            evaluate_1L_lgl = F)))
    return(sd_dbl)
}
#' Get secondary analysis descriptions
#' @description get_scndry_anlys_descs() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get secondary analysis descriptions. Function argument results_ls specifies the where to look for the required object. The function returns Secondary analysis descriptions (a character vector).
#' @param results_ls Results (a list)
#' @return Secondary analysis descriptions (a character vector)
#' @rdname get_scndry_anlys_descs
#' @export 
#' @importFrom purrr map_chr pluck
#' @importFrom ready4 get_from_lup_obj
#' @importFrom ready4use remove_labels_from_ds
#' @importFrom stringi stri_replace_last_fixed
#' @keywords internal
get_scndry_anlys_descs <- function (results_ls) 
{
    nbr_of_scndry_analyses_1L_int <- get_nbr_of_scndry_analyses(results_ls, 
        as_words_1L_lgl = F)
    if (nbr_of_scndry_analyses_1L_int > 0) {
        scndry_anlys_descs_chr <- 1:nbr_of_scndry_analyses_1L_int %>% 
            purrr::map_chr(~{
                secondary_ls <- results_ls$mdl_ingredients_ls %>% 
                  purrr::pluck(paste0("secondary_", .x))
                mdls_lup <- secondary_ls$mdls_lup
                predictors_chr <- mdls_lup$predrs_ls %>% unique() %>% 
                  purrr::map_chr(~{
                    .x %>% purrr::map_chr(~ready4::get_from_lup_obj(secondary_ls$dictionary_tb %>% 
                      ready4use::remove_labels_from_ds(), match_value_xx = .x, 
                      match_var_nm_1L_chr = "var_nm_chr", target_var_nm_1L_chr = "var_desc_chr", 
                      evaluate_1L_lgl = F)) %>% paste0(collapse = ", ") %>% 
                      stringi::stri_replace_last_fixed(",", " and")
                  })
                paste0(ifelse(nbr_of_scndry_analyses_1L_int == 
                  1, "The secondary analysis used ", paste0("Secondary Analysis ", 
                  LETTERS[.x], " used ")), ifelse(length(predictors_chr) == 
                  1, paste0(predictors_chr, " as a predictor."), 
                  paste0(predictors_chr, " as predictors.")))
            })
    }
    return(scndry_anlys_descs_chr)
}
#' Get selected mixed models
#' @description get_selected_mixed_mdls() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get selected mixed models. Function argument results_ls specifies the where to look for the required object. The function returns Mixed models (an output object of multiple potential types).
#' @param results_ls Results (a list)
#' @param collapse_1L_lgl Collapse (a logical vector of length one), Default: T
#' @return Mixed models (an output object of multiple potential types)
#' @rdname get_selected_mixed_mdls
#' @export 
#' @importFrom purrr pmap_chr
#' @importFrom stringi stri_replace_last
#' @keywords internal
get_selected_mixed_mdls <- function (results_ls, collapse_1L_lgl = T) 
{
    mixed_mdls_xx <- results_ls$ttu_lngl_ls$best_mdls_tb %>% 
        purrr::pmap_chr(~paste0(..1, " (", ..2, ")"))
    if (collapse_1L_lgl) {
        mixed_mdls_xx <- mixed_mdls_xx %>% paste0(collapse = ", ") %>% 
            stringi::stri_replace_last(fixed = ",", " and")
    }
    return(mixed_mdls_xx)
}
#' Get significant covariates
#' @description get_signft_covars() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get significant covariates. Function argument mdls_with_covars_smry_tb specifies the where to look for the required object. The function returns Signt covariates (a character vector).
#' @param mdls_with_covars_smry_tb Models with covariates summary (a tibble)
#' @param covar_var_nms_chr Covariate variable names (a character vector)
#' @return Signt covariates (a character vector)
#' @rdname get_signft_covars
#' @export 
#' @importFrom purrr map flatten flatten_chr
#' @keywords internal
get_signft_covars <- function (mdls_with_covars_smry_tb, covar_var_nms_chr) 
{
    signif_vars_chr <- mdls_with_covars_smry_tb$Significant %>% 
        purrr::map(~strsplit(.x, " ")) %>% purrr::flatten() %>% 
        purrr::flatten_chr() %>% unique()
    signt_covars_chr <- covar_var_nms_chr[covar_var_nms_chr %in% 
        signif_vars_chr]
    if (identical(signt_covars_chr, character(0))) 
        signt_covars_chr <- NA_character_
    return(signt_covars_chr)
}
#' Get table prediction model
#' @description get_table_predn_mdl() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get table prediction model. Function argument mdl_nm_1L_chr specifies the where to look for the required object. The function returns Table prediction (a model).
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param ingredients_ls Ingredients (a list)
#' @param analysis_1L_chr Analysis (a character vector of length one), Default: NULL
#' @return Table prediction (a model)
#' @rdname get_table_predn_mdl
#' @export 
#' @importFrom ready4 get_from_lup_obj
#' @importFrom stringr str_sub
#' @importFrom purrr pluck
#' @importFrom dplyr filter
#' @keywords internal
get_table_predn_mdl <- function (mdl_nm_1L_chr, ingredients_ls, analysis_1L_chr = NULL) 
{
    mdl_type_1L_chr <- get_mdl_type_from_nm(mdl_nm_1L_chr, mdl_types_lup = ingredients_ls$mdl_types_lup)
    tfmn_1L_chr <- ready4::get_from_lup_obj(ingredients_ls$mdl_types_lup, 
        match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
        target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F)
    if (is.null(analysis_1L_chr)) {
        fake_ds_tb <- ingredients_ls$fake_ds_tb
    }
    else {
        reference_1L_chr <- ifelse(analysis_1L_chr == "Primary Analysis", 
            "Primary", paste0("secondary_", which(LETTERS == 
                stringr::str_sub(analysis_1L_chr, start = -1))))
        fake_ds_tb <- ingredients_ls %>% purrr::pluck(reference_1L_chr) %>% 
            purrr::pluck("fake_ds_tb")
    }
    fake_ds_tb <- fake_ds_tb %>% add_tfd_var_to_ds(depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)
    table_predn_mdl <- make_shareable_mdl(fake_ds_tb = fake_ds_tb, 
        mdl_smry_tb = ingredients_ls$mdls_smry_tb %>% dplyr::filter(Model == 
            mdl_nm_1L_chr), depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr, 
        id_var_nm_1L_chr = ingredients_ls$id_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr, 
        mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = ingredients_ls$mdl_types_lup, 
        control_1L_chr = ready4::get_from_lup_obj(ingredients_ls$mdl_types_lup, 
            match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
            target_var_nm_1L_chr = "control_chr", evaluate_1L_lgl = F), 
        start_1L_chr = NA_character_, seed_1L_int = ingredients_ls$seed_1L_int)
    return(table_predn_mdl)
}
