#' Plot automatic linear model
#' @description plot_auto_lm() is a Plot function that plots data Specifically, this function implements an algorithm to plot automatic linear model. The function is called for its side effects and does not return a value.
#' @param mdl Model (a model)
#' @param which_dbl Which (a double vector), Default: 1:6
#' @param ncol_1L_int Number of columns (an integer vector of length one), Default: 3
#' @param label_size_1L_int Label size (an integer vector of length one), Default: 3
#' @return NULL
#' @rdname plot_auto_lm
#' @export 
#' @importFrom pacman p_load
#' @importFrom ggplot2 autoplot
#' @keywords internal
plot_auto_lm <- function (mdl, which_dbl = 1:6, ncol_1L_int = 3L, label_size_1L_int = 3) 
{
    pacman::p_load(char = "ggfortify")
    plt <- ggplot2::autoplot(mdl, which = which_dbl, ncol = ncol_1L_int, 
        label.size = label_size_1L_int)
    if (6 %in% which_dbl) 
        plt[which(which_dbl == 6)] <- plt[which(which_dbl == 
            6)] + ggtitle("Cook's vs Leverage")
    plt
}
#' Plot linear comparison
#' @description plot_lnr_cmprsn() is a Plot function that plots data Specifically, this function implements an algorithm to plot linear comparison. The function is called for its side effects and does not return a value.
#' @param data_tb Data (a tibble)
#' @param predn_ds_tb Prediction dataset (a tibble)
#' @param predr_var_nm_1L_chr Predictor variable name (a character vector of length one)
#' @param predr_var_desc_1L_chr Predictor variable description (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param depnt_var_desc_1L_chr Dependent variable description (a character vector of length one), Default: 'Total weighted utility score'
#' @return NULL
#' @rdname plot_lnr_cmprsn
#' @export 
#' @importFrom dplyr filter
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_line theme_bw labs
#' @keywords internal
plot_lnr_cmprsn <- function (data_tb, predn_ds_tb, predr_var_nm_1L_chr, predr_var_desc_1L_chr, 
    depnt_var_nm_1L_chr = "utl_total_w", depnt_var_desc_1L_chr = "Total weighted utility score") 
{
    data_tb <- data_tb %>% dplyr::filter(!is.na(!!rlang::sym(predr_var_nm_1L_chr)))
    ggplot2::ggplot(data_tb, ggplot2::aes(x = !!rlang::sym(predr_var_nm_1L_chr), 
        y = !!rlang::sym(depnt_var_nm_1L_chr))) + ggplot2::geom_point() + 
        ggplot2::geom_smooth(method = "loess", size = 1.5) + 
        ggplot2::geom_line(data = predn_ds_tb, ggplot2::aes(x = !!rlang::sym(predr_var_nm_1L_chr), 
            y = !!rlang::sym(depnt_var_nm_1L_chr)), col = "red") + 
        ggplot2::theme_bw() + ggplot2::labs(x = predr_var_desc_1L_chr, 
        y = depnt_var_desc_1L_chr)
}
#' Plot observed predicted density
#' @description plot_obsd_predd_dnst() is a Plot function that plots data Specifically, this function implements an algorithm to plot observed predicted density. The function is called for its side effects and does not return a value.
#' @param tfd_data_tb Transformed data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param depnt_var_desc_1L_chr Dependent variable description (a character vector of length one), Default: 'Total weighted utility score'
#' @param predd_val_var_nm_1L_chr Predicted value variable name (a character vector of length one), Default: 'Predicted'
#' @param cmprsn_predd_var_nm_1L_chr Comparison predicted variable name (a character vector of length one), Default: 'NA'
#' @return NULL
#' @rdname plot_obsd_predd_dnst
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_rug theme_bw theme labs
#' @importFrom ggalt geom_bkde
#' @importFrom viridis scale_fill_viridis
#' @keywords internal
plot_obsd_predd_dnst <- function (tfd_data_tb, depnt_var_nm_1L_chr = "utl_total_w", depnt_var_desc_1L_chr = "Total weighted utility score", 
    predd_val_var_nm_1L_chr = "Predicted", cmprsn_predd_var_nm_1L_chr = NA_character_) 
{
    if (is.na(cmprsn_predd_var_nm_1L_chr)) 
        cmprsn_predd_var_nm_1L_chr <- NULL
    args_ls <- list(predd_val_var_nm_1L_chr, cmprsn_predd_var_nm_1L_chr)
    tfd_data_tb %>% dplyr::mutate(Observed = !!rlang::sym(depnt_var_nm_1L_chr)) %>% 
        tidyr::gather(variable, value, !!!args_ls, Observed) %>% 
        ggplot2::ggplot(ggplot2::aes(x = value, fill = variable)) + 
        ggalt::geom_bkde(alpha = 0.5) + ggplot2::geom_rug() + 
        viridis::scale_fill_viridis(discrete = TRUE) + ggplot2::theme_bw() + 
        ggplot2::theme(legend.position = "bottom") + ggplot2::labs(x = depnt_var_desc_1L_chr, 
        y = "Density", fill = "")
}
#' Plot observed predicted scatter comparison
#' @description plot_obsd_predd_sctr_cmprsn() is a Plot function that plots data Specifically, this function implements an algorithm to plot observed predicted scatter comparison. The function is called for its side effects and does not return a value.
#' @param tfd_data_tb Transformed data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param depnt_var_desc_1L_chr Dependent variable description (a character vector of length one), Default: 'Total weighted utility score'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param args_ls Arguments (a list), Default: NULL
#' @param predd_val_var_nm_1L_chr Predicted value variable name (a character vector of length one), Default: 'Predicted'
#' @return NULL
#' @rdname plot_obsd_predd_sctr_cmprsn
#' @export 
#' @importFrom ggplot2 ggplot geom_point aes theme_bw xlim ylim scale_color_manual labs theme
#' @importFrom rlang exec sym
#' @keywords internal
plot_obsd_predd_sctr_cmprsn <- function (tfd_data_tb, depnt_var_nm_1L_chr = "utl_total_w", depnt_var_desc_1L_chr = "Total weighted utility score", 
    round_var_nm_1L_chr = "round", args_ls = NULL, predd_val_var_nm_1L_chr = "Predicted") 
{
    ggplot2::ggplot(tfd_data_tb) + rlang::exec(ggplot2::geom_point, 
        ggplot2::aes(x = !!rlang::sym(depnt_var_nm_1L_chr), y = !!rlang::sym(predd_val_var_nm_1L_chr), 
            col = !!rlang::sym(round_var_nm_1L_chr)), size = 1, 
        !!!args_ls) + ggplot2::theme_bw() + ggplot2::xlim(0, 
        1) + ggplot2::ylim(0, 1) + ggplot2::scale_color_manual(values = c("#D55E00", 
        "#56B4E9")) + ggplot2::labs(x = paste0("Observed ", depnt_var_desc_1L_chr), 
        y = paste0(predd_val_var_nm_1L_chr, " ", depnt_var_desc_1L_chr), 
        col = "") + ggplot2::theme(legend.position = "bottom")
}
#' Plot scatter plot comparison
#' @description plot_sctr_plt_cmprsn() is a Plot function that plots data Specifically, this function implements an algorithm to plot scatter plot comparison. The function is called for its side effects and does not return a value.
#' @param tfd_data_tb Transformed data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predd_val_var_nm_1L_chr Predicted value variable name (a character vector of length one), Default: 'Predicted'
#' @return NULL
#' @rdname plot_sctr_plt_cmprsn
#' @export 
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_bw geom_abline xlim ylim
#' @importFrom rlang sym
#' @keywords internal
plot_sctr_plt_cmprsn <- function (tfd_data_tb, depnt_var_nm_1L_chr = "utl_total_w", predd_val_var_nm_1L_chr = "Predicted") 
{
    tfd_data_tb %>% ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(depnt_var_nm_1L_chr), 
        y = !!rlang::sym(predd_val_var_nm_1L_chr))) + ggplot2::geom_point() + 
        ggplot2::geom_smooth(method = "loess", size = 1.5) + 
        ggplot2::theme_bw() + ggplot2::geom_abline(intercept = 0, 
        slope = 1, color = "red", linetype = "dashed", size = 1.5) + 
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
}
