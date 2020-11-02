plot_obsd_predd_dnst <- function(tfd_data_tb){
  tfd_data_tb %>%
    dplyr::mutate(Observed=aqol6d_total_w) %>%
    tidyr::gather(variable, value, Predicted, Observed) %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = variable)) +
    ggalt::geom_bkde(alpha = 0.5) +
    ggplot2::geom_rug() +
    viridis::scale_fill_viridis(discrete = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::labs( x="AQoL-6D utility score",y="Density", fill="")
}
plot_obsd_predd_sctr <- function(tfd_data_tb,
                                 dep_var_nm_1L_chr,
                                 dep_var_desc_1L_chr,
                                 args_ls){
  ggplot2::ggplot(tfd_data_tb) +
    rlang::exec(ggplot2::geom_point,
                ggplot2::aes(x = !!rlang::sym(dep_var_nm_1L_chr),
                             y = Predicted,
                             col = !!rlang::sym(round_var_nm_1L_chr)),
                size = 1,
                !!!args_ls) +
    ggplot2::theme_bw()  +
    ggplot2::xlim(0,1)  +
    ggplot2::ylim(0,1)  +
    ggplot2::scale_color_manual(values=c("#D55E00","#56B4E9"))  +
    ggplot2::labs(x = paste0("Observed ", dep_var_desc_1L_chr),
                  y = paste0("Predicted ", dep_var_desc_1L_chr),
                  col = "") +
    ggplot2::theme(legend.position="bottom")
}
