#' Write brm mdl plt file
#' @description write_brm_mdl_plt_fl() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write brm mdl plt file. The function returns Path to plot (a character vector of length one).
#' @param plt_fn Plt (a function)
#' @param fn_args_ls Function arguments (a list), Default: NULL
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param plt_nm_1L_chr Plt name (a character vector of length one)
#' @param grpx_fn Grpx (a function), Default: grDevices::png
#' @param units_1L_chr Units (a character vector of length one), Default: 'in'
#' @param width_1L_dbl Width (a double vector of length one), Default: 6
#' @param height_1L_dbl Height (a double vector of length one), Default: 6
#' @param rsl_1L_dbl Rsl (a double vector of length one), Default: 300
#' @return Path to plot (a character vector of length one)
#' @rdname write_brm_mdl_plt_fl
#' @export 
#' @importFrom grDevices png
#' @importFrom rlang exec
#' @keywords internal
write_brm_mdl_plt_fl <- function (plt_fn, fn_args_ls = NULL, path_to_write_to_1L_chr, 
    plt_nm_1L_chr, grpx_fn = grDevices::png, units_1L_chr = "in", 
    width_1L_dbl = 6, height_1L_dbl = 6, rsl_1L_dbl = 300) 
{
    path_to_plot_1L_chr <- paste0(path_to_write_to_1L_chr, "/", 
        plt_nm_1L_chr, ifelse(identical(grpx_fn, grDevices::png), 
            ".png", ".tiff"))
    rlang::exec(grpx_fn, !!!list(path_to_plot_1L_chr, units = units_1L_chr, 
        width = width_1L_dbl, height = height_1L_dbl, res = rsl_1L_dbl))
    plt <- rlang::exec(plt_fn, !!!fn_args_ls)
    print(plt)
    dev.off()
    return(path_to_plot_1L_chr)
}
#' Write brm model plts
#' @description write_brm_model_plts() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write brm model plts. The function returns Mdl plts paths (a list).
#' @param mdl_ls Mdl (a list)
#' @param tfd_data_tb Transformed data (a tibble)
#' @param mdl_nm_1L_chr Mdl name (a character vector of length one)
#' @param path_to_write_to_1L_chr Path to write to (a character vector of length one)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param dep_var_desc_1L_chr Dep var description (a character vector of length one), Default: 'AQoL-6D utility score'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param tfmn_fn Tfmn (a function), Default: function(x) {
#'    x
#'}
#' @param units_1L_chr Units (a character vector of length one), Default: 'in'
#' @param height_dbl Height (a double vector), Default: c(rep(6, 2), rep(5, 2))
#' @param width_dbl Width (a double vector), Default: c(rep(6, 2), rep(6, 2))
#' @param rsl_dbl Rsl (a double vector), Default: rep(300, 4)
#' @param args_ls Arguments (a list), Default: NULL
#' @param seed_1L_dbl Seed (a double vector of length one), Default: 23456
#' @return Mdl plts paths (a list)
#' @rdname write_brm_model_plts
#' @export 
#' @importFrom purrr map
#' @importFrom stats setNames
#' @keywords internal
write_brm_model_plts <- function (mdl_ls, tfd_data_tb, mdl_nm_1L_chr, path_to_write_to_1L_chr, 
    dep_var_nm_1L_chr = "aqol6d_total_w", dep_var_desc_1L_chr = "AQoL-6D utility score", 
    round_var_nm_1L_chr = "round", tfmn_fn = function(x) {
        x
    }, units_1L_chr = "in", height_dbl = c(rep(6, 2), rep(5, 
        2)), width_dbl = c(rep(6, 2), rep(6, 2)), rsl_dbl = rep(300, 
        4), args_ls = NULL, seed_1L_dbl = 23456) 
{
    set.seed(seed_1L_dbl)
    tfd_data_tb$Predicted <- predict(mdl_ls)[, 1] %>% tfmn_fn()
    plt_nms_chr <- paste0(mdl_nm_1L_chr, "_", c("coefs", "hetg", 
        "dnst", "sctr_plt"))
    mdl_plts_paths_ls <- purrr::map(1:4, ~{
        if (.x %in% c(1, 2)) {
            plt_fn <- function(mdl_ls, idx_1L_int) {
                plot(mdl_ls, ask = F, plot = F)[idx_1L_int]
            }
            fn_args_ls <- list(mdl_ls = mdl_ls, idx_1L_int = as.integer(.x))
        }
        else {
            if (.x == 3) {
                plt_fn <- plot_obsd_predd_dnst
                fn_args_ls <- list(tfd_data_tb = tfd_data_tb)
            }
            else {
                plt_fn <- plot_obsd_predd_sctr
                fn_args_ls <- list(tfd_data_tb = tfd_data_tb, 
                  dep_var_nm_1L_chr = dep_var_nm_1L_chr, dep_var_desc_1L_chr = dep_var_desc_1L_chr, 
                  args_ls = args_ls)
            }
        }
        write_brm_mdl_plt_fl(plt_fn, fn_args_ls = fn_args_ls, 
            path_to_write_to_1L_chr = path_to_write_to_1L_chr, 
            plt_nm_1L_chr = plt_nms_chr[.x], units_1L_chr = units_1L_chr, 
            width_1L_dbl = width_dbl[.x], height_1L_dbl = height_dbl[.x], 
            rsl_1L_dbl = rsl_dbl[.x])
    }) %>% stats::setNames(plt_nms_chr)
    return(mdl_plts_paths_ls)
}
#' Write results to comma separated variables file
#' @description write_results_to_csv() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write results to comma separated variables file. The function returns Datasets (a tibble).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param output_dir_1L_chr Output directory (a character vector of length one), Default: '.'
#' @return Datasets (a tibble)
#' @rdname write_results_to_csv
#' @export 
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map_dfr map_dfc map map_dbl walk2
#' @importFrom stats setNames
#' @importFrom dplyr mutate select everything
write_results_to_csv <- function (synth_data_spine_ls, output_dir_1L_chr = ".") 
{
    measurements_tb <- tibble::tibble(timepoint_nms_chr = synth_data_spine_ls$timepoint_nms_chr, 
        nbr_obs_dbl = synth_data_spine_ls$nbr_obs_dbl)
    var_summ_res_tb <- suppressMessages(purrr::map_dfr(1:length(synth_data_spine_ls$timepoint_nms_chr), 
        ~{
            idx_dbl <- .x
            suppressWarnings({
                synth_data_spine_ls[c(4:6)] %>% purrr::map_dfc(~.x[idx_dbl])
            }) %>% stats::setNames(c("Mean", "SD", "N_Missing")) %>% 
                dplyr::mutate(var_names_chr = synth_data_spine_ls$var_names_chr, 
                  timepoint_nms_chr = synth_data_spine_ls$timepoint_nms_chr[idx_dbl]) %>% 
                dplyr::select(timepoint_nms_chr, var_names_chr, 
                  dplyr::everything())
        }))
    cor_tb_ls <- synth_data_spine_ls$cor_mat_ls %>% purrr::map(~tibble::as_tibble(.x) %>% 
        stats::setNames(synth_data_spine_ls$var_names_chr) %>% 
        dplyr::mutate(var_names_chr = synth_data_spine_ls$var_names_chr) %>% 
        dplyr::select(var_names_chr, dplyr::everything())) %>% 
        stats::setNames(paste0(synth_data_spine_ls$timepoint_nms_chr, 
            "_correlations_tb"))
    var_class_pars_tb <- synth_data_spine_ls[7:9] %>% tibble::as_tibble() %>% 
        dplyr::mutate(min_dbl = purrr::map_dbl(min_max_ls, ~.x[1]), 
            max_dbl = purrr::map_dbl(min_max_ls, ~.x[2])) %>% 
        dplyr::select(var_names_chr, dplyr::everything(), -min_max_ls)
    output_ls <- list(measurements_tb = measurements_tb, var_summ_res_tb = var_summ_res_tb, 
        var_class_pars_tb = var_class_pars_tb) %>% append(cor_tb_ls)
    dss_tb <- tibble::tibble(ds_obj_nm_chr = names(output_ls), 
        title_chr = c("Brief summary table of the number of observations for which data was collected at each study timepoint.", 
            "Summary statistics (Mean, SD and Number Missing) for AQoL6D health utility and six mental health outcome measures for each study timepoint.", 
            "Brief information about the data structure (whether discrete and allowable range) of AQoL6D health utility and six mental health outcome variables.", 
            paste0("Correlation matrix for AQoL6D health utility and six mental health outcome measures at the ", 
                synth_data_spine_ls$timepoint_nms_chr, " study timepoint.")))
    purrr::walk2(output_ls, names(output_ls), ~write.csv(.x, 
        file = paste0(output_dir_1L_chr, "/", .y, ".csv"), row.names = F))
    return(dss_tb)
}
