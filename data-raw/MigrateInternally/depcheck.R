library(TTU)
# library(dplyr)
library(DependenciesGraphs)
# library(mvbutils)
# deps <- funDependencies("package:dplyr", "mutate")
# allFunctionEnv("package:TTU")
# test_tb<-linksForAll("package:TTU")
# graphfun <- mvbutils::foodweb(where = "package:TTU", descendents = F,
#                               plotting = F, ancestors = T, plotting = F)$funmat
pkg_depcy_ls <- DependenciesGraphs::envirDependencies("package:TTU")
# saveRDS(migrate_ls,
#         "migrate_ls.RDS")
# write_fns_to_split_destns <- function(pkg_depcy_ls,
#                                       pkg_1_core_fns_chr,
#                                       original_pkg_nm_1L_chr = get_dev_pkg_nm(),
#                                       pkg_1_nm_1L_chr = "package_1",
#                                       pkg_2_nm_1L_chr = "package_2",
#                                       tmp_dir_path_1L_chr = "data-raw/pkg_migration",
#                                       path_to_fns_dir_1L_chr = "data-raw/fns"){
#   utils::data("fns_dmt_tb", package = original_pkg_nm_1L_chr, envir = environment())
#   read_fns(path_to_fns_dir_1L_chr)
#   fns_for_pkg_1_chr <- get_all_depcys_of_fns(pkg_depcy_ls = pkg_depcy_ls,
#                                                fns_chr = pkg_1_core_fns_chr)
#   fns_for_pkg_2_chr <- setdiff(pkg_depcy_ls$Nomfun$label, fns_for_pkg_1_chr)
#   migrate_ls <- list(fns_for_pkg_1_chr = fns_for_pkg_1_chr,
#                      fns_for_pkg_2_chr = fns_for_pkg_2_chr)
#   if(!dir.exists(tmp_dir_path_1L_chr))
#     dir.create(tmp_dir_path_1L_chr)
#   new_dest_dir_chr <- purrr::map_chr(c(pkg_1_nm_1L_chr,pkg_2_nm_1L_chr),
#                                      ~ {
#                                        new_dir_1L_chr <- paste0(tmp_dir_path_1L_chr,
#                                                                 "/",
#                                                                 .x)
#                                        if(!dir.exists(new_dir_1L_chr))
#                                          dir.create(new_dir_1L_chr)
#                                        new_dir_1L_chr
#                                      })
#   migrate_ls %>%
#     purrr::walk2(new_dest_dir_chr,
#                  ~{
#                    fns_tb <- fns_dmt_tb %>%
#                      dplyr::filter(fns_chr %in% .x) %>%
#                      dplyr::select(fns_chr, file_nm_chr)
#                    file_nms_chr <- fns_tb$file_nm_chr %>% unique()
#                    new_dest_dir_1L_chr <- .y
#                    file_nms_chr %>%
#                      purrr::walk(~
#                                    {
#                                      tb <- fns_tb %>%
#                                        dplyr::filter(file_nm_chr == .x)
#                                      first_lgl_vec <- c(T,rep(F,nrow(tb)-1))
#                                      dest_path_1L_chr <- paste0(new_dest_dir_1L_chr,"/",.x)
#                                      purrr::walk(1:nrow(tb),
#                                                  ~
#                                                    {
#                                                      fn <- eval(parse(text=tb[[.x,1]]))
#                                                      fn_chr <- deparse(fn)
#                                                      sink(dest_path_1L_chr, append =  !first_lgl_vec[.x])
#                                                      writeLines(paste0(tb[[.x,1]]," <- ",fn_chr[1]))
#                                                      writeLines(fn_chr[2:length(fn_chr)])
#                                                      close_open_sinks()
#                                                    })
#                                    }
#                      )
#                  })
# }


write_fns_to_split_destns(pkg_depcy_ls,
                          pkg_1_core_fns_chr = c("add_utility_predn_to_ds",
                                                 "make_fake_ts_data",
                                                 "write_all_alg_outps",
                                                 "write_rprt",
                                                 "write_shareable_mdls",
                                                 "write_ts_mdls_from_alg_outp"),
                          original_pkg_nm_1L_chr = get_dev_pkg_nm(),
                          pkg_1_nm_1L_chr = "TTU",
                          pkg_2_nm_1L_chr = "map2aqol",
                          tmp_dir_path_1L_chr = "data-raw/pkg_migration",
                          path_to_fns_dir_1L_chr = "data-raw/fns")

