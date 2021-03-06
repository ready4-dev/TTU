% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_transform.R
\name{transform_tb_to_mdl_inp}
\alias{transform_tb_to_mdl_inp}
\title{Transform tibble to model input}
\usage{
transform_tb_to_mdl_inp(
  data_tb,
  depnt_var_nm_1L_chr = "utl_total_w",
  predr_vars_nms_chr,
  id_var_nm_1L_chr = "fkClientID",
  round_var_nm_1L_chr = "round",
  round_bl_val_1L_chr = "Baseline",
  drop_all_msng_1L_lgl = T,
  scaling_fctr_dbl = 0.01,
  tfmn_1L_chr = "NTF",
  ungroup_1L_lgl = F
)
}
\arguments{
\item{data_tb}{Data (a tibble)}

\item{depnt_var_nm_1L_chr}{Dependent variable name (a character vector of length one), Default: 'utl_total_w'}

\item{predr_vars_nms_chr}{Predictor variables names (a character vector)}

\item{id_var_nm_1L_chr}{Identity variable name (a character vector of length one), Default: 'fkClientID'}

\item{round_var_nm_1L_chr}{Round variable name (a character vector of length one), Default: 'round'}

\item{round_bl_val_1L_chr}{Round baseline value (a character vector of length one), Default: 'Baseline'}

\item{drop_all_msng_1L_lgl}{Drop all missing (a logical vector of length one), Default: T}

\item{scaling_fctr_dbl}{Scaling factor (a double vector), Default: 0.01}

\item{tfmn_1L_chr}{Transformation (a character vector of length one), Default: 'NTF'}

\item{ungroup_1L_lgl}{Ungroup (a logical vector of length one), Default: F}
}
\value{
Transformed for model input (a tibble)
}
\description{
transform_tb_to_mdl_inp() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform tibble to model input. Function argument data_tb specifies the object to be updated. Argument depnt_var_nm_1L_chr provides the object to be updated. The function returns Transformed for model input (a tibble).
}
