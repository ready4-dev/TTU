% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_add.R
\name{add_utl_predn_to_new_ds}
\alias{add_utl_predn_to_new_ds}
\title{Add utility prediction to new dataset}
\usage{
add_utl_predn_to_new_ds(
  data_tb,
  ingredients_ls,
  mdl_nm_1L_chr,
  analysis_1L_chr = NULL,
  deterministic_1L_lgl = T,
  force_min_max_1L_lgl = T,
  id_var_nm_1L_chr = NULL,
  model_mdl = NULL,
  new_data_is_1L_chr = "Simulated",
  predr_vars_nms_chr = NULL,
  round_var_nm_1L_chr = "Timepoint",
  round_bl_val_1L_chr = "BL",
  utl_cls_fn = NULL,
  utl_var_nm_1L_chr = NULL
)
}
\arguments{
\item{data_tb}{Data (a tibble)}

\item{ingredients_ls}{Ingredients (a list)}

\item{mdl_nm_1L_chr}{Model name (a character vector of length one)}

\item{analysis_1L_chr}{Analysis (a character vector of length one), Default: NULL}

\item{deterministic_1L_lgl}{Deterministic (a logical vector of length one), Default: T}

\item{force_min_max_1L_lgl}{Force minimum maximum (a logical vector of length one), Default: T}

\item{id_var_nm_1L_chr}{Identity variable name (a character vector of length one), Default: NULL}

\item{model_mdl}{Model (a model), Default: NULL}

\item{new_data_is_1L_chr}{New data is (a character vector of length one), Default: 'Simulated'}

\item{predr_vars_nms_chr}{Predictor variables names (a character vector), Default: NULL}

\item{round_var_nm_1L_chr}{Round variable name (a character vector of length one), Default: 'Timepoint'}

\item{round_bl_val_1L_chr}{Round baseline value (a character vector of length one), Default: 'BL'}

\item{utl_cls_fn}{Utility class (a function), Default: NULL}

\item{utl_var_nm_1L_chr}{Utility variable name (a character vector of length one), Default: NULL}
}
\value{
Updated (a tibble)
}
\description{
add_utl_predn_to_new_ds() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add utility prediction to new dataset. Function argument data_tb specifies the object to be updated. The function returns Updated (a tibble).
}
\keyword{internal}
