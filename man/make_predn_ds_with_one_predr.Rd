% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_make.R
\name{make_predn_ds_with_one_predr}
\alias{make_predn_ds_with_one_predr}
\title{Make prediction dataset with one predictor}
\usage{
make_predn_ds_with_one_predr(
  model_mdl,
  depnt_var_nm_1L_chr = "utl_total_w",
  tfmn_1L_chr = "NTF",
  predr_var_nm_1L_chr,
  predr_vals_dbl,
  predn_type_1L_chr = NULL
)
}
\arguments{
\item{model_mdl}{Model (a model)}

\item{depnt_var_nm_1L_chr}{Dependent variable name (a character vector of length one), Default: 'utl_total_w'}

\item{tfmn_1L_chr}{Transformation (a character vector of length one), Default: 'NTF'}

\item{predr_var_nm_1L_chr}{Predictor variable name (a character vector of length one)}

\item{predr_vals_dbl}{Predictor values (a double vector)}

\item{predn_type_1L_chr}{Prediction type (a character vector of length one), Default: NULL}
}
\value{
Prediction dataset (a tibble)
}
\description{
make_predn_ds_with_one_predr() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make prediction dataset with one predictor. The function returns Prediction dataset (a tibble).
}
\keyword{internal}
