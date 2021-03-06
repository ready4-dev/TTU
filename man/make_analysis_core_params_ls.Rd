% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_make.R
\name{make_analysis_core_params_ls}
\alias{make_analysis_core_params_ls}
\title{Make analysis core params}
\usage{
make_analysis_core_params_ls(
  ds_descvs_ls,
  mdl_smry_ls,
  output_format_ls = make_output_format_ls(),
  predictors_lup,
  control_ls = NULL,
  iters_1L_int = 4000L,
  prior_ls = NULL,
  seed_1L_int = 12345,
  use_fake_data_1L_lgl = F
)
}
\arguments{
\item{ds_descvs_ls}{Dataset descriptives (a list)}

\item{mdl_smry_ls}{Model summary (a list)}

\item{output_format_ls}{Output format (a list), Default: make_output_format_ls()}

\item{predictors_lup}{Predictors (a lookup table)}

\item{control_ls}{Control (a list), Default: NULL}

\item{iters_1L_int}{Iterations (an integer vector of length one), Default: 4000}

\item{prior_ls}{Prior (a list), Default: NULL}

\item{seed_1L_int}{Seed (an integer vector of length one), Default: 12345}

\item{use_fake_data_1L_lgl}{Use fake data (a logical vector of length one), Default: F}
}
\value{
Analysis core params (a list)
}
\description{
make_analysis_core_params_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make analysis core params list. The function returns Analysis core params (a list).
}
\keyword{internal}
