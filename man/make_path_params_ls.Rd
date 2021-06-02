% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_make.R
\name{make_path_params_ls}
\alias{make_path_params_ls}
\title{Make path params}
\usage{
make_path_params_ls(
  path_to_data_from_top_level_chr,
  path_from_top_level_1L_chr = NULL,
  path_to_current_1L_chr = NULL
)
}
\arguments{
\item{path_to_data_from_top_level_chr}{Path to data from top level (a character vector)}

\item{path_from_top_level_1L_chr}{Path from top level (a character vector of length one), Default: NULL}

\item{path_to_current_1L_chr}{Path to current (a character vector of length one), Default: NULL}
}
\value{
Path params (a list)
}
\description{
make_path_params_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make path params list. The function returns Path params (a list).
}
\keyword{internal}