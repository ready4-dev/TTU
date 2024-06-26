% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mthd_author.R
\name{author-TTUReports}
\alias{author-TTUReports}
\alias{author,TTUReports-method}
\alias{author-TTUProject}
\alias{author,TTUProject-method}
\title{Author and save files}
\usage{
\S4method{author}{TTUReports}(
  x,
  args_ls = NULL,
  consent_1L_chr = "",
  depnt_var_desc_1L_chr = NA_character_,
  depnt_var_min_val_1L_dbl = numeric(0),
  download_tmpl_1L_lgl = T,
  fl_type_1L_chr = ".eps",
  timepoint_new_nms_chr = NA_character_,
  type_1L_chr = "Report",
  what_1L_chr = NA_character_,
  x_labels_chr = character(0),
  ...
)

\S4method{author}{TTUProject}(
  x,
  consent_1L_chr = "",
  custom_args_ls = NULL,
  custom_fn = NULL,
  depnt_var_min_val_1L_dbl = numeric(0),
  digits_1L_int = 2L,
  download_tmpl_1L_lgl = T,
  fl_nm_1L_chr = "TTUProject",
  items_as_domains_1L_lgl = F,
  supplement_fl_nm_1L_chr = "TA_PDF",
  timepoint_new_nms_chr = NA_character_,
  type_1L_chr = "auto",
  what_1L_chr = "default",
  ...
)
}
\arguments{
\item{x}{An object of class TTUProject}

\item{args_ls}{Arguments (a list), Default: NULL}

\item{consent_1L_chr}{Consent (a character vector of length one), Default: ''}

\item{depnt_var_desc_1L_chr}{Dependent variable description (a character vector of length one), Default: 'NA'}

\item{depnt_var_min_val_1L_dbl}{Dependent variable minimum value (a double vector of length one), Default: numeric(0)}

\item{download_tmpl_1L_lgl}{Download template (a logical vector of length one), Default: T}

\item{fl_type_1L_chr}{File type (a character vector of length one), Default: '.eps'}

\item{timepoint_new_nms_chr}{Timepoint new names (a character vector), Default: 'NA'}

\item{type_1L_chr}{Type (a character vector of length one), Default: 'auto'}

\item{what_1L_chr}{What (a character vector of length one), Default: 'default'}

\item{x_labels_chr}{X labels (a character vector), Default: character(0)}

\item{...}{Additional arguments}

\item{custom_args_ls}{Custom arguments (a list), Default: NULL}

\item{custom_fn}{Custom (a function), Default: NULL}

\item{digits_1L_int}{Digits (an integer vector of length one), Default: 2}

\item{fl_nm_1L_chr}{File name (a character vector of length one), Default: 'TTUProject'}

\item{items_as_domains_1L_lgl}{Items as domains (a logical vector of length one), Default: F}

\item{supplement_fl_nm_1L_chr}{Supplement file name (a character vector of length one), Default: 'TA_PDF'}
}
\value{
No return value, called for side effects.

x (An object of class TTUProject)
}
\description{
author method applied to TTUReports

author method applied to TTUProject
}
