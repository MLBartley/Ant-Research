% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_prep_functions.R
\name{prep_troph_data}
\alias{prep_troph_data}
\title{Data Preparation for Ant Trophallaxis Data File}
\usage{
prep_troph_data(data, hours, delta_t)
}
\arguments{
\item{data}{Ant trophallaxis data file in .csv format, already loaded into R.}

\item{hours}{Total time ants are observed.}

\item{delta_t}{Time increments for (starting time) data to be binned into.}
}
\value{
This function will return one or three (depending on high or low density,
  respectively) file lists with amount of interaction starting within
  binned time chunk (delta_t).
}
\description{
The purpose of this function is to take in .csv files and
apply necessary changes to the data format. Changes include
removing rows/columns of NA, removing redundant entries,
sorting by entrance times, and binning data into smaller
chunks.
}
\examples{
 prep_troph_data(high4, 60)


}
