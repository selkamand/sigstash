% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_load.R
\name{sig_load}
\alias{sig_load}
\title{Load a Signature Collection}
\usage{
sig_load(dataset, return_df = FALSE)
}
\arguments{
\item{dataset}{a valid signature collection to load. Run \code{\link[=sig_available]{sig_available()}} to list available datasets.}

\item{return_df}{return a data.frame instead of a signature collection list (boolean)}
}
\value{
The specified signature collection, either as a list conforming to sigverse signature collection format, or, if \code{return_df = TRUE} as a data.frame.
}
\description{
Load a Signature Collection
}
\examples{

# List available datasets
sig_available()

# Load available datasets
sig_load("COSMIC_v3.3.1_SBS_GRCh38")

}
