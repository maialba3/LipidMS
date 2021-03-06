% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{cbs}
\alias{cbs}
\title{Total number of carbons and double bounds}
\usage{
cbs(mz, ppm, db, charge = 0)
}
\arguments{
\item{mz}{mz value to be matched}

\item{ppm}{mass error tolerance}

\item{db}{database}

\item{charge}{numeric value indicating the charge of the ion}
}
\value{
Character value or vector indicating structural information
(carbons:bounds)
}
\description{
This function matches mz values with neutral masses from a dataframe which
links masses and structures (carbons and double bounds) and extracts the
structural information. It is used by identification functions to look for
the structure of the previously chosen candidates.
}
\author{
M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
}
\keyword{internal}
