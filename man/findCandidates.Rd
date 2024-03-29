% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subfunctionsIdentification.R
\name{findCandidates}
\alias{findCandidates}
\title{Search of lipid candidates of a certain class}
\usage{
findCandidates(
  MS1,
  db,
  ppm,
  rt,
  adducts,
  rttol = 3,
  dbs,
  rawData = data.frame(),
  coelCutoff = 0
)
}
\arguments{
\item{MS1}{peaklist of the MS function. Data frame with 3 columns: mz, RT (in seconds) and int (intensity).}

\item{db}{database (i.e. pcdb, dgdb, etc.). Data frame with at least 2 columns: Mass (exact mass) and total (total number of carbons and double bound of the FA chains, i.e. "34:1").}

\item{ppm}{m/z tolerance in ppm.}

\item{rt}{rt range where the function will look for candidates. By default,
it will search within all RT range in MS1.}

\item{adducts}{character vector containing the expected adducts to search
for (i.e. "M+H", "M+Na", "M-H", etc.). See details.}

\item{rttol}{rt tolerance in seconds to match adducts.}

\item{dbs}{list of data bases required for the annotation. By default, dbs
contains the required data frames based on the default fragmentation rules.
If these rules are modified, dbs may need to be changed. If data bases have
been customized using \link{createLipidDB}, they also have to be modified
here.}

\item{rawData}{raw scans data. Output of \link{dataProcessing} function (MS1$rawData).}

\item{coelCutoff}{coelution score threshold between parent and fragment ions.
Only applied if rawData info is supplied.}
}
\value{
Data frame with the found candidates. It contains 6 columns: mz,
RT, int (from the peaklist data.frame), ppms, cb (total number of carbons and
double bounds of the FA chains) and adducts.
}
\description{
Search of lipid candidates from a peaklist based on a set of expected
adducts.
}
\details{
\link{findCandidates} looks for matches between the m/z of the
MS1 peaklist and the expected m/z of the candidates in the database for each
adduct. If several adducts are expected, results are combined.

Adducts allowed are contained in adductsTable data frame, which can be
modified if required (see \link{adductsTable}).
}
\author{
M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
}
