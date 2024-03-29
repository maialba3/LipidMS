% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idFunctionsPos.R
\name{idPOS}
\alias{idPOS}
\title{Lipids annotation for ESI+}
\usage{
idPOS(
  msobject,
  ppm_precursor = 5,
  ppm_products = 10,
  rttol = 5,
  coelCutoff = 0.8,
  lipidClasses = c("MG", "LPC", "LPE", "PC", "PCo", "PCp", "PE", "PEo", "PEp", "PG",
    "PI", "Sph", "SphP", "Cer", "AcylCer", "CerP", "SM", "Carnitines", "CE", "DG", "TG"),
  dbs,
  verbose = TRUE
)
}
\arguments{
\item{msobject}{an msobject returned by \link{dataProcessing}.}

\item{ppm_precursor}{mass tolerance for precursor ions. By default, 5 ppm.}

\item{ppm_products}{mass tolerance for product ions. By default, 10 ppm.}

\item{rttol}{total rt window for coelution between precursor and product
ions. By default, 5 seconds.}

\item{coelCutoff}{coelution score threshold between parent and fragment ions.
Only applied if rawData info is supplied. By default, 0.8.}

\item{lipidClasses}{classes of interest to run the identification functions.}

\item{dbs}{list of data bases required for annotation. By default, dbs
contains the required data frames based on the default fragmentation rules.
If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
and \link{assignDB}.}

\item{verbose}{print information messages.}
}
\value{
annotated msobject (list with several elements). The results element
is a data frame that shows: ID, lipid class, CDB (total number of carbons
and double bounds), FA composition (specific chains composition if it has
been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
confidenceLevel (Subclass, FA level, where chains are known but not their
positions, or FA position level), peakID, and Score (parent-fragment coelution 
score mean in DIA data or relative sum intensity in DDA of all fragments used 
for the identification); and the annotatedPeaklist element shows the original 
MS1 peaklist with the annotations on it.
}
\description{
Lipids annotation based on fragmentation patterns for LC-MS/MS DIA or DDA data
acquired in positive mode. This function compiles all functions written for
ESI+ annotations.
}
\examples{
\dontrun{
msobject <- idPOS(msobject)
}

}
\author{
M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
}
