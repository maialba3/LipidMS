% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subfunctionsIdentification.R
\name{chainFrags}
\alias{chainFrags}
\title{Search of chain specific fragments}
\usage{
chainFrags(coelfrags, chainfrags, ppm = 10, candidates, f = NULL, dbs)
}
\arguments{
\item{coelfrags}{coeluting fragments for each candidate. Output of
\link{coelutingFrags}.}

\item{chainfrags}{character vector containing the fragmentation rules for
the chain fragments. If it is an empty vector, chains will be calculated based
on the difference between the precursor and the other chain. See details.}

\item{ppm}{m/z tolerance in ppm.}

\item{candidates}{candidates data frame. If any chain needs to be calculated
based on the difference between the precursor and the other chain, this
argument will be required. Output of \link{chainFrags}.}

\item{f}{known chains. If any chain needs to be calculated
based on the difference between the precursor and the other chain, this
argument will be required. Output of \link{chainFrags}.}

\item{dbs}{list of data bases required for the annotation. By default, dbs
contains the required data frames based on the default fragmentation rules.
If these rules are modified, dbs may need to be changed. If data bases have
been customized using \link{createLipidDB}, they also have to be modified
here.}
}
\value{
List of data frames with the chain fragments found.
}
\description{
Search of specific fragments that inform about the chains structure.
}
\details{
The chainfrags argument must contain the fragmentation rules
which inform about the chains structure. For example, in the case of PG
subclass, the chain in sn1 position is identified by the lysoPG as M-H
resulting from the loss of the FA chain of sn2; and the chain in sn2 position
is identified as the free FA chain as M-H. These two fragments need to be
searched in two different steps: in the fist step we will look for lysoPGs
coeluting with the precursor using chainfrags = c("lysopg_M-H");
then, we will look for FA chains using chainfrags = c("fa_M-H"). This
information can be combined later using \link{combineChains} function.

To indicate the fragments to be searched, the class of lipid is writen
using the same names as the LipidMS databases without the "db" at the end
(i.e. pa, dg, lysopa, mg, CE, etc.), and the adduct has to be indicated as
it appears in the adductsTable, both parts separated by "_". In case some
chain needs to be searched based on a neutral loss, this can be defined using
"NL-" prefix, followed by the database and adduct. If this neutral loss
is employed to find the remaining chain, "cbdiff-" prefix allows to calculate
the difference in carbons and doubles bounds between the precursor and the
building block found. For example, "cbdiff-dg_M+H-H2O" will look for DG as
M+H-H2O and then, it will return the difference between their number of
carbons and double bounds and the ones from the precursor. Otherwise,
"NL-mg_M+H-H2O" will look for fragments coming from the loss of MGs.

In case these fragments identified as losses from the precursors are going to
be employed for the intensity rules, this same prefix has to be added.

If a chain is calculated based on the difference of total number of
carbons and double bounds between the precursor and a previously searched chain,
\code{chainfrags} argument must be a character vector c("") and candidates
data frame and chain fragments list must be provided.
}
\author{
M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
}
