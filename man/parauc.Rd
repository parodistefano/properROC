% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROCProper.R
\name{parauc}
\alias{parauc}
\title{Partial AUC}
\usage{
parauc(sens, onespec, x1, x2)
}
\arguments{
\item{sens}{Array of sensitivity values}

\item{onespec}{Array of 1 - specificity}

\item{x1}{First value of 1 - specificity}

\item{x2}{second First value of 1 - specificity}
}
\value{
pauc, the partial area under the empirical ROC curve
}
\description{
Calculate the partial AUC between two values of 1-specificity
(for internal use only: called by tprc()).
}
