% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twophase.R
\name{evalTwoPhase}
\alias{evalTwoPhase}
\title{evalTwoPhase functions}
\usage{
evalTwoPhase(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q, p, method)
}
\arguments{
\item{Outcome}{vector of length n}

\item{X}{matrix of size n x t (t standard predictors included in the model)}

\item{Z}{matrix of size n x m (m biomarkers included in the model); some subjects have NA values}

\item{Stratum}{vector of length n}

\item{Phase_ID}{vector of length n}

\item{namesX}{names of X; vector of length t}

\item{namesZ}{names of Z; vector of length m}

\item{q}{threshold value for TPR(q)}

\item{p}{threshold value for FPR(p)}

\item{method}{"Benchmark", "ML"}

\item{numBoot}{number of bootstrap samples}
}
\value{
calculate TPR, FPR, and AUC, with two-phase estimators implemented via function tps() in R package "osDesign" (Haneuse and others, 2011).
}
\description{
evalTwoPhase functions
}
\examples{
est_benchmark <- evalTwoPhase(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q=0.2, p=0.2, method="Benchmark")
}
