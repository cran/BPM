% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayPM.R
\name{BayPM}
\alias{BayPM}
\title{Bayesian Purity Model (BPM) Main functions.}
\usage{
BayPM(betaValue, TOPK = 500, tumorNum = NULL, filterProbes = FALSE,
  userProbes = NULL)
}
\arguments{
\item{betaValue}{A matrix,TCGA methlation array data. Each row: loci,
Tumor1,Tumor2,...,Normal1,Nomral2,...}

\item{TOPK}{A number. Number of DMCs/nonDMCs selected}

\item{tumorNum}{The number of tumor samples.
if NULL, the default number is half of column number of dataset.}

\item{filterProbes}{Logistic. defalut is FALSE. The code use all probes in betaValue.
If TRUE,
you can use default good probes provided in our code.
you can also provide your good probes in userProbes.}

\item{userProbes}{A number list. The row numbers in betaValue.
These rows are considered as good probes.}
}
\value{
tumor purity estimation of tumor samples
}
\description{
Bayesian Purity Model (BPM) Main functions.
}
\examples{
### need to install package "limma"
### source("https://bioconductor.org/biocLite.R");biocLite("limma");
BayPM(simUCEC,20,2);

}
