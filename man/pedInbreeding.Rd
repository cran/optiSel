\name{pedInbreeding}
\Rdversion{1.0}
\alias{pedInbreeding}
\title{Calculates Pedigree Based Inbreeding 
}
\description{                                                 
Calculates Pedigree Based Inbreeding 
}
\usage{
pedInbreeding(Pedig)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam. Missing parents are coded as NA. Consider preparing the pedigree with function \link{prePed} before.}
}

\details{
Computation of pedigree based inbreeding.
This function is a wrapper function for \code{\link[pedigree:calcInbreeding]{pedigree}} from package  \code{\link[pedigree:pedigree-package]{pedigree}} and is compatible with the data.table package. 
}

\value{
A data frame or data table with column \code{Indiv} containing the individual IDs and column \code{Inbr} containing the inbreeding coefficients.
}


\examples{
data(PedigWithErrors)
data(Phen)
keep  <- Phen$Indiv 
Pedig <- prePed(PedigWithErrors, keep=keep)
Res   <- pedInbreeding(Pedig)
mean(Res$Inbr[Res$Indiv \%in\% keep])
#[1] 0.01943394
}



\author{Robin Wellmann}
