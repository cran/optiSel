\name{summary.kinMatrices}
\Rdversion{1.1}
\alias{summary.kinMatrices}
\title{Calculates Genetic Parameters for Age Cohorts}

\description{
Computes for every age cohort several genetic parameters. These are average kinships, gene diversity at native loci, native effective size, and native genome equivalents.
}
\usage{
\method{summary}{kinMatrices}(object, Pedig, tlim=NULL, histNe=NULL, base=NULL, df=4, ...)}

\arguments{
\item{object}{This is typically the output of function \link{kinlist}. A named list containing kinship matrices of the individuals to be used for computing mean kinships. Typically individuals from different age cohorts are included, but only from the breed of interest.}
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam, and column \code{Born} containing for every individual the age cohort to which it belongs (typically year of birth).}
\item{tlim}{Numeric vector with 2 components giving the time span for which genetic parameters are to be computed.}
\item{histNe}{The historic effective size of the population assumed for the time before \code{tlim[1]}. The default is 3 times the average native effective size within the time span of interest. This parameter is needed only for computing the conditional gene diversity and the native genome equivalents.}
\item{base}{The base year in which individuals are assumed to be unrelated. The default is 25 generations before  \code{tlim[1]}. This parameter is needed only for computing the conditional gene diversity and the native genome equivalents.}
\item{df}{Smoothing parameter used for computing the native effective size. The default is \code{df=4}.}
\item{...}{further arguments passed to or from other methods}
}

\details{
This function computes for every age cohort several genetic parameters. These may include 
\tabular{ll}{
\code{pedIBD}     \tab The average pedigree based kinship within the age cohort.\cr
\code{pedIBDorM}  \tab The probability that 2 alleles chosen from the age cohort are \cr
                  \tab IBD or that at least one of them originates from a migrant.\cr
\code{pedIBDorMM} \tab The probability that 2 alleles chosen from the age cohort are \cr
                  \tab IBD or that both of them are from migrants.\cr
\code{pedIBDofN}    \tab The conditional probability that 2 alleles chosen from the \cr
                  \tab age cohort are IBD, given that both are from native founders.\cr
\code{condGD}   \tab The diversity at native alleles, i.e. the conditional probability \cr
                \tab that 2 alleles chosen from the age cohort are not IBD, given that \cr
                \tab both are from native founders. Individuals born in the base-year \cr
                \tab are assumed to be unrelated.\cr 
\code{Ne}   \tab The native effective size, quantifying how fast the diversity at native\cr
                \tab alleles is decreasing. The diversity may increase for a short time \cr
                \tab span, in which case the estimate would be NA. Use a smaller value \cr 
                \tab for parameter \code{df} to get a smooth estimate.\cr 
\code{NGE}   \tab The native genome equivalents. The parameter estimates the  number of \cr
             \tab unrelated individuals that would be needed to establish a hypothetical \cr 
             \tab new population that has the same \code{condGD} as the population under \cr                  \tab study. Individuals born in the base-year are assumed to be unrelated.\cr 
}
}

\value{
A data frame providing for every age cohort (some of) the genetic parameters mentioned above. 
}

\examples{
data(ExamplePed)
Pedig    <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
Kinships <- kinlist(pedIBD=pedIBD(Pedig), pedIBDatN=pedIBDatN(Pedig, thisBreed="Hinterwaelder"))
summary(Kinships, Pedig, tlim=c(1970,2005), histNe=150, base=1800)
}



\author{Robin Wellmann}
