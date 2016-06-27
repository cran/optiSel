\name{summary.kinwac}
\Rdversion{1.1}
\alias{summary.kinwac}
\title{Computes genetic parameters for all age cohorts}

\description{
Computes for every age cohort several genetic parameters. These are average kinships, genetic diversities, conditional gene diversity, native effective size, and native genome equivalents.
}
\usage{
\method{summary}{kinwac}(object, tlim=NULL, histNe=NULL, base=NULL, df=4, ...)}

\arguments{
\item{object}{This is typically the output of function \link{kinwac}. It is a data frame with column names "cohort", "used", "I", and the names for different notions of kinship. Column \code{cohort} is numeric and indicates for every individual the age cohort to which it belongs (typically year of birth). Column \code{used} is a logical vector and indicates for every individual whether it should be used for computing the mean kinships (typically \code{TRUE} for individuals belonging to the breed of interest). Column \code{I} contains for every individual the realized generation interval. The remaining columns correspond to different notions of kinship. Each column contains for every individual the average kinship of the individual with the birth cohort to which it belongs. }
\item{tlim}{A numeric vector with 2 components giving the time span for which genetic parameters are to be computed.}
\item{histNe}{The historic effective size of the population assumed for the time before \code{tlim[1]}. The default is 3 times the average native effective size within the time span of interest. This parameter is needed only for computing the conditional gene diversity and the native genome equivalents.}
\item{base}{The base year in which individuals are assumed to be unrelated. The default is 25 generations before  \code{tlim[1]}. This parameter is needed only for computing the conditional gene diversity and the native genome equivalents.}
\item{df}{Smoothing parameter used for computing the native effective size. The default is \code{df=4}.}
\item{...}{further arguments passed to or from other methods}
}

\details{
This function computes for every age cohort several genetic parameters, see Wellmann et al. (2012) for details. These may include 
\tabular{ll}{
\code{fA}    \tab The average pedigree based kinship within the age cohort.\cr
\code{fB}    \tab The probability that 2 alleles chosen from the age cohort are IBD or that at least one of them originates from a migrant.\cr
\code{fC}    \tab The probability that 2 alleles chosen from the age cohort are IBD or that both of them are from migrants.\cr
\code{fD}    \tab The conditional probability that 2 alleles chosen from the age cohort are IBD, given that both are from native founders.\cr
\code{Div.fA}   \tab The genetic diversity of the age cohort (= 1-\code{fA}). \cr 
\code{Div.fB}   \tab The probability that 2 alleles chosen from the age cohort are not IBD and are from native founders  (= 1-\code{fB}). \cr 
\code{Div.fC}   \tab The probability that 2 alleles chosen from the age cohort are not IBD and at least one of them is from a native founder  (= 1-\code{fC}). \cr 
\code{Div.fD}   \tab The diversity of native alleles, i.e. the conditional probability that 2 alleles chosen from the age cohort are not IBD, given that both are from native founders. Individuals with unknown parents are assumed to be unrelated  (= 1-\code{fD}). \cr 
\code{condGD}   \tab The diversity of native alleles, i.e. the conditional probability that 2 alleles chosen from the age cohort are not IBD, given that both are from native founders. Individuals born in the base-year are assumed to be unrelated.\cr 
\code{Ne}   \tab The native effective size, quantifying how fast the diversity of native alleles is decreasing.\cr 
\code{NGE}   \tab The native genome equivalents. The parameter estimates the  number of unrelated individuals that would be needed to establish a hypothetical new population that has the same \code{condGD} as the population under study. Individuals born in the base-year are assumed to be unrelated.\cr 
}
}

\value{
A list with component \code{Param} which is a data frame providing for every age cohort (some of) the genetic parameters mentioned above. 
}


\examples{
data(ExamplePed)
Pedig   <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
Kinships<- kin(Pedig, method=c("A", "D"), thisBreed="Hinterwaelder") 
Kin     <- kinwac(Kinships, Pedig=Pedig, use=Pedig$Breed=="Hinterwaelder")
sy      <- summary(Kin, tlim=c(1970,2005), histNe=150, base=1800)$Param
}



\author{Robin Wellmann}
