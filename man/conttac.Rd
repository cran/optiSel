\name{conttac}
\Rdversion{1.1}
\alias{conttac}
\title{Contributions To Age Cohorts
}
\description{
Computes genetic contributions from other breeds to all age cohorts
}
\usage{
conttac(cont, cohort, use=rep(TRUE,length(cohort)), mincont=0.05)}

\arguments{
\item{cont}{Data frame containing the genetic contributions of several ancestors or breeds to all individuals. This is typically the output of function \link{breedcont}.}
\item{cohort}{A numeric vector indicating for every individual the age cohort to which it belongs (typically year of birth).}
\item{use}{A logical vector indicating for every individual whether it should be included in an age cohort (typically \code{TRUE} for individuals belonging to the breed of interest).}
\item{mincont}{Contributions of breeeds with average contribution smaller than \code{mincont} will be summarized in one row}
}

\details{
The genetic contributions from other breeds to all age cohorts are computed. The genetic contribution from a breed is the fraction of genes in the gene pool originating from the respective breed.
}

\value{
Data frame containing the genetic contribution from every breed to every age cohort.
}


\examples{
data(ExamplePed)
Pedig      <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
cont       <- breedcont(Pedig, thisBreed="Hinterwaelder")
contByYear <- conttac(cont, cohort=Pedig$Born, use=Pedig$Breed=="Hinterwaelder", mincont=0.04)
round(contByYear,2)

barplot(contByYear, ylim=c(0,1), col=1:10, ylab="genetic contribution",
        legend=TRUE, args.legend=list(x="topleft",cex=0.6))

}



\author{Robin Wellmann}
