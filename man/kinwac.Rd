\name{kinwac}
\Rdversion{1.1}
\alias{kinwac}
\title{Kinships With Age Cohorts}

\description{
Computes for every individual it's mean kinship with the age cohort to which it belongs (excluding the kinship with itself). Different notions of kinship can be used simultaneously.
}
\usage{
kinwac(K, Pedig, use=NULL)}

\arguments{
\item{K}{A named list containing kinship matrices, typically the output of function \link{kin}.}
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam, and column \code{Born} containing for every individual the age cohort to which it belongs (typically year of birth).}
\item{use}{A logical vector indicating for every individual whether it should be used for computing the mean kinships (typically \code{TRUE} for individuals belonging to the breed of interest). The default means that all individuals should be used.}
}

\details{
The function computes for every individual it's mean kinship with the age cohort to which it belongs. This can be computed for different notions of kinship simultaneously. If an individual has high kinship with it's age cohort then it should not be extensively used for breeding if the aim is to conserve genetic diversity.
}

\value{
A data frame with column names \code{c("cohort", "used", "I", names(K))} containing for every individual and every kinship matrix provided in \code{K} the average kinship of the individual with the birth cohort to which it belongs.


The data frame has class attribute "kinwac".
}


\examples{
data(ExamplePed)
Pedig   <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
Kinships<- kin(Pedig, method=c("A", "D"), thisBreed="Hinterwaelder") 
Kin     <- kinwac(Kinships, Pedig=Pedig, use=Pedig$Breed=="Hinterwaelder")
Kin[1001:1010,]

summary(Kin, tlim=c(1970,2005), histNe=150, base=1800)
}



\author{Robin Wellmann}
