\name{genecont}
\Rdversion{1.1}
\alias{genecont}
\title{Calculation of genetic contributions. 
}
\description{
Genetic contributions of all ancestors to each individual
}
\usage{
genecont(ID, Sire, Dam, NAncestors=NA)
}


\arguments{
\item{ID}{Vector with the IDs of the individuals, ordered by date of birth.}
\item{Sire}{Vector with the IDs of the sires.}
\item{Dam}{Vector with the IDs of the dams.}
\item{NAncestors}{If this parameter is not \code{NA} then only the genetic contributions of the first \code{NAncestors} individuals in the pedigree to all other individuals are computed.}
}



\details{
This function calculates genetic contributions of all ancestors to each individual.
}

\value{
Lower triangular matrix with genetic contributions for each pair of individuals. Column i contains the genetic contribution of ancestor i to all individuals 
}


\examples{

data(ExamplePed)
Pedig <- prePed(ExamplePed)
cont  <- genecont(ID=Pedig$Indiv, Sire=Pedig$Sire, Dam=Pedig$Dam)

plot(Pedig$Born, cont[,"276000803611144"], pch=18, ylim=c(0,1))
Pedig["276000803611144",]

}



\author{Robin Wellmann}
