\name{completeness}
\Rdversion{1.1}
\alias{completeness}
\title{Computes completeness of the pedigree}

\description{
Computes completeness of the pedigree for the specified individuals.
}
\usage{
completeness(Pedig, keep, genNo=-1, type="MEAN")}

\arguments{
\item{Pedig}{Data frame containing the pedigree where the first 3 columns correspond to: Individual ID,  Sire, and Dam. More columns can be passed in the \code{Pedig} argument including a column named \code{Sex}.}
\item{keep}{Vector with IDs of individuals to be included. \code{"0"} will select all individuals without children.}
\item{genNo}{Vector of generation numbers at which completeness should be calculated. Default is -1, which calculates completeness at each generation.}
\item{type}{If \code{type="MEAN"} (default), mean completeness over all specified individuals is calculated. If \code{type="IND"}, completeness is calculated for each specified individual.}
}

\details{
The function computes the completeness of the pedigree for the specified individuals. This is a wrapper function for function \link[GENLIB]{gen.completeness} from package \link[GENLIB:GenlibR-package]{GENLIB} that allows computation of pedigree completeness with little R code.
}

\value{
When type is \code{"MEAN"} a vector is returned with the completeness for each generation at which the completeness is calculated.

When type is \code{"IND"}, a data frame is returned with one row for each individual included in vector \code{keep}. The first columns contain the pedigree completeness for each generation at which the completeness is calculated. The columns of data frame \code{Pedig} are appended thereafter.

}


\examples{

#Computes the pedigree completeness of the parents of
#Hinterwald cattle born in 2006 and 2007.

data(PedigWithErrors)
Pedig <- prePed(PedigWithErrors)

inCohort <- (Pedig$Born \%in\% (2006:2007)) & Pedig$Breed=="Hinterwaelder"

compl<-completeness(Pedig, keep=Pedig$Sire[inCohort], type="IND")
compl[1:20,]
matplot(t(compl[1:20,1:10]),type="l",xlab="Generation", ylab="Completeness")
 
compl.m<-completeness(Pedig, keep=Pedig$Sire[inCohort])
compl.f<-completeness(Pedig, keep=Pedig$Dam[inCohort])

plot(compl.m, type="b", ylim=c(0,100), ylab="completeness", xlab="Generation")
lines(compl.f, col="red", type="b")
legend(12,95,c("Males", "Females"), col=c("black","red"), lty=1, cex=0.8)
}

\references{Cazes P, Cazes MH. (1996) Comment mesurer la profondeur genealogique d'une ascendance? Population (French Ed) 51:117-140.}

\author{Robin Wellmann}
