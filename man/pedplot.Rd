\name{pedplot}
\Rdversion{1.1}
\alias{pedplot}
\title{Plot a pedigree}
\description{
Plot a pedigree
}
\usage{
pedplot(Pedig, affected=NULL, status=NULL, label="Indiv", ...)}

\arguments{
\item{Pedig}{Data frame containing the pedigree where the first 3 columns correspond to: Individual ID,  Sire, and Dam. Use \link{subPed} to ensure that the pedigree is in the correct format.}
\item{affected}{
Logical vector indicating for each individual if its symbol should be plotted in colour. The default \code{NULL} means that the individuals in column \code{keep} of data frame \code{Pedig} are plotted in colour (if present).
}
\item{status}{Logical vector indicating for each individual if its symbol in the plot should be crossed out. The default \code{NULL} means that animals from other breeds than those plotted in colour are crossed out.}
\item{label}{Character vector containing the columns of data frame \code{Pedig} to be used as labels.}
\item{...}{Options passed to the underlying function \link[kinship2]{plot.pedigree} from package \code{kinship2}.}
}

\details{
This function plots a pedigree. If data frame \code{Pedig} has logical column \code{keep} then the default values mean that the symbols of these animals are plotted in color and for animals from other breeds the symbol is crossed out.
}

\value{
An invisible list returned by the underlying function \link[kinship2]{plot.pedigree} from package \code{kinship2}.
}


\examples{
data(PedigWithErrors)

sPed  <- subPed(PedigWithErrors, keep="276000891974272", prevGen=3, succGen=2)
pedplot(sPed, mar=c(2,4,2,4), label=c("Indiv", "Born", "Breed"), cex=0.7)
}



\author{Robin Wellmann}
