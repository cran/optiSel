\name{prePed}
\Rdversion{1.1}
\alias{prePed}
\title{Prepares a pedigree
}
\description{
Prepares a pedigree by sorting and adding founders and pruning the pedigree.
}
\usage{
prePed(Pedig, keep=NULL, thisBreed=NA, lastNative=NA, addNum=FALSE, I=0)}

\arguments{
\item{Pedig}{Data frame containing the pedigree where the first 3 columns correspond to: Individual ID,  Sire, and Dam. More columns can be passed in the \code{Pedig} argument including columns named \code{Sex}, \code{Breed} (with breed names), and \code{Born} (with years of birth). Missing parents are coded as NA, 0, or "0".}
\item{keep}{A vector with IDs of individuals or \code{NULL}. If this parameter is not \code{NULL}, then only these individuals and their ancestors will be kept in the pedigree.}
\item{thisBreed}{Name of the breed.}
\item{lastNative}{Last year of birth for which individuals with unknown pedigree are considered native.}
\item{addNum}{If \code{TRUE}, then columns with IDs of individuals, sires, and dams in integer form will be added.
}
\item{I}{offset}
}

\details{

This function takes a pedigree, adds missing founders, and sorts the pedigree. If parameter \code{keep} contains IDs of individuals then only these individuals and their ancestors will be kept in the pedigree.

If the pedigree contains loops, then the loops will be broken by setting the parents of one animal in each loop to \code{NA}. 

If the pedigree contains columnn \code{Sex} then the sexes will be recoded as 1 for males and 2 for females. Missing sexes will be determined from pedigree structure if possible. 

If the pedigree contains column \code{Breed} then for ancestors with missing breed the breed name is estimated. If parameter \code{lastNative} is not \code{NA} then for each animal with one missing parent an imaginary founder is added to the pedigree in order to enable classifying the breed names of all founders as follows: 
In general animals with missing breed are assumed to have the same breed as most of their offspring. But there is one exception: For founders belonging to \code{thisBreed} who are born after \code{lastNative} the breed name will be set to \code{"unknown"}. Moreover for founders from \code{thisBreed} with unknown year of birth the breed name will be set to \code{"unknown"} if all their descendants are born after \code{lastNative+I}. 

 

}

\value{
A data frame containing the pedigree.
}


\examples{
data(PedigWithErrors)
Pedig <- prePed(PedigWithErrors)
}



\author{Robin Wellmann}
