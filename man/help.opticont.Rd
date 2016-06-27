\name{help.opticont}
\Rdversion{1.1}
\alias{help.opticont}
\title{Available objective functions and constraints.
}
\description{
Lists available objective functions and constraints for function \link{opticont}.
}
\usage{
help.opticont(K, phen)}

\arguments{
\item{K}{A named list containing one or more NxN Kinship matrices for the N breeding animals, e.g. the output of function \link{kin}.}
\item{phen}{Data frame with N rows and various columns, with animal IDs in column 1 and sex in column 2. The other columns may contain breeding values or migrant contributions.}
}

\details{
Lists available objective functions and constraints for function \link{opticont}.
}

\examples{
#see ?opticont
}

\author{Robin Wellmann}
