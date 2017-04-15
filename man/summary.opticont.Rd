\name{summary.opticont}
\Rdversion{1.1}
\alias{summary.opticont}
\title{Summary Statistics and Check of Validity for Optimum Contributions.
}
\description{
Computation of genetic parameters and check of validity for optimum genetic contributions of the selection candidates.
}
\usage{
\method{summary}{opticont}(object, ...)}

\arguments{
\item{object}{An object from class \code{opticont}, which is usually the output of function \link{opticont}.}
\item{...}{further arguments passed to or from other methods}
}

\details{
Computation of genetic parameters and check of validity for optimum genetic contributions of the selection candidates.
}

\value{
A data frame with one row and the following columns is returned:
\tabular{ll}{
\code{VarName} \tab Name of parameter x \cr
\code{method} \tab Method used for optimization \cr
\code{obj.fun} \tab Value of the objective function \cr
\code{valid} \tab TRUE if the side constraints are fulfilled by the optimum contributions. \cr
\code{KIN} \tab The mean kinship in the offspring (one component for each kinship matrix included\cr
           \tab in the call of function \link{opticont}).\cr
\code{ub.KIN} \tab The upper bound for the mean kinship in the offspring (one component for each\cr
              \tab kinship matrix included in the call of function \link{opticont}). \cr
\code{Div.KIN} \tab The mean genetic diversity (\code{=1-KIN}) in the offspring (one component \cr
              \tab for each kinship matrix included in the call of function \link{opticont}).\cr
\code{ubM} \tab Upper bound for the contribution of males \cr
\code{ubF} \tab Upper bound for the contribution of females \cr
\code{ContMales} \tab The total contribution of breeding males to the offspring (should be 0.5) \cr
\code{ContFemales} \tab The total contribution of breeding females to the offspring (should be 0.5) \cr
\code{minCont} \tab The minimum contribution of a breeding animal \cr
\code{maxContMale} \tab The maximum contribution of a breeding male \cr
\code{maxContFemale} \tab The maximum contribution of a breeding female \cr
\code{solver} \tab Algorithm used for optimization\cr  
\code{lb.VAR} \tab The lower bound for the mean value of \code{VAR} in the offspring (one \cr  
              \tab component for each numeric column \code{VAR} of data frame \code{phen} \cr
              \tab included in the call of function \link{opticont}). \cr
\code{meanVAR} \tab The mean value of \code{VAR} in the offspring (one component for each \cr 
              \tab numeric column \code{VAR} of data frame \code{phen} included in the call of \cr
              \tab function \link{opticont}).\cr
\code{ub.VAR} \tab The upper bound for the mean value of \code{VAR} in the offspring (one \cr 
              \tab component for each numeric column \code{VAR} of data frame \code{phen} \cr
              \tab included in the call of function \link{opticont}).  \cr
}

}


\examples{

#see ?opticont
}



\author{Robin Wellmann}
