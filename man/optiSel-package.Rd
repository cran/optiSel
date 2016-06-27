\name{optiSel-package}
\Rdversion{1.1}
\alias{optiSel-package}
\alias{optiSel}
\docType{package}
\title{Optimum Contribution Selection and Pedigree Analysis}
\description{
Optimum Contribution Selection and Pedigree Analysis
}
\details{
\tabular{ll}{
Package: \tab optiSel\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 27.06.2016\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

The following functions are provided:

\tabular{ll}{
\link{breedcont}   \tab Computes for every individual the genetic contribution from other breeds.\cr 
\link{completeness}   \tab Computes pedigree completeness.\cr 
\link{conttac}   \tab Computes genetic contributions to the age cohorts.\cr
\link{genecont}   \tab  Computes for every individual the genetic contribution from all other individuals.\cr 
\link{help.opticont}   \tab Lists available objective functions and constraints for function \link{opticont}.\cr 
\link{kin}        \tab Calculation of kinship matrices according to different definitions.\cr
\link{kinwac}        \tab Computes for every individual it's mean kinship with the age cohort to which it belongs.\cr
\link{opticont}   \tab Computes optimum genetic contributions of the breeding animals to the next generation.\cr 
\link{pedplot}   \tab Plots a pedigree.\cr 
\link{prePed}   \tab Prepares a pedigree by sorting, adding founders and pruning the pedigree.\cr 
\link{subPed}   \tab Creates a subset of a large pedigree.\cr 
\link{summary.kinwac}       \tab Computes for every age cohort several genetic parameters. These may include average kinships, genetic diversities, conditional gene diversity, native effective size, and native genome equivalents.\cr
\link{summary.opticont}   \tab Computes genetic parameters and checks for validity of optimum genetic contributions.\cr 
}




}


\examples{
#For optimum contribution selection see ?opticont

data(ExamplePed)
Pedig<-prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
Pedig[1:10,]


############################################
# Evaluation of                            #
#    - kinships                            #
#    - genetic diversities                 #
#    - native effective size               #
#    - native genome equivalents           #
############################################

Kinships <- kin(Pedig, method=c("A","B","C","D"), thisBreed="Hinterwaelder")
Kin      <- kinwac(Kinships, Pedig=Pedig, use=Pedig$Breed=="Hinterwaelder")
Kin[1000:1010,]

x  <- summary(Kin, tlim=c(1970,2005), histNe=150, base=1800)
xP <- x$Param

plot(xP$cohort, xP$Ne, type="l", ylim=c(0,150), 
     main="Native Effective Size",ylab="Ne",xlab="")

matplot(xP$cohort, xP[,c("fA", "fB", "fC", "fD")], 
        type="l",ylim=c(0,1),main="Kinships", xlab="Year", ylab="mean Kinship")
mtext("See ?summary.kinwac for definitions of these kinships", cex=0.7)
legend("topleft", legend = c("fA", "fB", "fC", "fD"),lty=1:4,col=1:4,cex=0.6)


matplot(xP$cohort,xP[,c("Div.fA", "Div.fB", "Div.fC", "Div.fD")], 
        type="l",ylim=c(0,1),main="Diversities", xlab="Year", ylab="mean Diversity")
mtext("See ?summary.kinwac for definitions of these diversities", cex=0.7)
legend("topleft", legend = c("Div.fA", "Div.fB", "Div.fC", "Div.fD"),lty=1:4,col=1:4,cex=0.6)

plot(xP$cohort, xP$condGD, type="l", ylim=c(0,1), 
     main="conditional Gene Diversity",ylab="condGD",xlab="")
mtext(paste("Base Year =", x$base, "  historic Ne =", x$histNe), cex=0.7)

plot(xP$cohort,xP$NGE,type="l",main="Native Genome Equivalents", 
     ylab="NGE",xlab="",ylim=c(0,7))
mtext(paste("Base Year =", x$base, "  historic Ne =", x$histNe), cex=0.7)

plot(xP$cohort, xP$fN, type="l",ylim=c(0,1), xlab="")
mtext("Probability that two alleles chosen at random are both native", cex=0.8)

############################################
# Genetic contributions from other breeds  #
############################################

cont <- breedcont(Pedig, thisBreed='Hinterwaelder')
contByYear <- conttac(cont, cohort=Pedig$Born, use=Pedig$Breed=="Hinterwaelder", mincont=0.04)
round(contByYear,2)

barplot(contByYear, ylim=c(0,1), col=1:10, ylab="genetic contribution",
        legend=TRUE, args.legend=list(x="topleft",cex=0.6))
abline(1,0)
}




\author{
Robin Wellmann

Maintainer: Robin Wellmann <r.wellmann@uni-hohenheim.de>
}
\references{
Wellmann, R., and Pfeiffer, I. (2009). Pedigree analysis for conservation of genetic diversity and purging. {Genet. Res.} { 91}: 209-219

Wellmann, R., and Bennewitz, J. (2011). Identification and characterization of hierarchical structures in dog breeding schemes, a novel method applied to the Norfolk Terrier. {Journal of Animal Science}. 89: 3846-3858

Wellmann, R., Hartwig, S., Bennewitz, J. (2012). Optimum contribution selection for conserved populations with historic migration; with application to rare cattle breeds. {Genetics Selection Evolution}. 44: 34


}

\keyword{ package }

