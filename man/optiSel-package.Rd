\name{optiSel-package}
\alias{optiSel-package}
\alias{optiSel}
\docType{package}
\title{
\packageTitle{optiSel}
}
\description{
\packageDescription{optiSel}
}
\details{

\bold{Optimum Contribution Selection}

The following functions can be used to compute optimum contributions:

\tabular{ll}{
\link{opticont}      \tab  Calculates optimum genetic contributions of selection candidates to the next generation,\cr 
\link{opticont4mb}   \tab  Calculates optimum genetic contributions of selection candidates to  the next generation\cr
                     \tab   using multi-breed genotype data, \cr 
\link{opticomp}      \tab  Calculates breed composition of a multi-breed population with maximum diversity,\cr 
}

Function \link{noffspring} can then be used to  compute the optimum numbers of offspring of selection candidates from their optimum contributions. The available constraint settings and objective functions for a particular data set can be determined with the following help functions:

\tabular{ll}{
\link{help.opticont}   \tab Displays available objective functions and constraints for function \link{opticont},\cr 
\link{help.opticont4mb}   \tab  Displays available objective functions and constraints for function  \link{opticont4mb}.\cr 
}

For checking the results use

\tabular{ll}{
\link{summary.opticont}   \tab Calculates genetic parameters of the offspring population, and \cr
                          \tab checks for validity of optimum genetic contributions.\cr 
}

\bold{Kinships}

For pairs of individuals the following kinships can be computed:

\tabular{ll}{
\link{pedIBD}    \tab  Calculates \bold{ped}igree based probability of alleles to be \bold{IBD} ("pedigree based kinship""),\cr
\link{segIBD}    \tab  Calculates \bold{seg}ment based probability of alleles to  be \bold{IBD}  ("segment based kinship"),\cr
\link{pedIBDatN} \tab  Calculates \bold{ped}igree based probability of alleles to be \bold{IBD} \bold{at} segments with \bold{N}ative origin,\cr
\link{segIBDatN} \tab  Calculates \bold{seg}ment based probability of alleles to be \bold{IBD}  \bold{at} segments with \bold{N}ative origin,\cr
\link{pedIBDorM} \tab  Calculates \bold{ped}igree based probability of alleles to be \bold{IBD} \bold{or} \bold{M}igrant alleles,\cr
\link{segIBDandN} \tab Calculates \bold{seg}ment based probability of alleles to be \bold{IBD}  \bold{and} have \bold{N}ative origin,\cr 
\link{segN}       \tab Calculates \bold{seg}ment based probability of alleles to have \bold{N}ative origin.\cr 
}

Results from these functions can be combined with function \link{kinlist} into a single R object, which can then be used as an argument
to function \link{opticont} or \link{opticont4mb}.  Function \link{sim2dis} can be used to convert a similarity matrix (e.g. a kinship matrix) into a dissimilarity matrix which is suitable for multidimensional scaling. 

\bold{Breed Composition}

The breed composition of crossbred individuals can be accessed with

\tabular{ll}{
\link{pedBreedComp} \tab Calculates \bold{ped}igree based the \bold{Breed} \bold{Comp}osition, which is the genetic contribution\cr 
                    \tab of each individual from other breeds and from native founders. The native contribution \cr 
                    \tab is the proportion of the genome not originating from other breeds.\cr 
\link{segBreedComp} \tab  Calculates \bold{seg}ment based the \bold{Breed} \bold{Comp}osition. The native contribution is the \cr 
                    \tab proportion of the genome belonging to segments that have low frequency in \cr
                    \tab other breeds.\cr
}

The native contributions obtained by the above functions can be contrained or maximized with functions \link{opticont} or \link{opticont4mb} to remove introgressed genetic material, or alternatively, they can be considered a quantitative trait and included in a selection index.

\bold{Haplotype frequencies} 

Frequencies of haplotype segments in particular breeds can be computed and plotted with 

\tabular{ll}{
\link{haplofreq} \tab  Calculates the maximum frequency each segment has in a set of reference breeds,\cr
                 \tab and the name of the breed in which the segment has maximum frequency. \cr
                 \tab Identification of native segments.\cr 
\link{freqlist} \tab Combines results obtained with function \link{haplofreq} for different reference breeds \cr
                \tab into a single R object which is suitable for plotting.\cr
\link{plot.HaploFreq} \tab Plots frequencies of haplotype segments in particular reference breeds.\cr
}

\bold{Inbreeding Coefficients and Genetic Contributions} 

The inbreeding coefficients and genetic contributions from ancestors can be computed with:

\tabular{ll}{
\link{pedInbreeding}  \tab  Calculates  \bold{ped}igree based  \bold{Inbreeding}.\cr 
\link{segInbreeding}  \tab  Calculates  \bold{seg}ment based  \bold{Inbreeding}, i.e. inbreeding based on \cr 
                      \tab  runs of homozygosity (ROH).\cr 
\link{genecont}   \tab  Calculates  \bold{gene}tic  \bold{cont}ributions  each individual has from all it's ancestors in \cr
\tab the pedigree.\cr 
}


\bold{Preparing and plotting pedigree data} 

There are some functions for preparing and plotting pedigree data

\tabular{ll}{
\link{prePed}      \tab \bold{pre}pares a \bold{Ped}igree by sorting, adding founders and pruning the pedigree,\cr 
\link{completeness}\tab Calculates pedigree \bold{completeness} in all ancestral generations,\cr
\link{summary.Pedig}\tab Calculates number of equivalent complete generations, number of fully \cr 
                    \tab traced  generations, number of maximum generations traced, index of \cr 
                    \tab pedigree completeness, inbreeding coefficients,\cr
\link{subPed}       \tab Creates a \bold{sub}set of a large \bold{Ped}igree,\cr 
\link{pedplot}     \tab Plots a pedigree,\cr 
}


\bold{Population Parameters} 


Finally, there are some functions for estimating population parameters:

\tabular{ll}{
\link{conttac}   \tab Calculates genetic \bold{cont}ributions of breeds \bold{t}o  \bold{a}ge \bold{c}ohorts,\cr
\link{kinwac}         \tab Calculates for every individual it's mean \bold{kin}ship \bold{w}ith the \bold{a}ge \bold{c}ohort to which \cr
    \tab it belongs,\cr
\link{summary.kinwac} \tab Calculates for every age cohort several genetic parameters. These may \cr
                      \tab include average kinships, genetic diversities, gene diversity at native loci,\cr
                      \tab the native effective size, and the native genome equivalent.\cr
}



\bold{Genotype File Format}

All functions reading phased genotype data assume that the files are in the following format:

Each file has a header and no row names. Cells are separated by blank spaces. The number of rows is equal to the number of markers from the respective chromosome and the markers are in the same order as in the \code{map}. There can be some extra columns on the left hand side containing no genotype data. The remaining columns contain genotypes of individuals written as two alleles separated by a character, e.g. A/B, 0/1, A|B, A B, or 0 1. The same two symbols must be used for all markers. Column names are the IDs of the individuals. If the blank space is used as separator then the ID of each individual should be repeated in the header to get a regular delimited file. The columns to be skipped and the individual IDs must have no white spaces.

Use function \link{read.indiv} to extract the IDs of the individuals from a genotype file. 





}
\author{
\packageAuthor{optiSel}

Maintainer: \packageMaintainer{optiSel}
}
\references{
de Cara MAR, Villanueva B, Toro MA, Fernandez J (2013). Using genomic tools to maintain diversity and fitness in conservation programmes. Molecular Ecology. {22}: 6091-6099

Wellmann, R., and Pfeiffer, I. (2009). Pedigree analysis for conservation of genetic diversity and purging. {Genet. Res.} { 91}: 209-219

Wellmann, R., and Bennewitz, J. (2011). Identification and characterization of hierarchical structures in dog breeding schemes, a novel method applied to the Norfolk Terrier. {Journal of Animal Science}. 89: 3846-3858

Wellmann, R., Hartwig, S., Bennewitz, J. (2012). Optimum contribution selection for conserved populations with historic migration; with application to rare cattle breeds. {Genetics Selection Evolution}. 44: 34

Wellmann, R., Bennewitz, J., Meuwissen, T.H.E. (2014) A unified approach to characterize and conserve adaptive and neutral genetic diversity in subdivided populations. {Genetics Selection Evolution}. 69, e16
}

\keyword{ package }

\examples{

#See ?opticont for optimum contribution selection 
#See ?opticont4mb for multi-breed optimum contribution selection 
#These examples demonstrate computation of some population genetic parameters.

data(ExamplePed)
Pedig <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
head(Pedig)


############################################
# Evaluation of                            #
#    - kinships                            #
#    - genetic diversities                 #
#    - native effective size               #
#    - native genome equivalent            #
############################################

Kinships <- kinlist(pedIBD=pedIBD(Pedig), pedIBDatN=pedIBDatN(Pedig, thisBreed="Hinterwaelder"))
Kin      <- kinwac(Kinships, Pedig=Pedig, use=Pedig$Breed=="Hinterwaelder")
Kin[1000:1010,]

x  <- summary(Kin, tlim=c(1970,2005), histNe=150, base=1800)
xP <- x$Param

plot(xP$cohort, xP$Ne, type="l", ylim=c(0,150), 
     main="Native Effective Size",ylab="Ne",xlab="")

matplot(xP$cohort, xP[,c("pedIBD", "pedIBDatN")], 
        type="l",ylim=c(0,1),main="Kinships", xlab="Year", ylab="mean Kinship")
abline(0,0)
legend("topleft", legend = c("pedIBD", "pedIBDatN"),lty=1:2,col=1:2,cex=0.6)


plot(xP$cohort, xP$condGD, type="l", ylim=c(0,1), 
     main="Diversity at native alleles",ylab="condGD",xlab="")
mtext(paste("Base Year =", x$base, "  historic Ne =", x$histNe), cex=0.7)

plot(xP$cohort,xP$NGE,type="l",main="Native Genome Equivalents", 
     ylab="NGE",xlab="",ylim=c(0,7))
mtext(paste("Base Year =", x$base, "  historic Ne =", x$histNe), cex=0.7)


############################################
# Genetic contributions from other breeds  #
############################################

cont <- pedBreedComp(Pedig, thisBreed='Hinterwaelder')
contByYear <- conttac(cont, Pedig$Born, use=Pedig$Breed=="Hinterwaelder", mincont=0.04, long=FALSE)
round(contByYear,2)

barplot(contByYear, ylim=c(0,1), col=1:10, ylab="genetic contribution",
        legend=TRUE, args.legend=list(x="topleft",cex=0.6))


######################################################
# Frequencies of haplotype segments in other breeds  #
######################################################

data(map)
data(Cattle)
dir   <- system.file("extdata", package="optiSel")
files <- file.path(dir, paste("Chr", 1:2, ".phased", sep=""))

Freq <- freqlist(
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Rotbunt",   minSNP=20),
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Holstein",  minSNP=20),
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Fleckvieh", minSNP=20)
  )

plot(Freq, ID=1, hap=2, refBreed="Rotbunt")


}
