\name{plot.HaploFreq}
\alias{plot.HaploFreq}

\title{Plots Frequencies of Haplotype Segments in Particular Breeds}

\description{
For a particular haplotype from \code{thisBreed} and each marker \code{m} the frequency of the segment containing marker \code{m} in a particular breed is plotted.
}
\usage{
\method{plot}{HaploFreq}(x, ID=1, hap=1, refBreed=NULL, Chr=NULL, ...)
}

\arguments{
\item{x}{This is either an R-Object obtained with function \link{haplofreq} or a list obtained with function \link{freqlist}. If \code{Freq} is a list then it has one component for each reference breed, whereby the name of the component is the name of the reference breed. Each component is a Mx(2N) matrix containing for every marker m and for each haplotype from \code{thisBreed} the frequency of the segment containing marker m in the respective breed.
}
\item{ID}{Either the ID of the animal from this breed to be plotted, or the position of the animal in the matrices.}
\item{hap}{Number of the haplotype to be plotted (\code{1} or \code{2})}
\item{refBreed}{Breed name. Frequencies the haplotype segments have in this reference breed will be plotted.
This parameter is only needed if \code{Freq} is a list containing the segment frequencies for several reference breeds. The default means that the frequencies will be plotted for the first breed in the list.}
\item{Chr}{Vector with chromosomes to be plotted. The default means that all chromosomes will be plotted.}
\item{...}{Arguments to be passed to methods, such as graphical parameters.}
}
\details{
 For a particular haplotype from \code{thisBreed} and each marker \code{m} from chromosomes \code{Chr} the frequency of the segment containing marker \code{m} in reference breed \code{refBreed} is plotted (red line), as well as the maximum frequency the segment has in one of the evaluated breeds (black line), and the maximum frequency a segment from \code{thisBreed} has in one of the evaluated breeds (grey area).
}
\examples{
data(map)
data(Cattle)
dir   <- system.file("extdata", package="optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")

Freq <- freqlist(
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Rotbunt",   minSNP=20),
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Holstein",  minSNP=20),
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Fleckvieh", minSNP=20)
  )

names(Freq)

plot(Freq, ID=1, hap=2, refBreed="Rotbunt")
plot(Freq, ID=1, hap=2, refBreed="Holstein", Chr=1)

Freq <- haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="others",   minSNP=20)
plot(Freq, ID=1, hap=2)

Freq <- haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Angler",   minSNP=20)
plot(Freq, ID=1, hap=2)

}
\author{Robin Wellmann}

