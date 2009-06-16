% --- Source file: bedplot.Rd ---
\name{bedplot}
\alias{bedplot}
\title{Creation of BED and BedGraph custom tracks}
\description{Generates matched sets of files for linkage or association
  statistics along a chromosome for viewing in the UCSC genome browser
  from an input file containing a table of marker names, physical
  positions and one or more statistical scores.}  
\usage{
 bedplot(bed.data)
}
\arguments{
  \item{bed.data}{File containing a table of marker names, physical position and
    scores. }
}
\details{
bed.data example:

\preformatted{
Marker   Position   TRAIT_ALL   
M1       144255     0.670
-        144305     0.640
M3       144355     0.590
-        144378     0.600
M2       144400     0.610
}

Bedplot creates two types of files: a BED.\# file containing a custom
BED annotation track and a BedGraph.\# file custom BedGraph annotation
track. These files have the same suffix as the input bed.data file. When
there are multiple scores in \code{bed.data} file, a matched pair of
BedGraph track BED track files is created for each score, labelled with
the score names, as well as the chromosome numbers,
e.g. \code{BedGraph.score1.#} and \code{BED.score1.#},
\code{BedGraph.score2.#} and \code{BED.score2.#}, etc.  
}
\value{TRUE or FALSE depending on whether runs successfully. }
\examples{
\dontrun{bedplot("bed.data.05")}
}
\keyword{utilities}


