% --- Source file: genomeplot.Rd ---
\name{genomeplot}
\alias{genomeplot}
\title{Creation of Genome Graph files}
\description{
The \code{genomeplot} function generates two formatted files, one
containing \dQuote{chromosome base} formatted genome data and the other
containing marker-specific results with dbSNP SNP IDs for displaying  
genome-wide data sets in the UCSC genome browser. }
\usage{
genomeplot(gg.data)
}
\arguments{
  \item{gg.data}{a file containing chromosome, marker, physical position and scores. }
}
\details{
gg.data example:

\preformatted{
Chromosome   Marker   Position   TRAIT_ALL   
         5   M1       0.000      0.670
         5   -        2.500      0.640
         5   M3       5.000      0.590
         5   -        6.500      0.600
         5   M2       8.000      0.610                   
         8   M4       0.000      0.670
         8   -        2.500      0.640
         8   M6       5.000      0.590
         8   -        6.500      0.600
         8   M5       8.000      0.610 
}

Two files are created, \dQuote{GG.positons.all} for the
\dQuote{chromosome base} format, and \dQuote{GG.markers.all} for the
marker-names based format. When there are multiple scores in
\code{gg.data} file, this results in matched pairs of files, one for
each score, labelled with the score names, e.g. GG.positions.score1.all,
and GG.markers.score1.all, GG.positions.score2.all and
GG.markers.score2.all, and so on.

}
\value{TRUE or FALSE depending on whether runs successfully. }
\examples{
\dontrun{genomeplot("GG.data.all")}
}
\keyword{utilities}


