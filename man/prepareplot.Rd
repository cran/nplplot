% --- Source file: bedplot.Rd ---
\name{prepareplot}
\alias{prepareplot}
\title{Prepare input data files for bedplot and genomeplot}
\description{The prepareplot function prepares input data files for
  \link{bedplot} and \link{genomeplot} functions from nplplot-formatted
  score files and a Mega2 annotated format map file with physical positions. }  
\usage{
 prepareplot(prefix, chrlist=1:24, mapfile, output="both")
}
\arguments{
  \item{prefix}{Prefix of the names of R table files, e.g. \dQuote{RMERLINDATA} for R table files 
  \dQuote{RMERLINDATA.01}, \dQuote{RMERLINDATA.02}, etc. Using
  \code{chrlist} below, it automatically finds R table files with the
  specified prefix and chromosome-specific extensions to convert. }
  
  \item{chrlist}{List of chromosome numbers to create plots for, default 1 through
  24. Chromosomes 23 and 24 produces files for the X chromosome X, 24
  denoting pseudo-autosomal markers on chromosome X. }
  
  \item{mapfile}{Mega2 annnotated format map file, containing
    marker names and and exactly one set of physical positions. }
  
  \item{output}{Which plotting function to generate data for,
    \dQuote{both} for both bedplot and genomeplot functions,
    \dQuote{bed} for generating input files for bedplot function,
    \dQuote{GG} for generating input file for genomeplot
    function. output is set to default \dQuote{both}. }
}

\details{
    mapfile example:

\preformatted{
Chromosome   Map.h.a   Name   Map.h.m   Map.h.f   Build52.p
5            0.0       M1     0.0       0.0       144255
5            5.0       M3     2.0       7.0       144355
5            8.0       M2     4.0       12.0      144400
8            0.0       M4     0.0       0.0       144255
8            5.0       M6     2.0       7.0       144355
8            8.0       M5     4.0       12.0      144400
}

The names of R table files should be linkage or association analysis
score files in nplplot-format with Mega2-style file names, i.e., having
a common specified prefix and 01-09, 11- 24, X, or XY as suffixes. The
list of suffixes are determined by the chromosome list. If this list
includes 23 or X, R table files with either the \dQuote{23}
suffix or \dQuote{X} suffix are accepted. If both files exist, the one
with the \dQuote{X} suffix is read in and the user warned. If the XY
chromosome is chosen, R table files can have either \dQuote{24} or
\dQuote{XY} as a suffix, with \dQuote{XY} suffixed file having
precedence.

The prepareplot function generates chromosome-specific formatted score
files \dQuote{bed.data.\#} for use by \link{bedplot} with the same
suffix as the R table file. If X chromosome is chosen, the output file
is named \dQuote{bed.data.23}. If XY chromosome is chosen, those records
on XY chromosome are included in \dQuote{bed.data.23} file. The output
file \dQuote{bed.data.\#} contains marker names and physical
positions followed by one or more score columns. The header is taken
from the input score file(s).

Prepareplot generates a combined file over all chromosomes
\dQuote{GG.data.all} for \link{genomeplot}. For pseudo-autosomal markers
denoted by chromosome XY or 24, these scores are assigned the X
chromosome. The output file \dQuote{GG.data.all} contains four or more
columns with headings. The first, second and third columns contain
chromosomes, marker names and physical positions respectively, followed
by one or more score columns with score names as headers.

}
\value{TRUE or FALSE depending on whether runs successfully. }
\examples{
\dontrun{prepareplot("RMERLINDATA", c(5,8), "map.all", "GG")}
}
\keyword{utilities}


