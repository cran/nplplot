\name{nplplot.multi}
\alias{nplplot.multi}
\title{Plotting linkage or association statistics for multiple results files}
\description{Wrapper function for the 'nplplot' function. Creates mutiple plots
from a list of plot files, with custom graphical parameters set by
header files. }
\usage{nplplot.multi(filenames, col=2, row=2, mode="l", output="screen", 
                headerfiles=NULL, lgnd="page", customtracks=FALSE,
                mega2mapfile=NULL, pagewidth=NULL, pageheight=NULL,
                topmargin=0.25, ...)
}
\arguments{
  \item{filenames}{ Vector of strings giving file names
  containing tables of linkage analysis results. See
  \code{\link{nplplot}} for a description of the file format. }

  \item{col}{ Integer indicating number of columns of plots to be drawn
  on a page.}

  \item{row}{ Integer signifying number of rows of plots will be drawn
  on a page.}

  \item{mode}{ 'p' or 'l' to denote 'portrait' or 'landscape' mode.}

  \item{output}{ String giving file name to save plots in. If
  set to 'screen', plots will be displayed and not saved. The file
  format id decided by the filename extension, '.pdf' for PDF, or
  '.ps' for postscript. If no extension is provided, or is not
  recognized, a PDF file will be produced with '.pdf' appended to the
  file name.}

  \item{headerfiles}{ Files containing R language commands to set
  various plot parameters, which are passed onto the
  \code{\link{nplplot}} command. The recommended use is to have one
  headerfile per plot file. For a list of parameters, consult the
  \code{\link{nplplot}} documentation. If the number of headerfiles is
  fewer than plot files, the last header file will be reused as many
  times as needed. If more headerfiles are provided than necessary, the
  last ones will be ignored.}

  \item{lgnd}{ TRUE, FALSE, 'page' or a list consisting of plot
    numbers. If a single value is given, TRUE causes legends to be drawn
    inside every plot, FALSE omits legends altogether, and 'page' causes
    a legend to be drawn inside the first plot on every page. If a list
    of numbers is provided, only plots corresponding to these numbers
    will have legends. }

  \item{customtracks}{ TRUE or FALSE. If set to TRUE, data files are
    created to draw custom tracks within the UCSC genome browser in BED
    format, as well as a combined data file to add a genome-wide track
    over all chromosomes present in the data. If set to TRUE, a
    \code{mega2mapfile} also needs to be supplied (see below). }

  \item{mega2mapfile}{ Mega2 annotated format map file containing
   physical positions for all SNPs and chromosomes present in the
   nplplot input data files. }

  \item{pagewidth}{ A number denoting width of the plot page in
  inches. If set to NULL, a width of 7.0 is used for the plot
  area. Assumes that a margin of 0.5 will be available around the plot
  area for axis annotations. }

  \item{pageheight}{ A number denoting height of the plot page in
  inches. If set to NULL, a height of 10.0 is used for the plot area. 
  Assumes that a margin of 0.5 will be available around the plot
  area for axis annotations. }

  \item{topmargin}{ A number denoting the width of the outside top
  margin of each plot. Since this contains marker names, it may need to
  be increased to accommodate long names.}

  \item{...}{ Further graphical parameters to be passed onto the
  'plot', 'lines' and 'points' commands within \code{nplplot}. }

}

\details{ This function is designed for use within the Mega2 software to
generate graphical output for some of the target analysis options,
namely Merlin, SimWalk2 and Allegro. It calls nplplot repeatedly to
create plots corresponding to each input file. The input arguments
control characteristics of all plots together, whereas the header
files allow customization within each plot. Thus, it is expected that
there should be as many header files as there are plot data files.

This function can also be used to create custom tracks within the UCSC
genome browser, as well as a genome-wide plot. To use this feature, make
sure that the names of the nplplot input data files each have a
"Mega2-style" chromosome extension (01 through 09 for chromosomes 1
thorugh 9, otherwise the chromosome number, or X for the human
sex-chromosome, 23). }

\value{TRUE or FALSE depending on whether all plot commands were successful.}

\seealso{\code{\link{nplplot}, \link{nplplot.old}}}
\examples{
datadir <- paste(system.file("data", package="nplplot"),
                 .Platform$file.sep, sep="")
f1 <- paste(datadir, "lods1.txt", sep="")
f2 <- paste(datadir, "lods2.txt", sep="")
h1 <- paste(datadir, "lods1header.R", sep="")
h2 <- paste(datadir, "lods2header.R", sep="")
nplplot.multi(c(f1, f2), col=1, row=2, output="screen", 
	     headerfiles=c(h1, h2), topmargin=0.5)
}

\keyword{hplot}
