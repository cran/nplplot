\name{nplplot}
\alias{nplplot}
\title{Plotting statistics along a chromosome}
\description{Plots linkage or association statistics along a chromosome,
  contained within a data frame or a file. Marker names are displayed
  along the top border.}
\usage{nplplot(plotdata=NULL, filename=NULL, yline=2.0, ymin=0, ymax=3.0,
         header=TRUE, yfix=FALSE, title=NULL, draw.lgnd=TRUE,
         xlabl="", ylabl="", lgndx=NULL, lgndy=NULL, lgndtxt=NULL,
         cex.legend = 0.7, cex.axis=0.7, tcl=1,
         bw=TRUE, my.colors=NULL, ltypes=NULL, ptypes=NULL,
         na.rm=TRUE, plot.width=0.0, ...) }

\arguments{
  \item{plotdata}{ A data frame containing marker names in the
  first column, marker map positions in the second column, and
  statistical scores in column 3 onwards. }

  \item{filename}{ A table format file containing the plot
  data as described above. }

  \item{header}{ TRUE or FALSE depending on whether the plotdata
  or file has a header line. }

  \item{yline}{ Y-value for displaying a horizontal cut-off
  line. If 'yfix' is set to TRUE and Y-line falls outside of \code{[ymin,
  ymax]}, then the cut-off line is omitted. }

  \item{ymin, ymax}{ Y-axis minimum and maximum values. If non-NULL
    values are provided, and \code{yfix} is set to TRUE, then the plot
    area will be cropped to these values. If \code{yfix} is set to
    FALSE, then \code{ymin} and \code{ymax} values are ignored. }

  \item{yfix}{ TRUE or FALSE to denote whether plot area should be
  cropped to the \code{ymin, ymax} values. This has no effect if
  \code{ymin, ymax} values are NULL. }

  \item{title}{ Used as the subtitle of the plot.}

  \item{xlabl}{ X-axis label. May interfere with the display of
    the subtitle provided as the \code{title} argument. }

  \item{ylabl}{ Y-axis label.}

  \item{draw.lgnd}{ TRUE or FALSE denoting whether a plot legend
  should be displayed. }

  \item{lgndx}{ X coordinate for the legend box, passed to the
  \code{\link{legend}} command. Ignored if \code{draw.legend}
  is set to FALSE. If set to NULL with \code{draw.legend} set to TRUE,
  the X-coordinate is automatically calculated. }

  \item{lgndy}{ Y coordinate for the legend box, passed to the
  \code{\link{legend}} command. Ignored if \code{draw.legend}
  is set to FALSE. If set to NULL with \code{draw.legend} set to TRUE,
  the Y-coordinate is automatically calculated. }

  \item{lgndtxt}{ Vector of strings to use in the legend. }

  \item{cex.legend}{ Character scaling for legend, passed as the
  \code{cex} argument to the  \code{\link{legend}} command. }

  \item{cex.axis}{ Character scaling for the axis, passed to the
  \code{\link{axis}} command for drawing the top border. }

  \item{tcl}{ Length of ticks for the top border, passed to the
  \code{\link{axis}} command. }

  \item{bw}{ TRUE or FALSE depending on whether plots should be
  drawn in color. If set to FALSE, then the colors defined by
  \code{my.colors} are used. }

  \item{my.colors}{ Vector of color specifications as described
  in the \code{\link{par}} command. Ignored if \code{bw} above
  is set to FALSE. If \code{bw} is to TRUE and my.colors is set to
  NULL, the \code{\link{rainbow}} palette will be used. }

  \item{ltypes}{ Vector of line types for the plots. Each
  non-zero line type is passed on to a \code{lines} command. Use 0 or
  'none' if a line is to be skipped. If NULL, no lines will be
  drawn. For line types see the \code{\link{par}} command. If set to
  "default", line-types 1 through the number of plots is used. }

  \item{ptypes}{ Vector of characters giving the point
  types, to be passed onto the \code{points} command. Use 'none' if no
  points are to be drawn for a score column. If NULL, no points will
  be displayed. If both the line-type and point-type specification for
  a results column is set to 'none', that column will not be plotted. }

  \item{na.rm}{ TRUE or FALSE depending on whether points with
  Y-coordinates set to NAs should be skipped. Setting \code{na.rm} to
  TRUE eliminates discontinuities in the plots. }

  \item{plot.width}{ A number giving the width of the plot in
  inches. This is used to decide whether some marker names should be
  dynamically hidden, if they are too close to each other along the top
  border. If set to 0, the default page-size is used to set the width. }

  \item{...}{ Further graphical parameters to be passed onto the
  'plot', 'lines' and 'points' commands. }
}

\details{ The nplplot function draws multiple curves within a single
plot by automatically calling 'plot', 'lines', and 'points' multiple
times, thus making it easy for the user to plot many columns of results
using a single plot command. It is intended for the display of linkage
and association analysis results such as LOD scores and P-values. It
allows the marker names to be displayed along the top border of the
plot, as well as a significance threshold line.

The input plot data has to be in a specific tabular format with
each column separated by white-space :

Here is an example:
\preformatted{
Marker   Position        score1  score2 score3      
d1s228   0.00            0.546   0.345  0.142                    
d1s429   1.00            0.346   0.335  0.252       
d1s347   2.00            0.446   0.245  0.342                    
}
This example file contains a header, therefore the \code{header} argument should
be set to TRUE.

Lines 2-4 contain scores at various marker positions. \bold{Missing
scores can be denoted with either "." or "NA"}. The position column
cannot have missing data. There can be any number of score columns
within a file and will be plotted as separate curves within the same
plot.  Each file is plotted as a separate plot.
}

\value{TRUE or FALSE depending on whether the plot data was successfully
  plotted.}

\seealso{\code{\link{nplplot.multi}, \link{nplplot.old}}}

\examples{
# plot with legend
par(omi=c(0.05, 0.05, 0.5, 0.05))
data(lods1, package="nplplot")
nplplot(plotdata=lods1, draw.lgnd=TRUE)

# plot without legend
data(lods2, package="nplplot")
nplplot(plotdata=lods2, draw.lgnd=FALSE)

# plotting from a data file
datadir <- paste(system.file("data", package="nplplot"), .Platform$file.sep, sep="")
nplplot(filename=paste(datadir, "lods2.txt.gz", sep=""))

}
\keyword{plotting}


%  LocalWords:  lgndtxt
