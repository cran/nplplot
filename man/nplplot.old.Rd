\name{nplplot.old}
\alias{nplplot.old}
\title{LOD score plotting (old version of nplplot)}
\description{Plots score curves contained within one or more
specified results files.}
\usage{nplplot.old(files, col=2, row=2, mode="p", output="screen", 
          yline=2.0, ymin=NULL, ymax=NULL, yfix=FALSE, batch=FALSE, 
          headerfiles=NULL, titles=NULL, xlabl="", ylabl="", 
          lgnd="page", lgndx=NULL, lgndy=NULL, bw=TRUE, na.rm=TRUE)}
\arguments{
  \item{files}{ List of files names (strings). Each file produces a
    separate plot.}
  
  \item{col}{ For multiple plots on a single page of pdf or postscript
    output, this item defines the number of columns of plots, and should
    be an integer greater than or equal to 1. Default is set to 2.}
  
  \item{row}{ For multiple plots on a page of pdf or postscript output,
    this defines the number of rows of plots, (value should be 1 or
    greater). Default value is set to 2. }
  
  \item{mode}{ Orientation for pdf or postscript output,
    "p" for portrait
    "l" for landscape. }
  
  \item{output}{ File name for saving plots; "screen", the default
    causes the plots to  be displayed on the screen. To produce a pdf
    file use the extension .pdf. To produce a postscript file, use the
    .ps file name extension. If no extension is given a pdf file is
    produced. }

  \item{yline}{ Y-value for displaying a horizontal cut-off line. }

  \item{ymin, ymax}{ Y-axis minimum and maximum values with
    default values NULL. If non-NULL values are provided, and
    \code{yfix} is set to TRUE, then the plot area will be cropped to
    these values.}
  
  \item{yfix}{ Set to TRUE or FALSE depending on whether ymin and ymax
    should be enforced across all plots irrespective of whether the plot
    data lie within these bounds. Ignored if ymin or ymax are set to
    NULL.}
  
  \item{batch}{ TRUE or FALSE, to determine whether the display
    screen should be closed. If nplplot is called within R, this should 
    be set to FALSE.}

  \item{headerfiles}{ List of file names, one for each data file
    specified above. Each header-file contains a string with column
    names corresponding to the columns in the data file. These column
    names are used in the plot legend. 
    If set to NULL (the default), nplplot uses the first item in 
    each column of a data file as plot legend.
    If a headerfile is provided, then nplplot will attempt to read in
    the first line of the datafile as data, so the user should be
    careful not to put in a headerline as well as a headerfile.}
  
  \item{titles}{ Array of strings denoting titles for each plot.
    If there are not enough titles, the last string is recycled for the
    remaining plots. Default is an empty string. }

  \item{xlabl}{ Array of strings, to use as the x-axis label on 
    each plot.}

  \item{ylabl}{ Array of strings to use as the y-axis label on 
    each plot. }

  \item{lgnd}{ TRUE, FALSE, "page" or a list of plot numbers denoting whether
    the legend should be drawn in all plots, none, first plot on a page,
    or specific plot numbers. Default "page". }

  \item{lgndx}{ NULL or a real value if a specific x-coordinate should
    be used to position the legend. Default NULL. }

  \item{lgndy}{ NULL or a real value if a specific y-coordinate should
    be used to position the legend. Default NULL. }

  \item{bw}{ TRUE or FALSE depending on whether plots should be
    drawn in color. A list of six colors are defined within nplplot,
    which are successively used to draw each curve, and reused as
    necessary. The order in which these colors are used is:
    magenta, lightblue, grey, navyblue, lightcyan and pink.
    The 7th color, reserved for black and white plots is black. }
    
  \item{na.rm}{ TRUE or FALSE depending on whether NAs should be removed
    prior to plotting the data. Including NAs will produce broken plots,
    when lines are drawn. This may be desirable in some cases, if
    missing data needs to be reported. }

}

\details{
Usually these results would be 
LOD scores, p-values, or log10(p-values). This is targetted
towards p-values or LOD scores obtained at various marker positions
from statistical analysis of genetic data.
A results file has to be in a specific tabular format with
each column separated by white-space :

A) First line = header line                                      

B) Next set of lines = any number of data lines

C) Final two lines = line type & point type definition.

Here is an example:
\preformatted{
marker   location        score1  score2 score3      
d1s228   0.00            0.546   0.345  0.142                    
d1s429   1.00            0.346   0.335  0.252       
d1s347   2.00            0.446   0.245  0.342                    
ltype    -99.99              1       2      3       
ptype    -99.99             15      16     17                  
}

In this example, line 1 column headers for the score columns
may be used as labels within the legend, as described in the 
usage of the "headerfile" argument. \bold{The first two headers 
are ignored.}

Lines 2-4 contain scores at various marker positions. \bold{Missing
scores can be denoted with either "." or "NA"}. The position column
cannot have missing data. There can be any number of score columns
within a file and will be plotted as separate curves within the same
plot.  Each file is plotted as a separate plot.

The last two lines give line types and point types for each curve. A
zero line or point type will not plot lines or points for that score
column respectively. For allowable ptype values, consult the R
documentation for "points". For line types, consult the documentation
on "par".

The names in the first column are used as axis labels on the top of
the plot border. Setting a name in the \bold{marker} column to "-"
will result in no label at that position.
}

\value{TRUE or FALSE depending on whether the input files were read 
  in successfully.}

\seealso{\code{\link{nplplot}, \link{nplplot.multi}}}
\examples{
\dontrun{nplplot.old("lod.1", output="lod.1.ps", batch=T, headerfiles="hdr.1")}
\dontrun{nplplot.old(c("lod.1", "lod.2"), col=1, row=2, headerfiles=c("hdr.1","hdr.2"))}
}
\keyword{plotting}

%  LocalWords:  ymin ymax nplplot headerfiles headerfile headerline
