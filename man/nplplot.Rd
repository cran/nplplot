\name{nplplot ver 3.0 rev 10}
\alias{nplplot}
\title{Non-parametric LOD score plotting}
\description{
Plots linkage score curves contained within one or more
specified linkage analysis results files. Usually these
results would be LOD scores, p-values, or log10(p-values).

A linkage results file has to be in a specific format with the
following three components:

A) First line = header line
B) Next set of lines = any number of data lines
C) Final two lines = line type & point type definition lines.

Here is an example:

\preformatted{
marker   location        score1  score2 score3      
d1s228   0.00            0.546   0.345  0.142                    
d1s429   1.00            0.346   0.335  0.252       
d1s347   2.00            0.446   0.245  0.342                    
ltype    -99.99              1       2      3       
ptype    -99.99             15      16     17                  
}

In this example, component A (the header line) is line1, component B
consists of lines 2-4; and component C (defining the line types & point
types) are the final two lines.

There can be any number of score columns within a file. Each column is
separated by white-space. The score columns within each file are
plotted inside the same plot. Each file is plotted as a separate
plot. Missing data in the \bold{score} columns can be denoted with either
"." or "NA". Missing labels in the \bold{marker} column can be denoted by a
"-".

The last two lines in this example refer to line type and point type
respectively. Both lines are necessary for nplplot to work correctly.
Setting "ltype" or "ptype" values to 0 prevents plotting curves and
points respectively.

Alternate plotting symbols can be selected by changing ptype and ltype
values in the score columns. For allowable ptype values, consult the R
documentation for "points". For line types, consult the documentation
on "par". }
\usage{
  nplplot(files, col=2, row=2, mode="p", output="screen", 
          yline=2.0, ymin=NULL, ymax=NULL, yfix=FALSE, batch=FALSE, 
          headerfiles=NULL, titles=NULL, xlabl="", ylabl="", 
          lgnd="page", lgndx=NULL, lgndy=NULL, bw=TRUE, na.rm=TRUE)
}
\arguments{
  \item{files}{List of files names (strings). Each file produces a
    separate plot.}
  
  \item{col}{For multiple plots on a page of postscript output,
    number of columns of plots, integer >= 1, defaults to 2.}
  \item{row}{For multiple plots on a page of postscript output,
    number of rows of plots, integer >= 1, defaults to 2. }
  \item{mode}{Orientation for postscript output,
    "p" for portrait
    "l" for landscape    }
  
  \item{output}{File name for saving plots; "screen", the default
    causes the plots to  be displayed on the screen.}

  \item{yline}{Y-value for displaying a horizontal cut-off line. }

  \item{ymin, ymax}{Y-axis minimum and maximum values with
    default values NULL. If non-NULL values are provided, then the
    plot area will be cropped to fit these values. Setting either of
    these values to NULL will enable plots to be expanded to fit the
    minimum and maximum y-values in the data files respectively.
  }
  \item{yfix}{Set to TRUE or FALSE depending on whether ymin and ymax
    should be enforced across all plots irrespective of whether the plot
    data lie within these bounds. Ignored if ymin or ymax are set to
    NULL.}
  
  \item{batch}{TRUE or FALSE, to determine whether the display
    screen should be closed. If nplplot is called within R, this should 
    be set to FALSE.}

  \item{headerfiles}{ List of file names, one for each data file
    specified above. Each header-file contains a string with column
    names corresponding to the columns in the data file. These column
    names are used in the plot legend. 
    If set to NULL (the default), nplplot tries to use the first item in 
    each column of a data file as plot legend.
    If a headerfile is provided, then nplplot will attempt to read in
    the first line of the datafile as data, so the user should be
    careful not to put in a headerline as well as a headerfile.}
    
  \item{titles}{Array of strings denoting titles for each plot.
    If there are not enough titles, the last string is recycled for the
    remaining plots. Default is an empty string.}

  \item{xlabl}{Array of strings, to use as the x-axis label on 
    each plot.}

  \item{ylabl}{Array of strings to use as the y-axis label on 
    each plot.}

  \item{lgnd}{TRUE, FALSE, "page" or a list of plot numbers denoting whether
    the legend should be drawn in all plots, none, first plot on a page,
    or specific plot numbers. Default "page"}

  \item{lgndx}{NULL or a real value if a specific x-coordinate should
    be used to position the legend. Default NULL.}

  \item{lgndy}{NULL or a real value if a specific y-coordinate should
    be used to position the legend. Default NULL.}

  \item{bw}{TRUE or FALSE depending on whether plots should be
    drawn in color. A list of colors are defined within nplplot. Then ltype
    values are used as indices into the list of color. These are
    1-magenta, 2-lightblue, 3-grey, 4=navyblue, 5-lightcyan, 6-pink,
    7-black. }
  \item{na.rm} {TRUE or FALSE depending on whether NAs should be removed
    prior to plotting the data. Including NAs will produce broken plots,
    when lines are drawn. This may be desirable in some cases, if
    missing data needs to be reported.
  }
}
\value{TRUE or FALSE depending on whether the input files were read 
	in successfully. }
\details{
  None.
}
\examples{
\dontrun{nplplot("lod.1", output="lod.1.ps", batch=T, headerfiles="hdr.1")}
\dontrun{nplplot(c("lod.1, lod.2"), col=1, row=2, headerfiles=c("hdr.1,hdr.2"))}
}
\keyword{hplot}


%  LocalWords:  ymin ymax nplplot headerfiles headerfile headerline
