\name{plot_fullsib_scan}

\alias{plot_fullsib_scan}

\title{Plot the QTL mapping scan}

\description{
  Plot the QTL mapping scan resulted from the \code{im_scan} and
  \code{cim_scan}. This is fully based on \code{plot.scanone} from
  \pkg{qtl} package.
} 

\usage{
\method{plot}{fullsib_scan}(x, x2, x3, x4, x5, lg, label.lg, cex.axis=1,
incl.mkr = c("none", "points", "name"), cex.incl, xlim, ylim, lty = 1,
lwd = 2, add=FALSE, gap = 50, col= c("black", "blue", "red", "orange",
"limegreen"), alternate.lgid = FALSE, \ldots)
}

\arguments{
  
  \item{x}{An object of class \emph{fullsib_scan}, output of \code{im_scan}
    and \code{cim_scan}.}
  
  \item{x2}{Second object of class \emph{fullsib_scan}, optional.}

  \item{x3}{Third object of class \emph{fullsib_scan}, optional.}

  \item{x4}{Fourth object of class \emph{fullsib_scan}, optional.}

  \item{x5}{Fifth object of class \emph{fullsib_scan}, optional.}
  
  \item{lg}{Indicates which linkage groups will be plotted, default is
    to plot all groups.}

  \item{label.lg}{for more than one linkage group plotted, a label for
    each group can be easily added to the plot, the length must be the
    same of the number of LG.}

  \item{cex.axis}{The same option for argument used in par function for
    graphical parameters. See \sQuote{?par} for details.}
  
  \item{incl.mkr}{indicates how the markers will be displayed on the
    plot. It can assume three options: \code{points}, \code{name},
    \code{none}, in which the first indicates the position of the markers
    with small triangules, the second indicates by its name (sometimes
    it is useful if only one group is plotted, otherwise it can be
    difficult to distinguish the names) and the last none additional
    information is given. Here, the default is to consider
    \code{incl.mkr = points}.}

  \item{cex.incl}{If incl.mkr is indicated as \code{points} or
    \code{name}, an integer can be defined to resize the symbols (for
    points) or text (for names) displayed at the graphic.}

  \item{xlim}{Limits for x-axis (optional).}
  
  \item{ylim}{Limits for y-axis (optional).}

  \item{lty}{Line types; vector of length 1 until 5, according to the
    number of object of class \emph{fullsib_scan} plotted.}
  
  \item{lwd}{Line widths; vector of length 1 until 5, according to the
    number of object of class \emph{fullsib_scan} plotted.}

  \item{add}{If TRUE, add to a current plot}

  \item{gap}{Gap separating chromosomes (in cM).}
  
  \item{col}{Line colors; vector of length 1 until 5, according to the
    number of object of class \emph{fullsib_scan} plotted.}
  
  \item{alternate.lgid}{If TRUE and more than one linkage group is
     plotted, alternate the placement of chromosome axis labels, so that
     they may be more easily distinguished.}

  \item{\dots}{Passed to the function \code{plot} when it
    is called.}

}

\value{None}

\details{
  This function is basically \code{plot.scanone} written by Karl Broman
  and adapted for this package. We also adapt the functions to receive
  more object to plot QTL profile.
}

\seealso{
  \code{\link[fullsibQTL]{im_scan}},
  
  \code{\link[fullsibQTL]{cim_scan}},
  
  \code{\link[fullsibQTL]{summary_fullsib_scan}},
  
  \code{par},
  
  \code{scanone} from \pkg{qtl} package,

  \code{plot.scanone} from \pkg{qtl} package,

}

\examples{

  data(QTLexample)

  fullsib <- create_fullsib(fullsib_data,
                            list(LG1.end, LG2.end, LG3.end, LG4.end), 
                            step=0, map.function="kosambi", condIndex = 3.5) 



  ## running im_scan:
  im1 <- im_scan(fullsib, lg="all", pheno.col=1, LOD=TRUE)
  ## results in LOD Score for all linkage groups

  im2 <- im_scan(fullsib, lg=1:4, pheno.col=1, LOD=FALSE)
  ## results in -log10(pvalue) for all linkage groups


  plot(im1) ##default
  ##customizing graphics
  plot(im1, lty=3, lwd=2, incl.mkr="points")
  plot(im1, lty=3, lwd=2, incl.mkr="name")
  plot(im1, lty=3, lwd=2, incl.mkr="name", cex.incl=0.5)
  plot(im1, im2, label.lg=c("1st LG", "2nd LG","3rd LG","4th LG"))
  plot(im1, im2, label.lg=c("1st LG", "2nd LG","3rd LG","4th LG"),
  alternate.lgid=TRUE)

  #add argument is another way of overlapping QTL results
  plot(im1,lg=2, incl.mkr="name")
  plot(im2,lg=2, incl.mkr="name", add=TRUE)

}


\author{
  Mainly writted by Karl Broman and present at \code{plot.scanone} and
  adapted by Rodrigo Gazaffi \email{rgazaffi@gmail.com}}

\keyword{plot}
