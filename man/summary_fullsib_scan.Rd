\name{summary_fullsib_scan}

\alias{summary_fullsib_scan}

\title{Summarize QTL mapping scan}

\description{
  Summarizes the QTL mapping scan resulted from the \code{im_scan} and
  \code{cim_scan}.
}

\usage{\method{summary}{fullsib_scan}(object, thr=0,\ldots)}

\arguments{
  \item{object}{An object of class \emph{fullsib_scan} (output of
    \code{im_scan} and \code{cim_scan}.}
  
  \item{thr}{threshold value to declare a QTL, i.e., only prints peaks
    with value higher than \sQuote{thr}.}
    
  \item{\ldots}{Further arguments, passed to other methods. Currently
    ignored.}
}

\value{
  It gets the object of class \emph{fullsib_scan} and simplifies in a
  matrix containing just the maximum peak of each linkage group that
  exceeds the threshold value defined by the user. 
}

\details{
  This function was designed to get only the maximum peak of each
  analysed group. For the cases that more than one peak is present
  should decided the position of the QTL examinating object of class
  \emph{fullsib_scan} by oneself. We assumed this, because to define a
  QTL can be very subjective in some situations.
}

\seealso{
  \code{\link[fullsibQTL]{im_scan}}, 

  \code{\link[fullsibQTL]{cim_scan}}, 

  \code{\link[fullsibQTL]{im_char}}, 

  \code{\link[fullsibQTL]{cim_char}}, 

  \code{\link[fullsibQTL]{plot_fullsib_scan}}, 

  \code{\link[fullsibQTL]{summary_fullsib_scan}}, 

  \code{scanone} from \pkg{qtl} package,

  \code{summary.scanone} from \pkg{qtl} package
}
\examples{
  data(QTLexample)

  fullsib <- create_fullsib(fullsib_data,
                            list(LG1.end, LG2.end, LG3.end, LG4.end), 
                            step=0, map.function="kosambi", condIndex = 3.5) 

  im1 <- im_scan(fullsib, lg="all", pheno.col=1, LOD=TRUE)
  summary(im1)
  summary(im1, thr=6.5)

}

\author{Rodrigo Gazaffi \email{rgazaffi@gmail.com}}
\keyword{utilities}
