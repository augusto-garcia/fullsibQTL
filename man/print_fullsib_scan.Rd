\name{print_fullsib_scan}

\alias{print_fullsib_scan}

\title{Prits object of class \emph{fullsib_scan}}

\description{
  Prints the object of class \emph{fullsib_scan} created by \code{im_scan}
  or \code{cim_scan}.
}

\usage{
\method{print}{fullsib_scan}(x, lg, \ldots)
}

\arguments{
  \item{x}{Object of class \emph{fullsib_scan} resulted from
  \code{im_scan} or \code{cim_scan}.}

  \item{lg}{Linkage group that will be printed. Default is to print all
  scanned groups.}

  \item{\ldots}{Further arguments, passed to other methods. Currently
  ignored.}
}

\value{None}

\seealso{
  \code{\link[fullsibQTL]{im_scan}}, 
  
  \code{\link[fullsibQTL]{cim_scan}}, 
  
  \code{\link[fullsibQTL]{plot_fullsib_scan}}, 
  
  \code{\link[fullsibQTL]{summary_fullsib_scan}}, 
  
  \code{scanone} from \pkg{qtl} package
}

\examples{

  data(QTLexample)
  fullsib <- create_fullsib(fullsib_data,
                            list(LG1.end, LG2.end, LG3.end, LG4.end), 
                            step = 0, map.function = "kosambi",
                            condIndex = 3.5)



  ## running im_scan:
  im1 <- im_scan(fullsib, lg = "all", pheno.col = 1, LOD = TRUE)

  im1
  print(im1, lg=2)
  print(im1, lg=c(1,4))
}

\keyword{utilities}


\author{Rodrigo Gazaffi, \email{rgazaffi@gmail.com}}