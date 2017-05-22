\name{print.fullsib.scan}

\alias{print.fullsib.scan}

\title{Prits object of class \emph{fullsib.scan}}

\description{
  Prints the object of class \emph{fullsib.scan} created by \code{im.scan}
  or \code{cim.scan}.
}

\usage{
\method{print}{fullsib.scan}(x, lg, \ldots)
}

\arguments{
  \item{x}{Object of class \emph{fullsib.scan} resulted from
  \code{im.scan} or \code{cim.scan}.}

  \item{lg}{Linkage group that will be printed. Default is to print all
  scanned groups.}

  \item{\ldots}{Further arguments, passed to other methods. Currently
  ignored.}
}

\value{None}

\seealso{
  \code{\link[fullsibQTL]{im.scan}}, 
  
  \code{\link[fullsibQTL]{cim.scan}}, 
  
  \code{\link[fullsibQTL]{plot.fullsib.scan}}, 
  
  \code{\link[fullsibQTL]{summary.fullsib.scan}}, 
  
  \code{scanone} from \pkg{qtl} package
}

\examples{

  data(QTLexample)
  fullsib <- create.fullsib(fullsib.data,
                            list(LG1.end, LG2.end, LG3.end, LG4.end), 
                            step = 0, map.function = "kosambi",
                            condIndex = 3.5)



  ## running im.scan:
  im1 <- im.scan(fullsib, lg = "all", pheno.col = 1, LOD = TRUE)

  im1
  print(im1, lg=2)
  print(im1, lg=c(1,4))
}

\keyword{utilities}


\author{Rodrigo Gazaffi, \email{rgazaffi@gmail.com}}