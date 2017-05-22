\name{plot.fullsib.perm}

\alias{plot.fullsib.perm}

\title{Plot the histogram of maximum values of permutation test}

\description{Plots the histogram of maximum values obtained with the
  permutation test. The graphic can be done for each QTL peak}

\usage{ \method{plot}{fullsib.perm}(x, peak=1, \ldots) }

\arguments{
  \item{x}{An object of class \emph{fullsib.perm}, output of
    \code{im.scan} or \code{cim.scan}.}
  
  \item{peak}{QTL peak to be used to plot the histogram. Default is to
    use the first peak. See details.}
  
  \item{\dots}{Passed to the function \code{plot} when it is called. See
  \code{hist} for more options to be included here.}


}

\details{
  When permutation test is done, an empirical distribuition is obtained
  for having a threshold for QTL mapping. However, we implemented the
  approach of Chen and Storey (2006) for having a relaxed threshold
  values. As a result, the object of class \emph{fullsib.perm} contains
  the distribution of two highest peaks, each one can be assessed using
  the argument. If \sQuote{peak=1} shows the classical permutation test
  present by Churchill and Doerge (1994).
}

\value{
  A histogram showing the distribution of maximum values obtained on
  each permutation test. \code{hist} function is used to do the plot.
}

\references{
  Chen, L., Storey, J.D. (2006) Relaxed Significance Criteria for
  Linkage Analysis. \emph{Genetics} 173: 2371-2381

  Churchill, G.A., Doerge, R.W. (1994) Empirical Threshold Values for
  Quantitative Trait Mapping. \emph{Genetics} 138: 963-971.
}

\seealso{
  \code{\link[fullsibQTL]{im.scan}},
  
  \code{\link[fullsibQTL]{cim.scan}},
  
  \code{\link[fullsibQTL]{summary.fullsib.perm}},

  \code{hist},  
}


\examples{
  \dontrun{
  data(QTLexample)

  fullsib <- create.fullsib(fullsib.data,
                            list(LG1.end, LG2.end, LG3.end, LG4.end), 
                            step=0, map.function="kosambi", condIndex = 3.5) 

  im.perm <- im.scan(fullsib, pheno.col=1, LOD=TRUE, n.perm=1000,
                     write.perm="perm_values.txt") 
  plot(im.perm, peak = 1) ## churchill and doerge (1994)
  plot(im.perm, peak = 2) ## chen and storey (2006)
  }
}

\author{Written by Karl Broman for \sQuote{plot.scanoneperm} from
  \pkg{qtl} package and adapted by Rodrigo Gazaffi  
  \email{rgazaffi@gmail.com}}

\keyword{plot}

