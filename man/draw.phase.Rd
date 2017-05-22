\name{draw.phase}

\alias{draw.phase}

\title{Represents the linkage phase between markers and QTL}

\description{
  Locate the QTL in a linkage group, representing the significative
  parental effects from the result of \code{im.char} and \code{cim.char}
}

\usage{
 draw.phase(fullsib, fschar, probs=0.05)
}

\arguments{
  \item{fullsib}{An object from class \emph{fullsib} or
  \emph{fullsib.cofactors}.} 
  \item{fschar}{An object from class \emph{fullsib.char}.}
  \item{probs}{probability used to identify significative parental
  effect.}

}

\value{
  It returns text output indicating in a linkage group linkage phase
  between markers and QTL.}

\author{ 
  Based on \code{print.sequence} from \pkg{onemap} and it was written by
  Marcelo Mollinari \email{mmollina@esalq.usp.br}. It was modified
  allowing the inclusion of QTL, by Rodrigo Gazaffi
  \email{rgazaffi@gmail.com} 
}

\seealso{
  \code{\link[fullsibQTL]{im.char}},
  
  \code{\link[fullsibQTL]{cim.char}},
  
  \code{\link[fullsibQTL]{get.segr}},
  
}

\examples{

  data(QTLexample)

  fullsib <- create.fullsib(fullsib.data,
                            list(LG1.end, LG2.end, LG3.end, LG4.end),
                            step=0,map.function="kosambi",condIndex=3.5)


  ###############################################
  ## cofactor selection using BIC (n.ind = 300)
  cofs.fs <- cof.selection(fullsib, pheno.col=1, k = log(300),
                           selection=1) 

  \dontrun{
  cim1 <- cim.scan(cofs.fs, pheno.col=1, ws = 22, LOD= TRUE, icim=FALSE)
  summary(cim1)
  plot(cim1)
  }

  qtl.lg3 <- cim.char(cofs.fs, pheno.col=1, ws=22, lg=3, pos="M38")
  draw.phase(cofs.fs, qtl.lg3, probs=0.05)

  qtl.lg4 <- cim.char(cofs.fs, pheno.col=1, ws=22, lg=4, pos="M52")
  draw.phase(cofs.fs, qtl.lg4, probs=0.05)
}
\keyword{utilities}
