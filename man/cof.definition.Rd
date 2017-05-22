\name{cof.definition}

\alias{cof.definition}
\alias{plot.fullsib.cofactors}

\title{Ad-hoc definition of cofactors for CIM analysis}

\description{  
  Defines any marker and/or position to be used as a cofactor. It was
  designed as an alternative way to include cofactor besides using
  \code{cof.selection} function.
}

\usage{
cof.definition(fullsib, pheno.col=1, addcovar=NULL, cof.pos,
thres.effect=1)

\method{plot}{fullsib.cofactors}(x, horizontal = FALSE, grid = FALSE,
cex.mrk = 1, cex.grp = 0.75,\ldots)

}

\arguments{
  \item{fullsib}{An object from class \emph{fullsib}.}
  
  \item{pheno.col}{Column number in the phenotype matrix (present in
    \emph{fullsib} object) which should be used as the phenotype.}

  \item{addcovar}{Additive covariates. If it is used, one must indicate the
    design matrix for those source of variation. It should be noted that
    additive covariates is included in the model as fixed effects, under
    ordinary linear regression.}

  \item{cof.pos}{character matrix, in each row represents a
    different cofactor and two columns indicating the linkage group and
    the name of the position to be considered. To define markers that
    are not place on markers, one should indicate as \code{NA}. See
    examples.}
  
  \item{thres.effect}{Threshold value to remove non significative effect
    for selected cofactors. Default is \code{thres.effect = 1}, none
    cofactors effect is removed from the analysis.}

  \item{x}{An object from class \emph{fullsib.cofactors}.}
  
  \item{horizontal}{if \code{TRUE}, indicates that the map should be plotted
    horizontally. Default is \code{FALSE}.}
  
  \item{grid}{if \code{TRUE}, displays a grid in the background. Default is
    \code{FALSE}.}
  
  \item{cex.mrk}{the magnification to be used for markers.}
  
  \item{cex.grp}{the magnification to be used for group axis
    annotation.}

  \item{\ldots}{Further arguments, passed to other methods. Currently
    ignored.}
}

\details{
  The standard procedure used to perform cofactor selection on this
  package is based on multiple linear regression using information
  criteria. However for any reason, one would like to have more
  flexibility to control the number and location of cofactors to be
  added in the model.

  One possible way of select cofactors is to perform an interval mapping
  analysis and with the result, one can perform CIM analysis, i.e.,
  first IM is done to detect QTL. After this, the position of mapped QTL
  can be added for CIM analysis, using \code{cof.definition}. The
  inclusion and exclusion of cofactor can be performed for some rounds,
  until one get a final QTL mapping profile.

  The method plot was designed for the user evaluated the cofactors
  saturation and dispersion along the genome. Helping with the decision
  of window size dimension and best selection options.
  
  Finally, with the development of \code{cof.definition} and
  \code{cof.selection}. We believe that user has enough flexibility of
  dealing with selection of markers to be used in CIM.
  
}

\value{
  An object of class \emph{fullsib.cofactors} returned, which has the
  same structure of an object of class \emph{fullsib} with the inclusion
  of an extra component (\sQuote{cofactor}) that is a list with the
  components: \code{names.cof}, \code{matrix.cof} and \code{trait.cof}.

  \code{names.cof} is a data frame showing which are the selected
  markers and their linkage groups. \code{NA} value is used for markers
  that are not placed on linkage groups.

  \code{matrix.cof} is a matrix with contains all the cofactors
  effect. During the CIM analysis columns are dropped in function of
  window size.

  \code{trait.cof} is the indication for which pheno.col was considered
  for cofactor selection step.
}

\author{Rodrigo Gazaffi, \email{rgazaffi@gmail.com}
  \code{plot.fullsib.cofactors} was based on \code{draw.map} that was
  written by Marcelo Mollinari \email{mmollina@esalq.usp.br} and modified
  by Rodrigo Gazaffi, allowing cofactors inclusion.
}

\seealso{
  \code{\link[fullsibQTL]{create.fullsib}},
  
  \code{\link[fullsibQTL]{cof.selection}}

  \code{\link[fullsibQTL]{cim.scan}}

  \code{\link[fullsibQTL]{cim.char}}
}
\examples{

  data(QTLexample)

  fullsib <- create.fullsib(fullsib.data,
                            list(LG1.end, LG2.end, LG3.end, LG4.end),
                            step=0,map.function="kosambi",condIndex=3.5)

 \dontrun{
  im1 <- im.scan(fullsib, pheno.col=1, LOD = TRUE)
  summary(im1)
  }

  ## using 4 QTL as cofactors (QTL peaks detected using im.scan)
  cofs <- matrix(c("1","M2",
                   "2","M27",
                   "3","M37",
                   "4","M52"),4,2, byrow=TRUE)

  cof.def <- cof.definition(fullsib, pheno.col=1, cof.pos=cofs)
  cof.def
  plot(cof.def)

  \dontrun{
  covar <- matrix(rep(c(1,-1), each=150), ncol=1)
 
  ##using 4 QTL as cofactors + 1 unlinked marker (just for illustration)
  cofs2 <- matrix(c("1","M2",
                    "2","M27",
                    "3","M37",
                    "4","M52",
                    NA, "M64"),5,2, byrow=TRUE)

  cof2.def <- cof.definition(fullsib, pheno.col=2, addcovar=covar, cof.pos=cofs2, thres.effect=1)
  cof2.def

  cof3.def <- cof.definition(fullsib, pheno.col=1, cof.pos=cofs2, thres.effect=0.05)
  cof3.def
  ### realize with thres.effect = 0.05, the 5th cofactor is removed,
  ### because it is non significative (never was selected on cof.selection)
  }

}
\keyword{utilities}
