\name{fullsibQTL-package}

\alias{fullsibQTL-package}

\alias{fullsibQTL}

\docType{package}

\title{
  Software for QTL mapping using full-sib cross or F1 population
}

\description{
  The package \pkg{fullsibQTL} was designed to perform QTL mapping in a
  diploid biparental cross (full-sib progeny or \eqn{F_1} population),
  where an integrated genetic map was previously obtained with molecular
  markers segregating for 1:1:1:1, 1:2:1, 3:1 or 1:1 fashion. We
  consider that genetic map was previously obtained with \pkg{onemap}
  package. To perform QTL mapping, we develop \pkg{fullsibQTL} which
  considers QTL mapping approach under mixture models and maximum
  lihelihood approach. The multipoint probabilities for QTL genotypes
  are obtained by Hidden Markov Models. We implemented Interval Mapping
  (IM) and Composite Interval Mapping (CIM). Briefly, after the mapping
  procedure, the software shows position on the genome, the estimable
  effects, including LOD Score to verify if they are significative, the
  QTL segregation pattern on the progeny, and its linkage phase with
  flanking markers. This methodology was present by Gazaffi et
  al. (submitted). 
}

\details{

  \tabular{ll}{
    Package: \tab fullsibQTL\cr
    Type: \tab Package\cr
    Version: \tab 1.0-0\cr
    Date: \tab 2011-09-30\cr %update
    License: \tab GNU GENERAL PUBLIC LICENSE (Version 3)\cr
  }
  
  Usage of \pkg{fullsibQTL} is described in a tutorial distributed with
  the package. One can look for it in directory \code{/doc} of the
  package distribution. 

  The most important functions are:
  \itemize{
    \item \code{\link[fullsibQTL]{create_fullsib}} creates the main
    object to perform QTL mapping, including the calculation for QTL
    genotypes using multipoint technology.

    \item \code{\link[fullsibQTL]{im_scan}} scan the genome or selected
    linkage groups for QTL using Interval Mapping method.
    
    \item \code{\link[fullsibQTL]{im_char}} for a given position, QTL is
    characterized, i.e., the genetic effects are estimated and tested if
    they are significative, which allows to estimate QTL segregation
    patterns and its linkage phases with flanking markers.

    \item \code{\link[fullsibQTL]{cim_scan}} scan the genome or selected
    linkage groups for QTL using COMPOSITE Interval Mapping method.

    \item \code{\link[fullsibQTL]{cim_char}} characterizes the QTL mapped
    for a given position as similar done in \code{im_char}, but using
    CIM analysis.

    \item \code{\link[fullsibQTL]{cof.selection}} performs cofactor
    selection using multiple linear regression using an information
    criteria (AIC or BIC). This function is a previous step for
    \code{cim_scan} or \code{cim_char}.

    \item \code{\link[fullsibQTL]{cof_definition}} allows to define any
    region of the genome to be considered as cofactors for CIM analysis,
    instead using of linear regression approach. This function is a
    previous step for \code{cim_scan} or \code{cim_char}.
     
    \item \code{\link[fullsibQTL]{draw_phase}} receives the object of
    class \emph{fullsib_char} (output from \code{im_char} or
    \code{cim_char}) and returns the linkage groups with QTL mapped,
    allowing to infer the linkage phase between QTL and markers.

    \item \code{\link[fullsibQTL]{get_segr}} using an object of class
    \emph{fullsib_char}, returns the QTL segregation for a mapped QTL.
    
    \item \code{\link[fullsibQTL]{r2_ls}} considering all QTL mapped for a
    single trait, it calculates the amount of phenotypic variation that
    all QTL explain jointly, and also its values for each QTL
    separately. Note: To obtain these estimations, we used least square 
    approach, once multiple interval mapping is not the focus of this
    package. 
  }
}

\note{
  The complete tutorial is distributed with the package, as well a
  script with the example.
}

\author{
  Rodrigo Gazaffi, \email{rgazaffi@gmail.com},
  Antonio Augusto F Garcia, \email{augusto.garcia@usp.br}
}

\references{

  Gazaffi, R, Mollinari, M, Pastina, MM, Margarido, GRA, Garcia, AAF
  (submitted) A model for quantitative trait loci mapping, linkage phase
  and segregation pattern estimation for outcrossing species. \emph{Tree
  Genetics and Genomes}

  Jiang, C, Zeng, Z-B (1997) Mapping quantitative trait loci with
  dominant and missing markers in various crosses from two inbred
  lines. \emph{Genetica} 101: 47-58.

  Margarido, GRA, Souza, AP, Garcia, AAF (2007) OneMap:software for
  genetic mapping in outcrossing species. \emph{Hereditas}144: 78-79.

  Wu, R, Ma, C-X, Painter, I, Zeng, Z-B (2002a) Simultaneous maximum
  likelihood estimation of linkage and linkage phases in outcrossing
  species. \emph{Theoretical Population Biology} 61: 349-363.

  Wu, R, Ma, C-X, Wu, SS, Zeng, Z-B (2002b). Linkage mapping of
  sex-specific differences. \emph{Genetical Research} 79: 85-96
}

\keyword{package}


