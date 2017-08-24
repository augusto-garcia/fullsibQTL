\name{example_QTLfullsib}
\alias{example_QTLfullsib}
\alias{example_QTLfullsib.raw}
\alias{LG1_final}
\alias{LG2_final}
\alias{LG3_final}
\alias{LG4_final}
\alias{im_perm}
\alias{cim_perm}
\alias{example_QTLfullsib}

\docType{data}

\title{Simulated data set to perform QTL mapping in a full-sib cross}

\description{
  Tutorial data set to follow the \pkg{fullsibQTL} tutorial.
}

\usage{data(example_QTLfullsib)}

\format{
  This data loads four objects of class \emph{sequence}
  (\dQuote{LG1_final}, \dQuote{LG2_final}, \dQuote{LG3_final} and
  \dQuote{LG4_final}), i.e., they are the ordered linkage group and the
  final result of \pkg{onemap} package).

  Also, on example directory of this package, one may find a example
  file \dQuote{example_QTLfullsib.raw} which contains 65 markers with several
  segregation pattern and two traits (\dQuote{pheno1} and
  \dQuote{pheno2}) for 300 individuals.

  If the file \dQuote{example_QTLfullsib.raw} was read with function
  \code{read_onemap} from \pkg{onemap} package (>=2.0.7) the object
  \dQuote{example_QTLfullsib} will be created.
  
  Additionally, there are two permutation objects for the IM and CIM analysis (respectivetely, \dQuote{im_perm} and \dQuote{cim_perm}) both created with this same data using 1000 permutations and the functions \dQuote{im_scan} and \dQuote{cim_scan} (respectivetely, \dQuote{set.seed} equal to 1234 and 4321). More information at the  \pkg{fullsibQTL} tutorial.
}

\details{
  To exemplify the functions present in this package, it is considered a
  simulation done in the article of Gazaffi et al. (2014).
  
  In a summary, \dQuote{example_QTLfullsib} has 65 markers, with different
  segregation patterns: 15 of type A.1, 1 A.4, 5 B1.5, 6 B2.6, 6 B3.7, 11
  C.8, 10 D1.9, 10 D2.14, 1 D2.18. The codification of the segregations
  fashion follows the notation of Wu et al. (2002a). Those markers are
  distributed in four linkage groups, in each one having 15 equally
  distributed with 10 cM and five more unlinked markers that were
  simulated to demonstrate some options specially for cofactors
  selection. The population size is 300 individuals.

  The genetic map represented by these markers can be viewed in the
  objects \dQuote{LG1_final}, \dQuote{LG2_final}, \dQuote{LG3_final} and
  \dQuote{LG4_final}, that are objects of class \emph{sequence}.

  There are two simulated phenotypes that are present in
  \dQuote{example_QTLfullsib.raw} file. \dQuote{pheno1} is the original
  phenotype with 300 individuals and the trait heritability is
  0.70. \dQuote{pheno2} is the same trait as pheno1, but with an
  extra source of variation, modelled by addcovar parameter, along the
  functions to perform QTL mapping.
}

\references{
  Gazaffi, R., Mollinari, M., Pastina, M.M., Margarido, G.R.A, Garcia,
  A.A.F. (2014) QTL mapping model for characterizing segregation
  pattern and linkage phase for outcrossing species. \emph{Tree Genetics
  and Genomes}, 10(4), 791-801
  
  Margarido, G. R. A., Souza, A.P. and Garcia, A. A. F. (2007)
  OneMap: software for genetic mapping in outcrossing species.
  \emph{Hereditas} 144: 78-79

  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002a)
  Simultaneous maximum likelihood estimation of linkage and linkage
  phases in outcrossing species. \emph{Theoretical Population Biology}
  61: 349-363

  Wu, R., Ma, C.-X., Wu, S. S. and Zeng, Z.-B. (2002b). Linkage
  mapping of sex-specific differences. \emph{Genetical Research} 79:
  85-96
}
  
\examples{
  data(example_QTLfullsib)
  ls()

  fullsib <- create_fullsib(example_QTLfullsib,
                            list(LG1_final, LG2_final, LG3_final, LG4_final),
                            step=0,map.function="kosambi",condIndex=3.5)

  ##im analysis
  im1 <- im_scan(fullsib, pheno.col=1)

  \dontrun{
  plot(im1)
  summary(im1)
  print(im1, lg=4)
  }

  qtl1 <- im_char(fullsib, pheno.col=1, lg=1, pos="M1")

  ##cofactor definition
  cofs <- matrix(c("1","M1",
                   "2","M27",
                   "3","M37",
                   "4","M52"),4,2, byrow=TRUE)
  (cof.def <- cof_definition(fullsib, pheno.col=1, cof.pos=cofs))

  ##cofactor selection
  cofs.fs <- cof_selection(fullsib, pheno.col = 1, k = log(300),
  selection = 1)
  cofs.fs
  plot(cofs.fs)

  ##cim analysis
  cim1 <- cim_scan(cofs.fs, pheno.col=1, ws = 22, LOD= TRUE, icim=FALSE)
  \dontrun{
  plot(cim1)
  summary(cim1)
  print(cim1, lg=4)
  }

  qtl.lg3 <- cim_char(cofs.fs, pheno.col=1, ws=22, lg=3, pos="M38")
  draw_phase(cofs.fs, qtl.lg3, probs=0.05)
  get_segr(qtl.lg3)

  qtl.lg4 <- cim_char(cofs.fs, pheno.col=1, ws=22, lg=4, pos="M52")
  draw_phase(cofs.fs, qtl.lg4, probs=0.05)
  get_segr(qtl.lg4)

  r2_ls(fullsib, pheno.col=1, lg=c(3,4), pos=c("M38","M52"),
  ls.estimation=TRUE)
}

\seealso{
  \code{\link[fullsibQTL]{create_fullsib}},
  \code{onemap}

}

\keyword{datasets}


\author{Rodrigo Gazaffi, \email{rgazaffi@gmail.com}}
