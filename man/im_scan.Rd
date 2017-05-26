\name{im_scan}

\alias{im_scan}

\title{QTL searching using Interval Mapping approach}

\description{
  Performs QTL mapping with Interval Mapping (IM) approach by mixture
  models. MLE are obtained using EM algorithm. Permutation test is also
  implemented for threshold determination.
}

\usage{
im_scan (fullsib, lg, pheno.col=1, addcovar=NULL, LOD=TRUE, maxit=1000,
         tol = 1e-04, n.perm, write.perm, verbose, verbose.scan) 
}

\arguments{
  \item{fullsib}{An object from class \emph{fullsib}.}
  
  \item{lg}{Vector indicating which linkage group will be
    scanned. Default is to analyse all groups.} 
  
  \item{pheno.col}{Column number in the phenotype matrix (present in
    \emph{fullsib} object) which should be used as the phenotype.}
  
  \item{addcovar}{Additive covariates. If it is used, one must indicate the
    design matrix for those source of variation. It should be noted that
    additive covariates is included in the model as fixed effects, under
    ordinary linear regression.}
  
  \item{LOD}{if \code{TRUE} indicates the mapping result as LOD Score,
    if \code{FALSE} QTL search is presents as \eqn{-\log_{10}(pvalue)}.}
  
  \item{maxit}{Maximum number of iteration in EM algorithm.}
  
  \item{tol}{Tolerance for determining convergence in EM algorithm. See
  details.}

  \item{n.perm}{If \code{n.perm} is missing or zero, usual genome scan
    is performed, otherwise, permutation test is done. The number of
    permutation to be used is defined for \code{n.perm} value}

  \item{write.perm}{Optional string (e.g., \dQuote{file.txt}) to create
    a text file containing the results obtained on the permutation
    test.}
  
  \item{verbose}{If \code{TRUE} display information during EM
    algorithm. For each position on the genome, it indicates the
    iteration of EM is performing, with log-likelihood, convergence,
    genetics effects, square root of variance. The 
    log-likelihood of the model under \eqn{H_0} model is also showed.}
  
  \item{verbose.scan}{If \code{TRUE} display a progress bar showing the
    percentage of genome scan is done, based on the number of linkage
    groups to be scanned. However, when permutation test is done, the
    option is coded as FALSE for having a clear output.}

}

\details{
  
  This function receives the object of class \emph{fullsib} and scan the
  genome for QTL, using interval mapping. For each position, EM
  algorithm is performed considering the following criteria for
  convergence:
  
  \deqn{conv = abs\left[\frac{(new.lk - old.lk)}{old.lk}\right]}
  
  in which, \eqn{old.lk} is the likelihood for \eqn{i^{th}} iteration and
  \eqn{new.lk} is the likelihood for \eqn{(i+1)^{th}}.
      
  In Chen and Storey's (Chen and Storey, 2006) article of having relaxed
  threshold values, they showed different way of finding a QTL
  peak. However to define a peak can be very trick and subjective, so
  for this function, peak is the maximum value of LOD or
  \eqn{-\log_{10}(pvalue)} can assume in each linkage group, for
  example, on fullsib_data (data set contained in this package) has only
  four groups, so we consider that we can define 4 peaks. In cases that
  dataset has more than 10 linkage groups, we restrict to storage until
  the 2nd highest peaks.
  
  The permutation test can be done for this \sQuote{n.perm} argument is
  higher then 0, but if one wants to repeat the result of permutation,
  please used \code{set.seed} function from \R before run \code{im_scan}.

  Individuals with missing phenotypes are dropped.

  
}

\value{
  This function can return two differents classes of objects
  \dQuote{fullsib_scan} and \dQuote{fullsib_perm}:

  If usual analysis is done (\sQuote{n.perm=0}), the function returns a
  matrix with 4 columns: The first two columns indicates the linkage
  groups number and position in centiMorgan the analysis has been done
  (\code{lg} and \code{pos.cM}, respectivaly). The third is the value
  of QTL mapping test, one can choose between \code{LOD} or
  \eqn{-\log_{10}(pvalue)}. The fourth indicates which model was
  considered based on the level of informativeness of the probabilities.
  These values can vary between 0 to 9. For further details about this
  see the \pkg{fullsibQTL} Vignette.  The names of each row is important
  for QTL characterization, because it will be used as argument for
  \code{cim_char} and \code{r2_ls} functions.

  For situation that \sQuote{nperm} is higher then zero, permutation
  test is performed and the function returns an object from
  \dQuote{fullsib_perm}. This object is a matrix which the number of row
  indicates the number of permutations and the number of column defines
  the number of highest peaks colected on each permutation. As described
  on details, we define peak the maximum value of either \code{LOD} or
  \eqn{-\log_{10}(pvalue)} of each linkage group.  They are
  ordered from the maximum to minimum value. If one determines the
  threshold with the distribuition of the highest peak (column one), one
  have the threshold as present by Churchill and Doerge (1994). If one
  choose the threshold with the second peak or lower is considering the
  approach by Chen and Storey (2006). The function
  \emph{summary_fullsib_perm} provides automatically these results.

}

\references{
  Chen, L., Storey, J.D. (2006) Relaxed Significance Criteria for
  Linkage Analysis. \emph{Genetics} 173: 2371-2381

  Churchill, G.A., Doerge, R.W. (1994) Empirical Threshold Values for
  Quantitative Trait Mapping. \emph{Genetics} 138: 963-971.

  Gazaffi, R., Mollinari, M., Pastina, M.M., Margarido, G.R.A, Garcia,
  A.A.F. (submitted) A model for quantitative trait loci mapping,
  linkage phase and segregation pattern estimation for outcrossing
  species. \emph{Tree Genetics and Genomes}

}

\seealso{
  \code{\link[fullsibQTL]{create_fullsib}},

  \code{\link[fullsibQTL]{im_char}},
   
  \code{\link[fullsibQTL]{plot_fullsib_scan}},

  \code{\link[fullsibQTL]{summary_fullsib_scan}},

  \code{\link[fullsibQTL]{plot_fullsib_perm}},
  
  \code{\link[fullsibQTL]{summary_fullsib_perm}}

} 

\examples{

  data(QTLexample)

  fsib <- create_fullsib(fullsib_data,
                         list(LG1.end, LG2.end, LG3.end, LG4.end), 
                         step=0, map.function="kosambi", condIndex = 3.5) 


  ## running im_scan:

  im1 <- im_scan(fsib, lg="all", pheno.col=1, LOD=TRUE)
  ## results in LOD Score for all linkage groups

  im2 <- im_scan(fsib, lg=1:4, pheno.col=1, LOD=FALSE)
  ## results in -log10(pvalue) for all linkage groups

  ## im_scan with using additive covariate
  covar <- matrix(rep(c(1,-1), each=150), ncol=1)
  im3 <- im_scan(fsib, lg=c(2,3), pheno.col=2, addcovar=covar, LOD=FALSE)

  \dontrun{
  ## RESULTS:
  im1
  print(im1,lg=3) ## printing only one linkage group

  ## detecting QTL peaks
  summary(im1)
  summary(im1, thr=6.5)


  ##plotting
  plot(im1)

  plot(im1,im2, col=c("black","blue"), label.lg=c("A","B","C","D"))
  abline(3,0) ##threshold
  legend("topleft",legend=c("LOD","-log10(pvalue)"),
          col=c("black","blue"), lwd=c(2,2), text.col=c("black","blue"))

  plot(im2,lg=2)
  plot(im2,lg=2, incl.mkr="points")
  plot(im2,lg=2, incl.mkr="name")
  }

  ## permutation test:
  \dontrun{
  im.perm <- im_scan(fsib, lg=c(2,4), pheno.col=1, LOD=TRUE,
                     n.perm=1000, write.perm="perm_values.txt")
  ## permutations done just using linkage groups 2 and 4
  summary(im.perm)
  plot(im.perm, peak=1)
}

}

\author{Rodrigo Gazaffi, \email{rgazaffi@gmail.com}}

\keyword{utilities}

