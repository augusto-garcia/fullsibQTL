\name{cim.scan}

\alias{cim.scan}

\title{QTL searching using Composite Interval Mapping}

\description{
  Performs QTL mapping using Composite Interval Mapping (CIM) approach,
  with mixture models. MLE are obtained using EM algorithm. Permutation test
  is also implemented for threshold determination
}

\usage{
cim.scan(fullsib, lg, pheno.col=1, ws = 10, LOD=TRUE, maxit=1000, tol=1e-04,
         n.perm, write.perm, icim=FALSE, verbose, verbose.scan) 
}

\arguments{
  \item{fullsib}{An object from class \emph{fullsib.cofactors}.}
  
  \item{lg}{Vector indicating which linkage group will be
    scanned. Default is to analyse all groups.}
  
  \item{pheno.col}{Column number in the phenotype matrix (present in
    \emph{fullsib.cofactors} object) which should be used as the phenotype.}
  
  \item{ws}{Window Size in cM. Default is 10 cM (i.e., 5 cM for each
    side.}
  
  \item{LOD}{if \code{TRUE} indicates the mapping result as LOD Score,
    if \code{FALSE} QTL search is presents as \eqn{-\log_{10}(pvalue)}.}

  \item{maxit}{Maximum number of iteration in EM algorithm.}
  
  \item{tol}{Tolerance for determining convergence in EM algorithm. See
    details.}
   
  \item{n.perm}{If \code{n.perm} is missing or zero, usual genome scan
    is performed, otherwise, permutation test is done. The number of
    permutation to be used is defined for \code{n.perm} value.}

  \item{write.perm}{Optional string (e.g., \dQuote{file.txt}) to create
    a text file containing the results obtained on the permutation
    test.}
  
  \item{icim}{if \code{TRUE}, icim approach proposed by Li et al. (2007)
    is done (extended for full sib progeny), if \code{FALSE} (default)
    traditional CIM is performed. See details for icim method.}

  \item{verbose}{If \code{TRUE} display information during EM
    algorithm. For each position on the genome, it indicates the 
    iteration of EM is performing, with log-likelihood, convergence,
    genetics effects, square root of variance. The log-likelihood of the
    model under \eqn{H_0} model is also showed.}
  
  \item{verbose.scan}{If \code{TRUE} display a progress bar showing the
    percentage of genome scan is done, based on the number of linkage
    groups to be scanned. However, when permutation test is done, the
    option is coded as FALSE for having a clear output.}

 
}

\details{
  The mapping procedure is performed similar as present in
  \code{im.scan}. However the inclusion of cofactors can locate QTL with
  more precision with higher statistical power in comparison of interval
  mapping approach.

  The usage of the cofactors on the CIM is function of window
  size. Here, we consider window size minimal distance that cofactors
  should be from the interval under analysis. For example, if one
  considers a window size of 20 cM, while CIM is being performed, any
  marker located closer than 10 cM for each side, will not be included
  in the analysis.

  If extra additive covariates is used, one need to specify them, on
  \code{cof.selection} or \code{cof.definition} step.
  
  Li et al. (2007) proposed an adaption for CIM first developed by Zeng
  (1994). In their metodology, first the phenotype should be regressed
  to remove the variation due to cofactors and intercept. On a second
  moment, the residual is used to map QTL, similar an interval mapping
  approach. The method was proposed for inbred lines, and we included
  this method on CIM done in fullsib context.

  For convergence during EM iteration we use the following criterion:
  \deqn{conv = abs\left[\frac{(new.lk - old.lk)}{old.lk}\right]}
  in which, \eqn{old.lk} is the likelihood for \eqn{i^{th}} iteration and
  \eqn{new.lk} is the likelihood for \eqn{(i+1)^{th}}. If convergence
  value was higher than \code{tol} argument, iterative process continues, if
  not iteration is stopped.
  
  We also implemented permutation test on \code{cim.scan}, including the
  approach proposed by Chen and Storey (Chen and Storey, 2006). In this
  case, to define peak we used the maximum value of LOD or
  \eqn{-\log_{10}(pvalue)} in each linkage group. In fullsib.data
  (example dataset) there are four linkage groups, so we consider four
  peaks, if a dataset has more than 10 linkage groups, we restrict to
  storage until the 10th highest peaks.

  Individuals with missing phenotypes are dropped.

}

\value{

  This function can return two differents classes of objects
  \dQuote{fullsib.scan} and \dQuote{fullsib.perm}:

  If usual analysis is done (\sQuote{n.perm=0}), the function returns a
  matrix with 4 columns: The first two columns indicates the linkage
  groups number and position in centiMorgan the analysis has been done
  (\code{lg} and \code{pos.cM}, respectivaly). The third is the value
  of QTL mapping test, one can choose between \code{LOD} or
  \eqn{-\log_{10}(pvalue)}. The fourth indicates which model was
  considered based on the level of informativeness of the probabilities.
  These values can vary between 0 to 9. For further details about this
  see the \pkg{fullsibQTL} Vignette. The names of each row is important
  for QTL characterization, because it will be used as argument for
  \code{cim.char} and \code{r2.ls} functions.

  For situation that \sQuote{nperm} is higher then zero, permutation
  test is performed and the function returns an object from
  \dQuote{fullsib.perm}. This object is a matrix which the number of row
  indicates the number of permutations and the number of column defines
  the number of highest peaks colected on each permutation. As described
  on details, we define peak the maximum value of either \code{LOD} or
  \eqn{-\log_{10}(pvalue)} of each linkage group.  They are
  ordered from the maximum to minimum value. If one determines the
  threshold with the distribuition of the highest peak (column one), one
  have the threshold as present by Churchill and Doerge (1994). If one
  choose the threshold with the second peak or lower is considering the
  approach by Chen and Storey (2006). The function
  \code{summary.fullsib.perm} or just \code{summary} provides
  automatically these results.
  
}

\references{
  Chen, L., Storey, J.D. (2006) Relaxed Significance Criteria for
  Linkage Analysis. \emph{Genetics} 173: 2371-2381
  
  Gazaffi, R., Mollinari, M., Pastina, M.M., Margarido, G.R.A, Garcia,
  A.A.F. (submitted) A model for quantitative trait loci mapping,
  linkage phase and segregation pattern estimation for outcrossing
  species. \emph{Tree Genetics and Genomes}
  
  Li, H., Ye, G.;  Wang, J. (2007) A Modified Algorithm for the
  Improvement of Composite Interval Mapping. \emph{Genetics} 175: 361-374
  
}

\author{
  Rodrigo Gazaffi, \email{rgazaffi@gmail.com}
}

\seealso{
  \code{\link[fullsibQTL]{create.fullsib}},

  \code{\link[fullsibQTL]{im.scan}},

  \code{\link[fullsibQTL]{cim.char}},

  \code{\link[fullsibQTL]{cof.selection}},

  \code{\link[fullsibQTL]{cof.definition}}
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

  cim1 <- cim.scan(cofs.fs, pheno.col=1, ws = 22, LOD= TRUE, icim=FALSE)

  \dontrun{
  covar <- matrix(rep(c(1,-1), each=150), ncol=1)
  cofs.fs <- cof.selection(fullsib, pheno.col=2, addcovar=covar, k=log(300))
  cim2 <- cim.scan(cofs.fs, pheno.col=2, ws = 22, LOD= TRUE, icim=TRUE) 
  ## cim with covariate, one just indicate the additive covar on
  ## cof.selection or cof.definition

  ##permutation test:
  cim.perm <- cim.scan(cofs.fs, pheno.col=1, ws = 22, LOD = FALSE,
                       n.perm=1000) 
  summary(im.perm) # threshold values
  summary(im.perm, alpha=0.10)
  }

}

\keyword{utilities}
