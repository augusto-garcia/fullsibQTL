 
#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: im_scan.R                                                     #
#                                                                     # 
# Contains:                                                           #
# im_scan                                                             #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# Updated by Rodrigo Amadeu                                           #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
# im_scan is based on scanone function from qtl package               #
#                                                                     #
# First version: 09/30/2011                                           #
# Last  version: 08/15/2017                                           #
# License: GPL-3                                                      #
#                                                                     #
#######################################################################


#######################################################################
#                                                                     #
# Function: im_scan                                                   #
#                                                                     #
# Performs QTL mapping on the genome, using Interval Mapping approach #
# This function is based on scanone from qtl pkg written by K. Broman #
#                                                                     #
# It does two basic analysis:                                         #
# if n.perm argument is missing or is 0, the QTL mapping is done.     #
# one can also include covariable in the model, as a fixed a fixed    #
# effect and using ordinary linear regression.                        #
# The result is an object from fullsib_scan class and also a matrix   #
# for xtable users.                                                   #
#                                                                     #
# If n.perm > 0, this function can perform permutation test           #
# The result is an object from fullsib_perm class.                    #
# EM algorithm is done in C code                                      #
#######################################################################

## -------------------------
## im_scan function

#' QTL searching using Interval Mapping approach
#' 
#' Performs QTL mapping with Interval Mapping (IM) approach by mixture models.
#' MLE are obtained using EM algorithm. Permutation test is also implemented
#' for threshold determination.
#' 
#' 
#' This function receives the object of class \emph{fullsib} and scan the
#' genome for QTL, using interval mapping. For each position, EM algorithm is
#' performed considering the following criteria for convergence:
#' 
#' \deqn{conv = abs\left[\frac{(new.lk - old.lk)}{old.lk}\right]}
#' 
#' in which, \eqn{old.lk} is the likelihood for \eqn{i^{th}} iteration and
#' \eqn{new.lk} is the likelihood for \eqn{(i+1)^{th}}.
#' 
#' In Chen and Storey's (Chen and Storey, 2006) article of having relaxed
#' threshold values, they showed different way of finding a QTL peak. However
#' to define a peak can be very trick and subjective, so for this function,
#' peak is the maximum value of LOD or \eqn{-\log_{10}(pvalue)} can assume in
#' each linkage group, for example, on fullsib_data (data set contained in this
#' package) has only four groups, so we consider that we can define 4 peaks. In
#' cases that dataset has more than 10 linkage groups, we restrict to storage
#' until the 2nd highest peaks.
#' 
#' The permutation test can be done for this \sQuote{n.perm} argument is higher
#' then 0, but if one wants to repeat the result of permutation, please used
#' \code{set.seed} function from before run \code{im_scan}.
#' 
#' Individuals with missing phenotypes are dropped.
#' 
#' @param fullsib An object from class \emph{fullsib}.
#' @param lg Vector indicating which linkage group will be scanned. Default is
#' to analyse all groups.
#' @param pheno.col Column number in the phenotype matrix (present in
#' \emph{fullsib} object) which should be used as the phenotype.
#' @param addcovar Additive covariates. If it is used, one must indicate the
#' design matrix for those source of variation. It should be noted that
#' additive covariates is included in the model as fixed effects, under
#' ordinary linear regression.
#' @param LOD if \code{TRUE} indicates the mapping result as LOD Score, if
#' \code{FALSE} QTL search is presents as \eqn{-\log_{10}(pvalue)}.
#' @param maxit Maximum number of iteration in EM algorithm.
#' @param tol Tolerance for determining convergence in EM algorithm. See
#' details.
#' @param n.perm If \code{n.perm} is missing or zero, usual genome scan is
#' performed, otherwise, permutation test is done. The number of permutation to
#' be used is defined for \code{n.perm} value
#' @param write.perm Optional string (e.g., \dQuote{file.txt}) to create a text
#' file containing the results obtained on the permutation test.
#' @param verbose If \code{TRUE} display information during EM algorithm. For
#' each position on the genome, it indicates the iteration of EM is performing,
#' with log-likelihood, convergence, genetics effects, square root of variance.
#' The log-likelihood of the model under \eqn{H_0} model is also showed.
#' @param verbose.scan If \code{TRUE} display a progress bar showing the
#' percentage of genome scan is done, based on the number of linkage groups to
#' be scanned. However, when permutation test is done, the option is coded as
#' FALSE for having a clear output.
#' 
#' @return This function can return two differents classes of objects
#' 
#' \dQuote{fullsib_scan} and \dQuote{fullsib_perm}:
#' 
#' If usual analysis is done (\sQuote{n.perm=0}), the function returns a matrix
#' with 4 columns: The first two columns indicates the linkage groups number
#' and position in centiMorgan the analysis has been done (\code{lg} and
#' \code{pos.cM}, respectivaly). The third is the value of QTL mapping test,
#' one can choose between \code{LOD} or \eqn{-\log_{10}(pvalue)}. The fourth
#' indicates which model was considered based on the level of informativeness
#' of the probabilities.  These values can vary between 0 to 9. For further
#' details about this see the \pkg{fullsibQTL} Vignette.  The names of each row
#' is important for QTL characterization, because it will be used as argument
#' for \code{cim_char} and \code{r2_ls} functions.
#' 
#' For situation that \sQuote{nperm} is higher then zero, permutation test is
#' performed and the function returns an object from \dQuote{fullsib_perm}.
#' This object is a matrix which the number of row indicates the number of
#' permutations and the number of column defines the number of highest peaks
#' colected on each permutation. As described on details, we define peak the
#' maximum value of either \code{LOD} or \eqn{-\log_{10}(pvalue)} of each
#' linkage group.  They are ordered from the maximum to minimum value. If one
#' determines the threshold with the distribuition of the highest peak (column
#' one), one have the threshold as present by Churchill and Doerge (1994). If
#' one choose the threshold with the second peak or lower is considering the
#' approach by Chen and Storey (2006). The function \emph{summary_fullsib_perm}
#' provides automatically these results.
#' @author Rodrigo Gazaffi, \email{rgazaffi@@gmail.com}
#' 
#' @seealso 
#' \code{\link[fullsibQTL]{create_fullsib}}
#' \code{\link[fullsibQTL]{im_char}}
#' \code{\link[fullsibQTL]{plot.fullsib_scan}}
#' \code{\link[fullsibQTL]{summary.fullsib_scan}}
#' \code{\link[fullsibQTL]{plot.fullsib_perm}}
#' \code{\link[fullsibQTL]{summary.fullsib_perm}}
#' 
#' @references 
#' Chen, L., Storey, J.D. (2006) Relaxed Significance Criteria for
#' Linkage Analysis. \emph{Genetics} 173: 2371-2381
#' 
#' Churchill, G.A., Doerge, R.W. (1994) Empirical Threshold Values for
#' Quantitative Trait Mapping. \emph{Genetics} 138: 963-971.
#' 
#' Gazaffi, R.; Margarido, G. R.; Pastina, M. M.; Mollinari, M.; Garcia, A. A.
#' F. (2014) A model for Quantitative Trait Loci Mapping, Linkage Phase, and
#' Segregation Pattern Estimation for a Full-Sib Progeny. \emph{Tree Genetics &
#' Genomes} 10(4): 791-801
#' 
#' @keywords utilities
#' 
#' @examples
#'   data(example_QTLfullsib)
#' 
#'   fsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ),
#'                              step = 0, map.function = "kosambi", condIndex = 3.5 )
#' 
#' 
#'   ## running im_scan:
#' 
#'   im1 <- im_scan( fsib, lg = "all", pheno.col = 1, LOD = TRUE )
#'   ## results in LOD Score for all linkage groups
#' 
#'   im2 <- im_scan( fsib, lg = 1:4, pheno.col = 1, LOD = FALSE )
#'   ## results in -log10(pvalue) for all linkage groups
#' 
#'   ## im_scan with using additive covariate
#'   covar <- matrix( rep( c( 1, -1 ), each = 150 ), ncol = 1 )
#'   im3 <- im_scan( fsib, lg = c( 2, 3 ), pheno.col = 2, addcovar = covar, LOD = FALSE )
#' 
#'   \dontrun{
#'   ## RESULTS:
#'   im1
#'   print( im1, lg = 3 ) ## printing only one linkage group
#' 
#'   ## detecting QTL peaks
#'   summary( im1 )
#'   summary( im1, thr = 6.5 )
#' 
#'   ##plotting
#'   plot( im1 )
#' 
#'   plot( im1, im2, col = c( "black", "blue" ), label.lg = c( "A", "B", "C", "D" ) )
#'   abline( 3, 0 ) ##threshold
#'   legend( "topleft",legend = c( "LOD", "-log10(pvalue)" ),
#'            col = c( "black", "blue" ), lwd = c( 2, 2 ), text.col = c( "black", "blue" ) )
#' 
#'   plot( im2, lg = 2 )
#'   plot( im2, lg = 2, incl.mkr = "points" )
#'   plot( im2, lg = 2, incl.mkr = "name" )
#'   }
#' 
#'   ## permutation test:
#'   \dontrun{
#'   im.perm <- im_scan( fsib, lg = c( 2, 4 ), pheno.col = 1, LOD = TRUE,
#'                       n.perm = 1000, write.perm = "perm_values.txt" )
#'   ## permutations done just using linkage groups 2 and 4
#'   summary( im.perm )
#'   plot( im.perm, peak = 1 )
#' }
#' 

im_scan <- function ( fullsib, lg, pheno.col = 1, addcovar = NULL, LOD = TRUE,
                      maxit = 1000, tol = 1e-04, n.perm, write.perm,
                      verbose, verbose.scan)
{

  ##check fullsib object
  if (!any(class(fullsib) == "fullsib")) 
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib'")

  ##checking lg argument
  if(missing(lg))
    lg <- "all"
  if(length(lg)==1 && lg == "all")
    lg <-  seq(1, length(fullsib$map), 1)
  else{
    check.lg <- match(lg, seq(1, length(fullsib$map), 1))
    lg.na <- which(is.na(check.lg))
    if(length(lg.na) > 0)
      stop("lg must be 'all' or any numbers between 1 to ",
           length(fullsib$map))
  }

  ##checking pheno.col
  if (is.na(match(pheno.col, 1:ncol(fullsib$pheno))) )
      stop("pheno.col must be a number between 1 and",
           ncol(fullsib$pheno))

  ##checking maxit
  if (maxit < 1) 
    stop("maxit argument must be >= 1")

  ##check tol
  if (tol <= 0) 
    stop("tol argument must be > 0")

  ##check verbose
  if (missing(verbose)) 
    verbose <- FALSE
  
  if (missing(verbose.scan)) 
    verbose.scan <- TRUE

  if (missing(write.perm))
    logical.perm <- FALSE
  else {
    logical.perm <- TRUE
    if(LOD == T)
      header.names <- matrix(c("lg", "pos.cM", "LOD", "model"),ncol=4)
    else
      header.names <- matrix(c("lg", "pos.cM", "-log10(pval)", "model"),ncol=4)
    write.table(header.names, file=write.perm, sep=";",
                col.names=FALSE, row.names=FALSE, dec=".")
  }
  
  if (!missing(n.perm) && n.perm > 0) {
    m.peak <- min(c(2, length(fullsib$probs),length(lg)))
    results <- matrix(NA, ncol = m.peak, nrow = n.perm)

    cat("Permutation test: This operation may take a VERY long time\n")
    pb <- txtProgressBar(style=3)
    setTxtProgressBar(pb, 0)

    for (i in 1:n.perm) {

      pheno.idx <- which(!is.na(fullsib$pheno[,pheno.col]))
      fullsib$pheno[pheno.idx, pheno.col] <- fullsib$pheno[sample(pheno.idx),pheno.col]

      if(verbose == TRUE)  cat("\nPermutation: ", i, "\n")

      temp <- im_scan(fullsib, lg=lg, pheno.col=pheno.col, addcovar=addcovar, maxit=maxit, tol=tol, LOD=LOD, verbose=verbose,
                      n.perm=0, verbose.scan=FALSE)

     ##if NaN is returned substitute for 0 avoinding the program crashes
      if(length(is.na(temp[,3])) > 0)
        temp[which(is.na(temp[,3])),3] <- 0

      if(logical.perm == TRUE){
        write.table(temp, file=write.perm, sep=";", dec=".", append=T, col.names=FALSE)
      }
      
      results[i, ] <- rev(sort(tapply(temp[,3], temp[,"lg"], max)))[1:m.peak]

      setTxtProgressBar(pb, i/n.perm)
    }
    close(pb)
    colnames(results) <- paste("peak", 1:ncol(results), sep=".")
    class(results) <- c("fullsib_perm", "matrix")
    return(results)
  }


  if (is.null(addcovar) == TRUE) {
    n.ac <- 0
    ac <- NULL
  }
  else {
    ac <- as.matrix(addcovar)
    n.ac <- ncol(ac)
  }
  
  ##keep individuals with phenotype (not NA)
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]
  if(n.ac > 0) ac <- ac[pheno.index,]

  ##preliminary information to send for C...
  if(n.ac > 0)
    lm1 <- lm(pheno ~ ac)
  else
    lm1 <- lm(pheno ~ 1)
  ac <- model.matrix(lm1)
  n.ac <- ncol(ac)
  Zls <- t(solve(crossprod(ac), t(ac)))
  resid0 <- lm1$resid
  gamma <- coef(lm1)
  
  ##results will be in scan.results
  scan.results <- NULL

  if(verbose.scan == TRUE){   
    cat("QTL mapping for ", length(lg)," groups\n")
    count.lg <- 1
    pb <- txtProgressBar(style=3)
    setTxtProgressBar(pb, 0)
  }
  for (LGs in lg) { ##loop for linkage groups
    dim_probs <- dim(fullsib$probs[[LGs]]$cond.prob)
    if (is.null(dim_probs) == TRUE)
      next
    ##just a check, group with 2 markers with rf=0.0 falls here...
    else {
      if(verbose == TRUE) cat("\n\nPositions to scan: ", length(fullsib$probs[[LGs]]$newmap$dist),"\n")
      n.ind <- get("dim_probs")[1]
      n.pos <- get("dim_probs")[2]
      n.gen <- get("dim_probs")[3]

      z <- .C("R_scan_qtl",
              as.integer(length(pheno.index)),
              as.integer(n.pos), 
              as.integer(fullsib$probs[[LGs]]$colin$geno.class),
              as.integer(fullsib$probs[[LGs]]$colin$type), 
              as.double(fullsib$probs[[LGs]]$cond.prob[pheno.index, ,]), 
              as.double(ac),
              as.integer(n.ac),
              as.double(gamma),
              as.double(Zls),
              as.double(pheno),
              result = as.double(rep(0, n.pos)),#ln_likelihood
              as.integer(maxit), 
              as.double(tol),
              #as.integer(verbose))$result
              as.integer(verbose), PACKAGE = "fullsibQTL")$result
      if(LOD == TRUE)
        z <- z / log(10)
      else
        z <- -log10(pchisq(2*(z),
                           fullsib$probs[[LGs]]$colin$geno.class - 1,
                           lower.tail = FALSE))
      pos.cM <- fullsib$probs[[LGs]]$newmap$dist
      scan.results <- rbind(scan.results,
                            cbind(lg = rep(LGs, length(pos.cM)),
                                  pos.cM = pos.cM, z,
                                  fullsib$probs[[LGs]]$colin$type))
    }
    if(verbose.scan == TRUE){
      setTxtProgressBar(pb, count.lg/length(lg))
      count.lg <- count.lg + 1
    }
  }
  if(verbose.scan == TRUE) close(pb)
  
  if(LOD == TRUE)
    colnames(scan.results) <- c("lg", "pos.cM", "LOD", "model")
  else
    colnames(scan.results) <- c("lg", "pos.cM","-log10(pval)", "model")
  class(scan.results) <- c("fullsib_scan", "matrix")
  scan.results
}
