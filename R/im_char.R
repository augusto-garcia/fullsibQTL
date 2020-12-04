 #######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: im_char.R                                                     #
#                                                                     # 
# Contains:                                                           #
# im_char                                                             #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# Updated by Rodrigo Amadeu                                           #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
#                                                                     #
# First version: 09/30/2011                                           #
# Last  version: 08/15/2017                                           #
# License: GPL-3                                                      #
#                                                                     #
#######################################################################


#######################################################################
#                                                                     #
# Function: im_char                                                   #
#                                                                     #
# This function receives the position that qtl caracterization is done#
# i.e., one estimate the genetic effects of QTL and test if they are  #
# significative, and segregation. Linkage phase is obtained by the    #
# intepretation of QTL effects signals and better view with function  #
# draw_phase that used the result of this function...                 #
# all the MLE and tests are done in C code                            #
#######################################################################

## -------------------------
## im_char function

#' QTL characterization using Interval Mapping
#' 
#' Considering a detected QTL on a given position, this function allows to
#' estimate the genetic effects and also provides a series os estatistical test
#' allowing to infer QTL segregation patterns and its linkage phases with
#' flanking markers.
#' 
#' The method implemented in this package consists first to scan the genome for
#' a QTL (done with \code{im_scan}) and than on a second step, the QTL is
#' characterized (done with \code{im_char}).
#' 
#' In our model, QTL can segregate in one of four possible types: 1:1:1:1,
#' 1:2:1, 3:1 or 1:1. The QTL segregation may vary as function of the number of
#' significative effects and their magnitude between each other. To understand
#' with details, check Gazaffi et al. (2014) and \pkg{fullsibQTL} Vignette. To
#' help the user with the understanding of mapped QTL, the function
#' \code{get_segr} and \code{draw_phase} were devolped to infer directly the
#' segregation pattern and the linkage phase between QTL and markers,
#' respectively.
#' 
#' @param fullsib An object from class \emph{fullsib}.
#' @param pheno.col Column number in the phenotype matrix (present in
#' \emph{fullsib} object) which should be used as the phenotype.
#' @param addcovar Additive covariates. If it is used, one must indicate the
#' design matrix for those source of variation. It should be noted that
#' additive covariates is included in the model as fixed effects, under
#' ordinary linear regression.
#' @param lg Indicates which linkage group contains the region to be
#' characterized.
#' @param pos String representing the name of the locus that will be studied.
#' This name is found as the row label for the matrix returned by
#' \code{cim_scan}.
#' @param maxit Maximum number of iteration in EM algorithm.
#' @param tol Tolerance for determining convergence in EM algorithm. See
#' details.
#' @param verbose If \code{TRUE} display information during EM algorithm. It
#' indicates the iteration of EM is performing, with log-likelihood,
#' convergence, genetics effects, square root of variance. The log-likelihood
#' of the model under \eqn{H_0} model is also showed.
#' 
#' @return An object of class \emph{fullsib_char} is returned, consisting in a
#' matrix of one column and 15 rows. The information provided are for each row
#' are: linkage group that was analysed, position (cM),
#' \eqn{-\log_{10}(pvalue)}, LOD Score, model intercept (\code{mu}), QTL effect
#' for parent \eqn{P} (\code{alpha_p}) and its LOD Score (\code{LOD_H1}), QTL
#' effect for parent \eqn{Q} (\code{alpha_q}) and its LOD Score
#' (\code{LOD_H2}), dominance effect (\code{delta_pq} and its LOD Score
#' (\code{LOD_H3}). Additional tests are also performed to identify the
#' segregation pattern (\code{LOD_H4}, \code{LOD_H5} and \code{LOD_H6}). For
#' getting the segregation estimated \code{get_segr} function should be used.
#' The last row indicated the \code{model} used in the analysis. See
#' \pkg{fullsibQTL} Vignette for details.
#' 
#' @author Rodrigo Gazaffi, \email{rgazaffi@@gmail.com}
#' 
#' @seealso 
#' \code{\link[fullsibQTL]{im_scan}}
#' \code{\link[fullsibQTL]{cim_char}}
#' \code{\link[fullsibQTL]{draw_phase}}
#' \code{\link[fullsibQTL]{get_segr}}
#' 
#' @references
#' Gazaffi, R.; Margarido, G. R.; Pastina, M. M.; Mollinari, M.; Garcia, A. A.
#' F. (2014) A model for Quantitative Trait Loci Mapping, Linkage Phase, and
#' Segregation Pattern Estimation for a Full-Sib Progeny. \emph{Tree Genetics &
#' Genomes} 10(4): 791-801
#' 
#' @keywords utilities
#' 
#' @examples
#'   data( "example_QTLfullsib" )
#' 
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ),
#'                              step = 3, map.function = "kosambi", condIndex = 3.5)
#' 
#'   \dontrun{
#'   im1 <- im_scan( fullsib, pheno.col = 1 )
#'   summary( im1 )
#'   }
#' 
#'   qtl1 <- im_char( fullsib, pheno.col = 1, lg = 1, pos = "loc3" )
#'   qtl2 <- im_char( fullsib, pheno.col = 1, lg = 2, pos = "M27" )
#'   qtl3 <- im_char( fullsib, pheno.col = 1, lg = 3, pos = "loc63" )
#'   qtl4 <- im_char( fullsib, pheno.col = 1, lg = 4, pos = "loc57" )
#' 
#'   \dontrun{
#'   covar <- matrix( rep( c( 1, -1 ), each = 150 ), ncol = 1 )
#'   im2 <- im_scan( fullsib, pheno.col = 2, addcovar = covar )
#'   summary( im2 )
#' 
#'   qtl5 <- im_char( fullsib, pheno.col = 1, addcovar = covar, lg = 2, pos = "M27" )
#'   }
#' 
#' @importFrom stats lm model.matrix pchisq
#' @export
im_char <- function( fullsib, pheno.col = 1, addcovar = NULL, lg, pos,
                     maxit = 1000, tol = 1e-10, verbose = FALSE ){
  
  ##checking input arguments
  if (!any(class(fullsib) == "fullsib"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib'")
  
  if (is.na(match(pheno.col, 1:ncol(fullsib$pheno))) )
    stop("pheno.col must be a number between 1 and",ncol(fullsib$pheno))
  
  if(missing(lg))
    stop("lg argument is missing")
  if(!(lg %in% seq(1,length(fullsib$map),1)))
    stop("lg must be an integer between 1 and ", length(fullsib$map))
  
  if(missing(pos))
    stop("pos argument is missing")
  ##it can be a number from 1 to last position of the group, good to loops...
  if(is.numeric(pos) && !(pos %in% seq(1,length(fullsib$probs[[lg]]$newmap$dist),1)))
    stop("pos argument may be an integer between 1 and ",
         length(fullsib$probs[[lg]]$newmap$dist))
  
  if(is.character(pos)){
    pos.arg <- match.arg(pos, names(fullsib$probs[[lg]]$newmap$dist))
    if(!(pos.arg %in% names(fullsib$probs[[lg]]$newmap$dist)))     
      stop("pos argument may be an string that contain the name of the loci: please check it")
    pos <- which(names(fullsib$probs[[lg]]$newmap$dist) == pos.arg)
  }


  if(is.null(addcovar)==TRUE){
    n.ac <- 0
    ac <- NULL
  }
  else{
    ac <- as.matrix(addcovar)
    n.ac <- ncol(ac)
  }
  
  if (maxit < 1) 
    stop("maxit must be >= 1")
  
  if (tol <= 0) 
    stop("tol must be > 0")
  
  
  ##keeping non NA elements
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]
  
  if(n.ac > 0) ac <- ac[pheno.index,]
  
  if(n.ac > 0)
    lm1 <- lm(pheno ~ ac)
  else
    lm1 <- lm(pheno ~ 1)

  ac <- model.matrix(lm1)
  n.ac <- ncol(ac)
  Zls <- t(solve(crossprod(ac), t(ac)))
  resid0 <- lm1$resid
  gamma <- coef(lm1)

  ##getting the flanking markers or closest markers
  name.pos <- names(fullsib$probs[[lg]]$newmap$dist[pos])
  if(substring(name.pos, 1, 3) == "loc"){
    mkrs <- which(substring(names(fullsib$probs[[lg]]$newmap$dist),1,3) != "loc")
    flank <- which(cumsum(pos < mkrs) == 1)
    mk.rgt <- names(fullsib$probs[[lg]]$newmap$dist)[ mkrs[flank] ]
    mk.lft <- names(fullsib$probs[[lg]]$newmap$dist)[ mkrs[flank - 1] ]
    closest.mkr <-  paste(mk.lft,mk.rgt, sep="-")
  }
  else
    closest.mkr <- name.pos

  z <- .C("R_char_qtl",
          as.integer(length(pheno.index)),
          as.integer(fullsib$probs[[lg]]$colin$geno.class[pos]),
          as.integer(fullsib$probs[[lg]]$colin$type[pos]),      
          as.double(fullsib$probs[[lg]]$cond.prob[pheno.index,pos,]),
          as.double(ac),
          as.integer(n.ac),
          as.double(gamma),
          as.double(Zls),
          as.double(pheno),        
          result=as.double(rep(0, 12)), #this is LOD
          as.integer(maxit),
          as.double(tol),
          as.integer(verbose), PACKAGE="fullsibQTL")$result
  ##as.integer(verbose))$result 
  z <- as.matrix(z)

  pval <- -log10(pchisq(z[1]*4.61,
                        fullsib$probs[[lg]]$colin$geno.class[pos] - 1,
                        lower.tail = FALSE)) 

  ##convert LOD to pvalue: H04, H05 and H06
  z[9:11] <- pchisq(z[9:11]*4.61, 1, lower.tail = FALSE)
  z <- rbind(lg, fullsib$probs[[lg]]$newmap$dist[pos], pval, z)
 
  rownames(z) <- c("LG", "pos", "-log10(pval)","LOD_Ha", "mu",
                   "alpha_p", "LOD_H1","alpha_q","LOD_H2", "delta_pq",
                   "LOD_H3", "H4_pvalue", "H5_pvalue", "H6_pvalue", "model")
  
  colnames(z) <- closest.mkr
  structure(z, class = c("fullsib_char", "matrix"))
}#end function
