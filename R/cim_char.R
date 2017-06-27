#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: cim_char.R                                                    #
#                                                                     # 
# Contains: cim_char                                                  #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
#                                                                     #
# First version: 2011/30/09                                           #
# Last  version: 2011/30/09                                           #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: cim_char                                                  #
#                                                                     #
# This function receives the position that qtl caracterization is done#
# i.e., one estimate the genetic effects of QTL and test if they are  #
# significative, and segregation. Linkage phase is obtained by the    #
# intepretation of QTL effects signals and better view with function  #
# draw_phase that used the result of this function.                   #
# all the MLE and tests are done in C code                            #
#######################################################################

## -------------------------
## cim_char function

#' QTL characterization using Composite Interval Mapping
#' 
#' Considering a detected QTL on a given position, this function allows to
#' estimate the genetic effects and also provides a series os estatistical test
#' allowing to infer QTL segregation patterns and its linkage phases with
#' flanking markers.
#' 
#' The method implemented in this package consists first to scan the genome for
#' a QTL (done with \code{cim_scan}) and than on a second step, the QTL is
#' characterized (done with \code{cim_char}).
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
#' For convergence during EM iteration we use the following criterion:
#' \deqn{conv = abs\left[\frac{(new.lk - old.lk)}{old.lk}\right]} in which,
#' \eqn{old.lk} is the likelihood for \eqn{i^{th}} iteration and \eqn{new.lk}
#' is the likelihood for \eqn{(i+1)^{th}}. If convergence value was higher than
#' \code{tol} argument, iterative process continues, if not iteration is
#' stopped.
#' 
#' If \code{icim} is TRUE, first a linear regression is done to remove the
#' cofactors, intercept and additive covariates (if defined) effects from the
#' phenotype. The residual of this analysis is used as a \sQuote{new}
#' phenotype. The QTL analysis is done similar an interval mapping.
#' 
#' @param fullsib An object from class \emph{fullsib_cofactors}.
#' @param pheno.col Column number in the phenotype matrix (present in
#' \emph{fullsib_cofactors} object) which should be used as the phenotype.
#' @param ws Window Size in cM. Default is 10 cM (i.e., 5 cM for each side.
#' @param lg Integer indicating which linkage group will be studied.
#' @param pos String representing the name of the locus that will be studied.
#' This name is found as the row label for the matrix returned by
#' \code{cim_scan}.
#' @param maxit Maximum number of iteration in EM algorithm.
#' @param tol Tolerance for determining convergence in EM algorithm. See
#' details.
#' @param icim if \code{TRUE}, icim approach proposed by Li et al. (2007) is
#' done (extended for full sib progeny), if \code{FALSE} (default) traditional
#' CIM is performed. See details for icim method.
#' @param verbose If \code{TRUE} display information during EM algorithm. It
#' indicates the iteration of EM is performing, with log-likelihood,
#' convergence, genetics effects, square root of variance. The log-likelihood
#' of the model under \eqn{H_0} model is also showed
#' @return An object of class \emph{fullsib_char} is returned, consisting in a
#' matrix of one column and 15 rows. The information provided are for each row
#' are: linkage group that was analysed, position (cM),
#' \eqn{-\log_{10}(pvalue)}, LOD Score, model intercept (\code{mu}), QTL effect
#' for parent \eqn{P} (\code{alpha_p}) and its LOD Score (\code{LOD_H1}), QTL
#' effect for parent \eqn{Q} (\code{alpha_q}) and its LOD Score
#' (\code{LOD_H2}), dominance effect (\code{delta_pq} and its LOD Score
#' (\code{LOD_H3}). Additional tests are also performed to identify the
#' segregation pattern (\code{LOD_H4}, \code{LOD_H5} and \code{LOD_H6}).
#' 
#' For getting the segregation estimated \code{get_segr} function should be
#' used. The last row indicated the \code{model} used in the analysis. See
#' \pkg{fullsibQTL} Vignette for details.
#' @author Rodrigo Gazaffi, \email{rgazaffi@@gmail.com}
#' @seealso 
#' \code{\link[fullsibQTL]{cim_scan}}
#' \code{\link[fullsibQTL]{im_char}}
#' \code{\link[fullsibQTL]{draw_phase}}
#' \code{\link[fullsibQTL]{get_segr}}
#' 
#' @references
#' 
#' Gazaffi, R.; Margarido, G. R.; Pastina, M. M.; Mollinari, M.; Garcia, A. A.
#' F. (2014) A model for quantitative trait loci mapping, linkage phase, and
#' segregation pattern estimation for a full-sib progeny. \emph{Tree Genetics &
#' Genomes} 10(4): 791-801
#' 
#' Li, H., Ye, G.; Wang, J. (2007) A Modified Algorithm for the Improvement of
#' Composite Interval Mapping. \emph{Genetics} 175: 361-374
#' @keywords utilities
#' @examples
#'   data(example_QTLfullsib)
#' 
#'   fullsib <- create_fullsib(example_QTLfullsib,
#'                             list(LG1_final, LG2_final, LG3_final, LG4_final),
#'                             step=0,map.function="kosambi",condIndex=3.5)
#' 
#' 
#'   ###############################################
#'   ## cofactor selection using BIC (n.ind = 300)
#'   cofs.fs <- cof_selection(fullsib, pheno.col=1, k = log(300),
#'                            selection=1) 
#' 
#'   \dontrun{
#'   cim1 <- cim_scan(cofs.fs, pheno.col=1, ws = 22, LOD= TRUE, icim=FALSE)
#'   summary(cim1)
#'   }
#' 
#'   qtl <- cim_char(cofs.fs, pheno.col=1, ws=22, lg=3, pos="M38")
#' 

cim_char <- function(fullsib, pheno.col=1, ws = 10, lg, pos,
                     maxit=1000,tol=1e-10, icim=FALSE, verbose=FALSE)
{
  ##checking input arguments
  if (!any(class(fullsib) == "fullsib_cofactors"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib_cofactors': execute 'cof.selection' or 'cof_definition' first")
  
  if(ws < 0)
    stop("window size must be non negative number")
  ws <- ws / 2 #similar to qtl
  
  if (is.na(match(pheno.col, 1:ncol(fullsib$pheno))) )
    stop("pheno.col must be a number between 1 and", ncol(fullsib$pheno))
  
  if(missing(lg))
    stop("lg argument is missing")
  if(!(lg %in% seq(1,length(fullsib$map),1)))
    stop("lg must be an integer between 1 and ", length(fullsib$map))
  
  if(missing(pos))
    stop("pos argument is missing")
  if(is.character(pos)){
    pos.arg <- match.arg(pos, names(fullsib$probs[[lg]]$newmap$dist))
    if(!(pos.arg %in% names(fullsib$probs[[lg]]$newmap$dist)))  
      stop("pos argument may be an string that contain the name of the loci: please check this name")
    pos <- which( names(fullsib$probs[[lg]]$newmap$dist) == pos.arg)
  } 
  
  if(pheno.col != fullsib$cofactors$trait.cof)
    stop("pheno.col argument is different than pheno.col used to obtain cofactors")
  
  ##keeping non NA elements
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]
  if(icim == TRUE)
    pheno.icim <- pheno
  
  ##obtaining the model that will be used for that position
  cim.models <- model_select_lg(fullsib, ws = ws, LG=lg)

  for(mods in 1:length(cim.models)){
    find.pos <- which(pos %in% cim.models[[mods]]$position)
    if(length(find.pos) > 0)
      flag.mkr <- mods
  }

  if(cim.models[[flag.mkr]]$mkr2rm[1] == "none")
    cim.cof.matrix <- fullsib$cofactors$matrix.cof[pheno.index,]
  else{
    mkr2rm <-do.call("c",lapply(cim.models[[flag.mkr]]$mkr2rm,
           function(x) which(colnames(fullsib$cofactors$matrix.cof) == x)))

    cim.cof.matrix <- fullsib$cofactors$matrix.cof[pheno.index,-mkr2rm]
    
  }

  
  unar <- matrix(1,nrow(cim.cof.matrix))
  colnames(unar) <- "intercept"
  cim.cof.matrix <- cbind(unar,cim.cof.matrix)
  
  if(icim == FALSE){
    ##traditional analysis
    Zls <- t(solve(crossprod(cim.cof.matrix), t(cim.cof.matrix)))
    gamma <- coef(lm(pheno ~ cim.cof.matrix -1))
  }
  else{
    ##icim proposed -  remove cofactors from phenotype
    pheno.rm <- lm(pheno.icim ~ cim.cof.matrix - 1)
    pheno <- pheno.rm$residuals
    gamma <- coef(pheno.rm)[1]#for intercept
    Zls <- t(solve(crossprod(unar), t(unar)))
    cim.cof.matrix <- unar
  }
  
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
  
  dim_probs <- dim(fullsib$probs[[lg]]$cond.prob)
  z <- .C("R_char_qtl",
          as.integer(length(pheno.index)),        
          as.integer(fullsib$probs[[lg]]$colin$geno.class[pos]),
          as.integer(fullsib$probs[[lg]]$colin$type[pos]),      
          as.double(fullsib$probs[[lg]]$cond.prob[pheno.index,pos,]),
          as.double(cim.cof.matrix),
          as.integer(ncol(cim.cof.matrix)),
          as.double(gamma),#
          as.double(Zls), #
          as.double(pheno),       
          result=as.double(rep(0, 12)),
          as.integer(maxit),
          as.double(tol),
          as.integer(verbose),PACKAGE="fullsibQTL")$result
          ##as.integer(verbose))$result 
  z <- as.matrix(z)


  pval <-  -log10(pchisq(z[1]*4.61,  fullsib$probs[[lg]]$colin$geno.class[pos] - 1, lower.tail = FALSE))
  #print(z[9:11])
  z[9:11] <- pchisq(z[9:11]*4.61, 1, lower.tail = FALSE)
  #print(z[9:11])
  
  z <- rbind(lg, fullsib$probs[[lg]]$newmap$dist[pos], pval, z)

  rownames(z) <- c("LG", "pos", "-log10(pval)", "LOD_Ha", "mu",
                   "alpha_p", "LOD_H1", "alpha_q", "LOD_H2", "delta_pq",
                   "LOD_H3", "H4_pvalue", "H5_pvalue", "H6_pvalue", "model")
  colnames(z) <- closest.mkr

  structure(z, class = c("fullsib_char", "matrix"))
}