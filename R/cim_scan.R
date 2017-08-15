#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: cim_scan.R                                                    #
#                                                                     # 
# Contains: cim_scan                                                  #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# Updated by Rodrigo Amadeu                                           #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
#                                                                     #
# First version: 09/30/2011 (american date format)                    #
# Last  version: 08/15/2017 (american date format)                    #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: cim_scan                                                  #
#                                                                     #
# Performs QTL mapping on the genome using Composite Interval Mapping #
# This function is based on scanone from qtl pkg written by K. Broman #
#                                                                     #
# It does two basic analysis:                                         #
# if n.perm argument is missing or is 0, the QTL mapping is done.     #
# The result is an object from fullsib_scan class and also a matrix   #
# for xtable users.                                                   #
#                                                                     #
# If n.perm > 0, this function can perform permutation test           #
# The result is an object from fullsib_perm class.                    #
# EM algorithm is done in C code                                      #
#######################################################################

## -------------------------
## cim_scan function

#' QTL searching using Composite Interval Mapping
#' 
#' Performs QTL mapping using Composite Interval Mapping (CIM) approach, with
#' mixture models. MLE are obtained using EM algorithm. Permutation test is
#' also implemented for threshold determination
#' 
#' The mapping procedure is performed similar as present in \code{im_scan}.
#' However the inclusion of cofactors can locate QTL with more precision with
#' higher statistical power in comparison of interval mapping approach.
#' 
#' The usage of the cofactors on the CIM is function of window size. Here, we
#' consider window size minimal distance that cofactors should be from the
#' interval under analysis. For example, if one considers a window size of 20
#' cM, while CIM is being performed, any marker located closer than 10 cM for
#' each side, will not be included in the analysis.
#' 
#' If extra additive covariates is used, one need to specify them, on
#' \code{cof_selection} or \code{cof_definition} step.
#' 
#' Li et al. (2007) proposed an adaption for CIM first developed by Zeng
#' (1994). In their metodology, first the phenotype should be regressed to
#' remove the variation due to cofactors and intercept. On a second moment, the
#' residual is used to map QTL, similar an interval mapping approach. The
#' method was proposed for inbred lines, and we included this method on CIM
#' done in fullsib context.
#' 
#' For convergence during EM iteration we use the following criterion:
#' \deqn{conv = abs\left[\frac{(new.lk - old.lk)}{old.lk}\right]} in which,
#' \eqn{old.lk} is the likelihood for \eqn{i^{th}} iteration and \eqn{new.lk}
#' is the likelihood for \eqn{(i+1)^{th}}. If convergence value was higher than
#' \code{tol} argument, iterative process continues, if not iteration is
#' stopped.
#' 
#' We also implemented permutation test on \code{cim_scan}, including the
#' approach proposed by Chen and Storey (Chen and Storey, 2006). In this case,
#' to define peak we used the maximum value of LOD or \eqn{-\log_{10}(pvalue)}
#' in each linkage group. In fullsib_data (example dataset) there are four
#' linkage groups, so we consider four peaks, if a dataset has more than 10
#' linkage groups, we restrict to storage until the 10th highest peaks.
#' 
#' Individuals with missing phenotypes are dropped.
#' 
#' @param fullsib An object from class \emph{fullsib_cofactors}.
#' @param lg Vector indicating which linkage group will be scanned. Default is
#' to analyse all groups.
#' @param pheno.col Column number in the phenotype matrix (present in
#' \emph{fullsib_cofactors} object) which should be used as the phenotype.
#' @param ws Window Size in cM. Default is 10 cM (i.e., 5 cM for each side.
#' @param LOD if \code{TRUE} indicates the mapping result as LOD Score, if
#' \code{FALSE} QTL search is presents as \eqn{-\log_{10}(pvalue)}.
#' @param maxit Maximum number of iteration in EM algorithm.
#' @param tol Tolerance for determining convergence in EM algorithm. See
#' details.
#' @param n.perm If \code{n.perm} is missing or zero, usual genome scan is
#' performed, otherwise, permutation test is done. The number of permutation to
#' be used is defined for \code{n.perm} value.
#' @param write.perm Optional string (e.g., \dQuote{file.txt}) to create a text
#' file containing the results obtained on the permutation test.
#' @param icim if \code{TRUE}, icim approach proposed by Li et al. (2007) is
#' done (extended for full sib progeny), if \code{FALSE} (default) traditional
#' CIM is performed. See details for icim method.
#' @param verbose If \code{TRUE} display information during EM algorithm. For
#' each position on the genome, it indicates the iteration of EM is performing,
#' with log-likelihood, convergence, genetics effects, square root of variance.
#' The log-likelihood of the model under \eqn{H_0} model is also showed.
#' @param verbose.scan If \code{TRUE} display a progress bar showing the
#' percentage of genome scan is done, based on the number of linkage groups to
#' be scanned. However, when permutation test is done, the option is coded as
#' FALSE for having a clear output.
#' @return
#' 
#' This function can return two differents classes of objects
#' \dQuote{fullsib_scan} and \dQuote{fullsib_perm}:
#' 
#' If usual analysis is done (\sQuote{n.perm=0}), the function returns a matrix
#' with 4 columns: The first two columns indicates the linkage groups number
#' and position in centiMorgan the analysis has been done (\code{lg} and
#' \code{pos.cM}, respectivaly). The third is the value of QTL mapping test,
#' one can choose between \code{LOD} or \eqn{-\log_{10}(pvalue)}. The fourth
#' indicates which model was considered based on the level of informativeness
#' of the probabilities.  These values can vary between 0 to 9. For further
#' details about this see the \pkg{fullsibQTL} Vignette. The names of each row
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
#' approach by Chen and Storey (2006). The function \code{summary.fullsib_perm}
#' or just \code{summary} provides automatically these results.
#' @author Rodrigo Gazaffi, \email{rgazaffi@@gmail.com}
#' @seealso \code{\link[fullsibQTL]{create_fullsib}},
#' 
#' \code{\link[fullsibQTL]{im_scan}},
#' 
#' \code{\link[fullsibQTL]{cim_char}},
#' 
#' \code{\link[fullsibQTL]{cof_selection}},
#' 
#' \code{\link[fullsibQTL]{cof_definition}}
#' @references Chen, L., Storey, J.D. (2006) Relaxed Significance Criteria for
#' Linkage Analysis. \emph{Genetics} 173: 2371-2381
#' 
#' Gazaffi, R.; Margarido, G. R.; Pastina, M. M.; Mollinari, M.; Garcia, A. A.
#' F. (2014) A model for Quantitative Trait Loci Mapping, Linkage Phase, and
#' Segregation Pattern Estimation for a Full-Sib Progeny. \emph{Tree Genetics &
#' Genomes} 10(4): 791-801
#' 
#' Li, H., Ye, G.; Wang, J. (2007) A Modified Algorithm for the Improvement of
#' Composite Interval Mapping. \emph{Genetics} 175: 361-374
#' 
#' Zeng, Z. B. (1994) Precision Mapping of Quantitative Trait Loci.  \emph{Proc
#' Natl Acad Sci U S A} 90: 10,972-10,976
#' @keywords utilities
#' @examples
#' 
#' 
#'   data( "example_QTLfullsib" )
#' 
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ),
#'                              step = 0 , map.function = "kosambi", condIndex = 3.5 )
#' 
#' 
#'   ###############################################
#'   ## cofactor selection using BIC (n.ind = 300)
#'   cofs.fs <- cof_selection( fullsib, pheno.col = 1, k = log( 300 ),
#'                             selection = 1 ) 
#' 
#'   cim1 <- cim_scan( cofs.fs, pheno.col = 1, ws = 22, LOD = TRUE, icim = FALSE )
#' 
#'   \dontrun{
#'   covar <- matrix( rep( c( 1, -1 ), each = 150 ), ncol = 1 )
#'   cofs.fs <- cof_selection( fullsib, pheno.col = 2, addcovar = covar, k = log(300) )
#'   cim2 <- cim_scan( cofs.fs, pheno.col = 2, ws = 22, LOD = TRUE, icim = TRUE ) 
#'   ## cim with covariate, one just indicate the additive covar on
#'   ## cof_selection or cof_definition
#' 
#'   ##permutation test:
#'   cim.perm <- cim_scan( cofs.fs, pheno.col = 1, ws = 22, LOD = FALSE,
#'                         n.perm=1000) 
#'   summary( cim.perm ) # threshold values
#'   summary( cim.perm, alpha = 0.10 )
#'   }
#' 

cim_scan <- function(fullsib, lg, pheno.col=1, ws = 10, LOD=TRUE,
                     maxit=1000, tol=1e-04, n.perm=0, write.perm,
                     icim=FALSE, verbose, verbose.scan)
{

  ##checking arguments
  if (!any(class(fullsib) == "fullsib_cofactors"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib_cofactors': execute 'cof_selection' or 'cof_definition' first")

  if(missing(lg))
    lg <- "all"
  if(length(lg)==1 && lg == "all")
    lg <-  seq(1, length(fullsib$map), 1)
  else{
    check.lg <- match(lg, seq(1, length(fullsib$map), 1))
    lg.na <- which(is.na(check.lg))
    if(length(lg.na) > 0)
      stop("lg must be 'all' or any numbers between 1 to ", length(fullsib$map))
  }
  
  if (is.na(match(pheno.col, 1:ncol(fullsib$pheno))) )
    stop("pheno.col must be a number between 1 and", ncol(fullsib$pheno))
    
  if(ws < 0)
    stop("window size must be non negative number")
  ##similar to Rqtl
  ws <- ws / 2 

  if (maxit < 1) 
    stop("maxit argument must be >= 1")
  
  if (tol <= 0) 
    stop("tol argument must be > 0")

  if (missing(verbose)) 
    verbose <- FALSE

  if (missing(verbose.scan)) 
    verbose.scan <- TRUE

  if (missing(write.perm))
    logical.perm <- FALSE
  else{
    logical.perm <- TRUE
    if(LOD == T)
      header.names <- matrix(c("lg", "pos.cM", "LOD", "model"),ncol=4)
    else
      header.names <- matrix(c("lg", "pos.cM", "-log10(pval)", "model"),ncol=4)
    if(!is.character(write.perm))
      write.perm <- as.character(write.perm) 
    write.table(header.names, file=write.perm, sep=";",
                col.names=FALSE, row.names=FALSE, dec=".")
  }
  if(is.logical(icim) == FALSE)
    stop("icim argument must be FALSE or TRUE")
  
  if(pheno.col != fullsib$cofactors$trait.cof)
    stop("pheno.col argument is different than pheno.col used to obtain cofactors")
  

  ##permutation
  if (!missing(n.perm) && n.perm > 0) {
    ##results <- matrix(ncol = 1, nrow = n.perm)
    m.peak <- min(c(10, length(fullsib$probs), length(lg)))
    results <- matrix(NA, ncol = m.peak, nrow = n.perm)

    cat("Permutation test: This operation may take a VERY long time\n")
    pb <- txtProgressBar(style=3)
    setTxtProgressBar(pb, 0)

    for (i in 1:n.perm) {
      ##necessary to make sure the same ind will be used in permutation
      pheno.idx <- which(!is.na(fullsib$pheno[,pheno.col]))
      fullsib$pheno[pheno.idx, pheno.col] <- fullsib$pheno[sample(pheno.idx),pheno.col]

      if(verbose == TRUE)  cat("\nPermutation: ", i, "\n")

      temp <- cim_scan(fullsib, lg=lg, pheno.col=pheno.col, ws=ws, maxit=maxit, tol=tol,
                       LOD=LOD, n.perm=0, icim=icim, verbose=verbose, verbose.scan=FALSE)
      
      if(logical.perm == TRUE){
        write.table(temp, file=write.perm, sep=";", dec=".", append=T, col.names=FALSE)
      }

      ##if NaN is returned substitute for 0, so avoids the program crashes
      if(length(is.na(temp[,3])) > 0)
        temp[which(is.na(temp[,3])),3] <- 0
      
      results[i, ] <- rev(sort(tapply(temp[,3], temp[,"lg"], max)))[1:m.peak]
      
      setTxtProgressBar(pb, i/n.perm)
    }
    close(pb)
    colnames(results) <- paste("peak", 1:ncol(results), sep=".")
    class(results) <- c("fullsib_perm", "matrix")
    return(results)
  }
  
  
  ##removing NA elements
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]
  #print(pheno.index)
  if(icim == TRUE) pheno.icim <- pheno

  ##preparing for window size
  cim.models <- model_select(fullsib, ws = ws)
  scan.final <- NULL

  if(verbose.scan == TRUE){   
    cat("QTL mapping for ", length(lg)," groups\n")
    count.lg <- 1
    pb <- txtProgressBar(style=3)
    setTxtProgressBar(pb, 0)
  }

  for (LGs in lg){ #scan to the genome

    if(verbose == T)
      cat("\nLinkage Group ", LGs, "under analysis:\n")
    dim_probs <- dim(fullsib$probs[[LGs]]$cond.prob)
    if(is.null(dim_probs)==TRUE)
      next
    ##just a check, group with 2 markers with rf=0.0 falls here...
    else {
      n.pos <- get("dim_probs")[2]    
      scan.results <- rep(NA, n.pos)

      for (mods in 1: length(cim.models[[LGs]])){
      ##for (mods in 2: length(cim.models[[LGs]])){
        #print(cbind(fullsib$pheno[,pheno.col],fullsib$cofactors$matrix.cof[,1:5]))
        if(cim.models[[LGs]][[mods]]$mkr2rm[1] == "none")
          cim.cof.matrix <- fullsib$cofactors$matrix.cof[pheno.index,]
        else{
          #print(fullsib$cofactors$matrix.cof)
          mkr2rm <- do.call("c",lapply(cim.models[[LGs]][[mods]]$mkr2rm,
               function(x) which(colnames(fullsib$cofactors$matrix.cof) == x)))

          cim.cof.matrix <- fullsib$cofactors$matrix.cof[pheno.index,-mkr2rm]#[,-mkr2rm]
        }
        #print(cim.cof.matrix)

        if(is.null(dim(cim.cof.matrix)))
            dim(cim.cof.matrix) <- c(length(cim.cof.matrix),1)
                
        ##finding the matrix used as cofactor
        
        unar <- matrix(1,nrow(cim.cof.matrix))
        colnames(unar) <- "intercept"
        cim.cof.matrix <- cbind(unar,cim.cof.matrix)


        if(icim == FALSE){
          ##traditional analysis
          Zls <- t(qr.solve(crossprod(cim.cof.matrix),t(cim.cof.matrix),tol=1E-30))
          #Zls <- t(solve(crossprod(cim.cof.matrix),t(cim.cof.matrix)))
          gamma <- coef(lm(pheno ~ cim.cof.matrix -1)) #-1 because intercept is already in cim.cof.matrix
        }
        else{
          ##icim proposed -  remove cofactors from phenotype
          pheno.rm <- lm(pheno.icim ~ cim.cof.matrix - 1)
          pheno <- pheno.rm$residuals
          gamma <- coef(pheno.rm)[1] #just the intercept
          Zls <- t(solve(crossprod(unar), t(unar)))
          cim.cof.matrix <- unar
        }
        ##selecting the positions to scan
        #print(dim(Zls))
        cM2scan <- cim.models[[LGs]][[mods]]$position
        if(verbose == TRUE){
          cat("\n\nPositions to scan: ", length(cM2scan),"\n")
          cat(cM2scan,"\n")
        }
        #print(cM2scan)
        ##cM2scan <-  160
        #cM2scan <- 17
        
        geno.class.tmp <- fullsib$probs[[LGs]]$colin$geno.class[cM2scan]
        colin.type.tmp <- fullsib$probs[[LGs]]$colin$type[cM2scan]
        cond.prob.tmp <-  fullsib$probs[[LGs]]$cond.prob[pheno.index,cM2scan,]
        #print(fullsib$probs[[LGs]]$cond.prob[pheno.index, 1, ])
        #print(fullsib$probs[[LGs]]$cond.prob[, 1, ])
        ##print(gamma)
        ##doing the CIM for a lg
        ##print(cim.cof.matrix)
        z <- .C("R_scan_qtl",
                as.integer(length(pheno.index)),
                as.integer(length(cM2scan)),             
                as.integer(geno.class.tmp), 
                as.integer(colin.type.tmp),
                as.double(cond.prob.tmp),#
                as.double(cim.cof.matrix),#
                as.integer(ncol(cim.cof.matrix)),
                as.double(gamma),
                as.double(Zls), #
                as.double(pheno),         #
                result=as.double(rep(0, length(cM2scan))),#ln(log-lk)
                as.integer(maxit),
                as.double(tol),
                as.integer(verbose), PACKAGE="fullsibQTL")$result
        ##as.integer(verbose))$result     
        scan.results[cM2scan] <- z
      }
      ##print(scan.results[cM2scan])
      if(LOD ==TRUE)
        scan.results <- scan.results/log(10)
      else
        ##to check: ta dando bosta: acho que tem q dar o famoso sinal de menos na frente do log10...
        scan.results <- -log10(pchisq(2*scan.results, (fullsib$probs[[LGs]]$colin$geno.class - 1), lower.tail=FALSE))
                                        
      scan.results <- cbind(lg = rep(LGs, length(fullsib$probs[[LGs]]$newmap$dist)),
                            pos.cM = fullsib$probs[[LGs]]$newmap$dist,
                            scan.results, fullsib$probs[[LGs]]$colin$type)
      scan.final <-rbind(scan.final, scan.results)
    }
    if(verbose.scan == TRUE){
      setTxtProgressBar(pb, count.lg/length(lg))
      count.lg <- count.lg + 1
    }
  }
  if(verbose.scan == TRUE) close(pb)

  if(LOD == TRUE)
    colnames(scan.final) <- c("lg", "pos.cM", "LOD", "model")
  else
    colnames(scan.final) <- c("lg", "pos.cM", "-log10(pval)", "model")
  class(scan.final) <- c("fullsib_scan", "matrix")
  scan.final
}
