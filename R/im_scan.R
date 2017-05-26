 
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
# copyright (c) 2011, Rodrigo Gazaffi                                 #
# im_scan is based on scanone function from qtl package               #
#                                                                     #
# First version: 09/30/2011                                           #
# Last  version: 09/30/2011                                           #
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

im_scan <- function (fullsib,lg, pheno.col=1, addcovar=NULL, LOD=TRUE,
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


  
