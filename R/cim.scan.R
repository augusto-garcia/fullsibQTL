#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: cim.scan.R                                                    #
#                                                                     # 
# Contains: cim.scan                                                  #
#                                                                     #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
#                                                                     #
# First version: 09/30/2011 (american date format)                    #
# Last  version: 09/30/2011 (american date format)                    #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: cim.scan                                                  #
#                                                                     #
# Performs QTL mapping on the genome using Composite Interval Mapping #
# This function is based on scanone from qtl pkg written by K. Broman #
#                                                                     #
# It does two basic analysis:                                         #
# if n.perm argument is missing or is 0, the QTL mapping is done.     #
# The result is an object from fullsib.scan class and also a matrix   #
# for xtable users.                                                   #
#                                                                     #
# If n.perm > 0, this function can perform permutation test           #
# The result is an object from fullsib.perm class.                    #
# EM algorithm is done in C code                                      #
#######################################################################

cim.scan <- function(fullsib, lg, pheno.col=1, ws = 10, LOD=TRUE,
                     maxit=1000, tol=1e-04, n.perm=0, write.perm,
                     icim=FALSE, verbose, verbose.scan)
{

  ##checking arguments
  if (!any(class(fullsib) == "fullsib.cofactors"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib.cofactors': execute 'cof.selection' or 'cof.definition' first")

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

      temp <- cim.scan(fullsib, lg=lg, pheno.col=pheno.col, ws=ws, maxit=maxit, tol=tol,
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
    class(results) <- c("fullsib.perm", "matrix")
    return(results)
  }
  
  
  ##removing NA elements
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]
  #print(pheno.index)
  if(icim == TRUE) pheno.icim <- pheno

  ##preparing for window size
  cim.models <- model.select(fullsib, ws = ws)
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
  class(scan.final) <- c("fullsib.scan", "matrix")
  scan.final
}
