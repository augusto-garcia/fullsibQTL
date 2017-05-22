#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: cof.definition.R                                              #
#                                                                     # 
# Contains:                                                           #
# cof.definition                                                      #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
#                                                                     #
# First version: 09/30/2011                                           #
# Last  version: 09/30/2011                                           #
# License: GPL-3                                                      #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: cof.definition.R                                          #
#                                                                     # 
# Allows that one defines their own regions considered as cofactors,  #
# for example, if one wants to considers cofactor by peak  or some    #
# ad-hoc decision                                                     #
#######################################################################

cof.definition <- function(fullsib, pheno.col=1, addcovar=NULL,
                           cof.pos, thres.effect=1){

  ##checking arguments
  if (!any(class(fullsib) == "fullsib"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib'\n")

  ##checking arguments
  if (is.na(match(pheno.col, 1:ncol(fullsib$pheno))) ) 
   stop("pheno.col must be a number between 1 and", ncol(fullsib$pheno), "\n")
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]


  ##checking arguments
  if(!is.null(addcovar)){
    if(!is.matrix(addcovar)) stop("addcovar argument must be a matrix\n")
    else{
      covar <- addcovar
      colnames(covar) <- rep("covar",ncol(covar))
    }
  }   
  
  ##checking arguments
  if(missing(cof.pos))
    stop("cof.pos argument is missing\n")

  ##checking arguments
  if(is.null(dim(cof.pos)))
    stop("cof.pos should be a matrix\n")

  ##checking arguments
  if((thres.effect < 0 ) || (thres.effect > 1))
    stop("thres.effect argument should be a value between 0 and 1\n")

  ##obtaing the cofactors matrix
  matrix.cof <- NULL
  for (i in 1:nrow(cof.pos)){
    lg <- as.numeric(cof.pos[i,1])
    pos <-  cof.pos[i,2]
    if(is.na(lg)){
      ##unlinked marker
      pos.num <- match(pos, colnames(get(fullsib$map[[1]]$data.name)$geno))
      tmp <- unlink2cof(pos.num, fullsib$map[[1]]$data.name, pheno.index)[[1]]
      colnames(tmp) <- rep(pos, ncol(tmp))
      
    }
    else{
      ##marker on the map
      pos.num <- match(pos, names(fullsib$probs[[lg]]$newmap$dist))
      if(is.na(pos.num))
        stop("Markers ", dQuote(cof.pos[i,2]), " was not found: Please check it\n")
      tmp <- prob2cof(fullsib$probs[[lg]]$cond.prob[, pos.num,],
                      fullsib$probs[[lg]]$colin$type[ pos.num] )
      tmp <- as.matrix(tmp[pheno.index,])
      colnames(tmp) <- rep(pos, ncol(tmp))
    }
    matrix.cof <- cbind(matrix.cof, tmp)
  }

  ## check if some effects could not be estimated: coef==NA
  check.na <- which(is.na(coef(lm(pheno ~ matrix.cof)))) 
  if(length(check.na) > 0) matrix.cof <- matrix.cof[,-check.na]
  
  ## check for non-significative effect (pvalue higher then thres.effect),
  ## if happens, remove from the matrix and cof.list
  label.matrix <- colnames(matrix.cof)

  check.na <- which(is.na(coef(lm(pheno ~ matrix.cof))))
  if(length(check.na) > 0) matrix.cof <- matrix.cof[,-(check.na - 1)]

  ## Elimination of non significative mkrs, using thres.effects argument
  rm.ns.effects <-
    which(summary(lm(pheno ~ matrix.cof))$coefficients[-1,"Pr(>|t|)"] < thres.effect)
  
  matrix.cof <- as.matrix(matrix.cof[,rm.ns.effects])
  colnames(matrix.cof) <- label.matrix[rm.ns.effects]

  mk.kept <- match(cof.pos[,2], unique(colnames(matrix.cof)))
  mk.kept1 <- which(is.na(mk.kept))

  ##checking
  if(length(mk.kept) == length(mk.kept1))
    warning("None cofactors included are significative, if 'thres.effect' used is ", thres.effect)

  else {
    if(length( which(is.na(mk.kept))) > 0)
      cof.pos <- cof.pos[ -which(is.na(mk.kept)), ]
    
    
    if(is.null(dim(cof.pos)))
      cof.list <- data.frame(lg = as.numeric(cof.pos[1]),
                             mkr= cof.pos[2])  
    else     
      cof.list <- data.frame(lg = as.numeric(cof.pos[,1]),
                             mkr= cof.pos[,2])
    
    matrix.cof2 <- matrix(NA, ncol=ncol(matrix.cof), nrow=nrow(fullsib$pheno))
    matrix.cof2[pheno.index,] <- matrix.cof
    colnames(matrix.cof2) <- colnames(matrix.cof)

    if(!is.null(addcovar)){
      matrix.cof2 <- cbind(covar, matrix.cof2)
      colnames(matrix.cof2)[1:ncol(addcovar)] <- "covar"
    }
    
    
    fullsib$cofactors <- list(names.cof=cof.list,
                              matrix.cof = matrix.cof2,
                              trait.cof = pheno.col)
    
      class(fullsib) <- c("fullsib", "fullsib.cofactors")
  }
  fullsib
}##end-function
