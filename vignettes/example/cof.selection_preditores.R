#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: cof.selection.R                                               #
#                                                                     # 
# Contains:                                                           #
# cof.selection                                                       #
# plot.fullsib.cofactors                                              #
# print.fullsib.cofactors                                             #
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
# Function: cof.selection                                             #
#                                                                     #
# Performs the selection of significative markers that is used as     #
# cofactors in CIM analysis. It receives all the information, prepares#
# the formulas and model selection is done by step function. Finally  #
# if deseared  a second selection is done ocnsidering the p-values os #
# selected effects, which are not significant considering a threshold #
# value indicated by user is removed from the model.                  #
#######################################################################

cof.selection <- function(fullsib, pheno.col=1, addcovar=NULL, k=2,
                          n.cofactor=10, stoppage.df=NULL, trace=0,
                          thres.effect=1, selection=1)
{
    
    ## checking arguments
    if (!any(class(fullsib) == "fullsib"))
        stop(sQuote(deparse(substitute(fullsib))),
             " is not an object of class 'fullsib'")
    
    
    ## checking arguments
    if(is.na(match(selection, c(0,1,2))) )
        stop("selection argument must be 0, 1 or 2")
    if( (selection == 0) || (selection == 2))
        if(is.null(fullsib$unlinked.mkr))
            stop("selection argument can be 0 or 2, if unlinked markers are included in fullsib")
    
    
    ## checking arguments
    if((thres.effect < 0 ) || (thres.effect > 1))
        stop("thres.effect argument should be a probability between 0 and 1")

    ##checking arguments and preparing phenotype for analysis
    if (is.na(match(pheno.col, 1:ncol(fullsib$pheno))) )
        stop("pheno.col must be a number between 1 and", ncol(fullsib$pheno))
    pheno <- fullsib$pheno[,pheno.col]
    pheno.index <- which(!is.na(pheno))
    pheno <- pheno[pheno.index]
    
    ## checking arguments
    if(!is.null(addcovar)){
        if(!is.matrix(addcovar)) stop("addcovar argument must be a matrix")
        else
            if(ncol(addcovar) == 1) covar <- as.matrix(addcovar[pheno.index])
            else covar <- addcovar[pheno.index,]
    }
              
    ## genetic predictors for stepwise precedure (which are in the map)
    total.lg <-  length(fullsib$map)
    tmp.cof.sel <- vector("list", total.lg + 1)
    names(tmp.cof.sel) <- paste("LG", 1:(total.lg + 1),sep="")

    if(selection > 0){ #1 or 2
        for ( i in 1:total.lg){
            mkr.num <- which(substring(names(fullsib$probs[[i]]$newmap$dist),1,3) != "loc")
            tmp.cof.sel[[i]] <- vector("list", length(mkr.num))
            names(tmp.cof.sel[[i]]) <- paste("mkr", fullsib$map[[i]]$seq.num, "cof", sep=".")
            for ( j in 1:length(tmp.cof.sel[[i]]) )
                tmp.cof.sel[[i]][[j]] <-
                    prob2cof(fullsib$probs[[i]]$cond.prob[ pheno.index , mkr.num[j],],
                             fullsib$probs[[i]]$colin$type[ mkr.num[j] ])  
        }
    }
    
    ##geentic preditors for markers that are not in the map...
    if(selection != 1){ #0 or 2
        mk.unmap <- fullsib$unlinked.mkr
        obj.name <- fullsib$map[[1]]$data.name
        if(!is.null(fullsib$unlinked.mkr))
            tmp.cof.sel[[total.lg+ 1]] <- unlink2cof(mk.unmap, obj.name, pheno.index)
    }
    
    cof.pred <<-tmp.cof.sel
    ######################################################
    cof.all <<- vector("list", total.lg + 1)
    for ( i in 1:total.lg){
        chr.pos <- length(fullsib$probs[[i]]$newmap$dist)
        cof.all[[i]] <<- vector("list", chr.pos)
        for (j in 1:chr.pos)
            cof.all[[i]][[j]] <<- prob2cof(fullsib$probs[[i]]$cond.prob[ pheno.index , j,], 
                                           fullsib$probs[[i]]$colin$type[ j ])  
    }
    ######################################################
    switch(EXPR=as.character(selection),
           "0" = seq.mk <- total.lg + 1,
           "1" = seq.mk <- 1:total.lg,
           "2" = seq.mk <- 1:(total.lg+1)
           )
    
    ## preparing the formulas for stepwise selection
    formula <- NULL
    for (i in seq.mk)
        formula <- c(formula, paste("tmp.cof.sel",
                                    names(tmp.cof.sel)[i], names(tmp.cof.sel[[i]]),
                                    sep="$", collapse="+"))
    formula <- paste(formula, collapse="+")
    
    if(is.null(addcovar)){
        max.form <- paste("pheno", formula, sep=" ~ ")
        min.form <- paste("pheno", "1", sep=" ~ ")
    }
    else{
        max.form <- paste("pheno", paste("covar",formula, sep=" + "), sep=" ~ ")
        min.form <- paste("pheno", "covar", sep=" ~ ")
    }
    
    ##
    cof.max <<-max.form
    cof.min <<-min.form
    ##
    
    ## doing the stepwise selection
    max.form <- as.formula(max.form)
    min.form <- as.formula(min.form)
    
    check.DF <- 0
    check.DF.old <- -1
    n.iter <- 1
    
    get.env <- environment()
    environment(step) <- get.env


  ##if is NULL add until the maximum number of individuals
  ##here the algorithm stops considering Inf. Crit. or n.cofactor
  if(is.null(stoppage.df))
    stoppage.df <- length(pheno.index)
  
  if (trace == 0) cat("Number of Cofactors selected: ")
  while( (check.DF < stoppage.df) && (n.iter <= n.cofactor) ){
      
    min.form <- step(lm(min.form),scope=list(upper=max.form, lower=min.form),
                     direction = "both", k=k, step=1, trace = trace)
    check.DF <- sum(anova(min.form)$Df) - df.residual(min.form)
    min.form <- as.formula(min.form)
    ##To avoid an infinity loop, if the model do not select anything more
    if (check.DF == check.DF.old) break
    check.DF.old <- check.DF
    if (trace == 0) cat( n.iter, "... ")
    n.iter <-  n.iter + 1
  }
    if (trace == 0) cat("done\n")
  cofactors <- lm(min.form)
  
  matrix.cof <- model.matrix(cofactors)
  names.cof <- names(coef(cofactors))
  
  if(!is.null(addcovar))  {
    rm.covar <- which(names.cof == "covar")
    matrix.cof <- matrix.cof[,-rm.covar]
    names.cof <- names.cof[-rm.covar]
    #print(rm.covar)
    #print(head(matrix.cof))
    #print(names.cof)
  }

  if(length(names.cof) == 1){
    ##check point: if just intercept no markers was selected
    warning("  Any markers were selection as cofactors, considering information criteria as ", sQuote(deparse(substitute(k))))
  }
   else {
     
    rm.int <- which(names.cof == "(Intercept)")
    matrix.cof <- matrix.cof[,-rm.int]
    names.cof <- names.cof[-rm.int]

    names.cof <- do.call("rbind",strsplit(names.cof, split="$", fixed=TRUE))[,2:3]

    if(is.null(dim(names.cof)))
      get.cofs <- do.call("rbind",
                          strsplit(names.cof[2], split=".",fixed=TRUE))[,2]
    else
      get.cofs <- do.call("rbind",
                          strsplit(names.cof[,2], split=".",fixed=TRUE))[,2]

    get.cofs <- as.numeric(get.cofs)
    cof.index <- which(!duplicated(get.cofs))

    selected.cofs <- colnames(get(fullsib$map[[1]]$data.name)$geno)[get.cofs]
    extract.cofs <- selected.cofs[which(!duplicated(selected.cofs))]


    if(is.null(dim(names.cof))){
      extract.lg <- as.numeric(substring(names.cof[1],3,
                                         nchar(names.cof[1])))
    }
    else{
      extract.lg <- as.numeric(substring(names.cof[cof.index,1],3,
                                         nchar(names.cof[cof.index,1])))
    }
    cof.list <- data.frame(lg = extract.lg, mkr= extract.cofs)
    cof.list <- cof.list[order(cof.list[,1], cof.list[,2]),]
    cof.unmap <- which(cof.list$lg  == (length(fullsib$map) + 1))

    if(length(cof.unmap) > 0)
      cof.list$lg[cof.unmap] <- NA

    matrix.cof <- as.matrix(matrix.cof)

    colnames(matrix.cof) <- selected.cofs
    
    ##check if some effects could not be estimated by singularity problems "NA"
    check.na <- which(is.na(coef(lm(pheno ~ matrix.cof))))
    if(length(check.na) > 0) matrix.cof <- matrix.cof[,-(check.na - 1)]
    ## check.na - 1 because intercept column is not used in this matrix
    labels.cofs <- colnames(matrix.cof)
    
    ## check for non-significative effect,
    ## if there is remove from the matrix and cof.list if necessary
    keep.effects <-
      which(summary(lm(pheno ~ matrix.cof))$coefficients[-1,"Pr(>|t|)"] < thres.effect)

    if(length(keep.effects) == 0) {
      warning("  Any markers were selection as cofactors, considering information criteria as ", thres.effect)
    }
    else{
      matrix.cof <- as.matrix(matrix.cof[,keep.effects])

      matrix.cof2 <- matrix(NA, ncol=ncol(matrix.cof), nrow=nrow(fullsib$pheno))
      matrix.cof2[pheno.index,] <- matrix.cof
      matrix.cof <- matrix.cof2
      colnames(matrix.cof) <- labels.cofs[keep.effects]
      keep.effects <-which(as.character(cof.list[,2]) %in% colnames(matrix.cof))
      if(length(keep.effects) > 1)
        cof.list <- cof.list[keep.effects,]

      if(!is.null(addcovar)){
        matrix.cof <- cbind(covar, matrix.cof)
        colnames(matrix.cof)[1:length(rm.covar)] <- "covar"
      }
      
      fullsib$cofactors <- list(names.cof =  cof.list,
                                matrix.cof = matrix.cof,
                                trait.cof = pheno.col)
      class(fullsib) <- c("fullsib", "fullsib.cofactors")
    }
  }
  fullsib
} #end function
