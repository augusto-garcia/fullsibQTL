#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: r2.ls.R                                                       #
#                                                                     # 
# Contains: r2.ls                                                     #
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
# Function: r2.ls                                                     #
#                                                                     #
# It provides an estimation of R2.LS, i.e., part of phenotypic        #
# variation explained by the QTL.                                     #
# Attention: these estimations are done using linear square aproach   #
# and not using MLE, it should be interpredted with some attention    #
#                                                                     #
#######################################################################

r2.ls <- function(fullsib, lg, pos, pheno.col=1, addcovar = NULL, ls.estimation=FALSE){

  ##checking arguments
  if (!any(class(fullsib) == "fullsib")) 
    stop(sQuote(deparse(substitute(fullsib))), " is not an object of class 'fullsib'")

  if (missing(lg))
    stop("lg argument must be indicated")
  
  if (missing(pos))
    stop("pos argument must be indicated")

  if (length(lg) != length(pos))
    stop("lg and pos arguments must be same length")

  if (is.null(addcovar))
    full.mat <- NULL
  else{
    full.mat <- addcovar
    colnames(full.mat) <- "addcovar"
  }

  ##creating genotype matrix
  for(i in 1:length(lg)){
    chr <- lg[i]
    position <- which(colnames(fullsib$probs[[ chr ]]$cond.prob) == pos[i])

    tmp.mat <- prob2cof(fullsib$probs[[ chr ]]$cond.prob[,position,],
    as.character(fullsib$probs[[ chr ]]$colin$type[position]))
    colnames(tmp.mat) <- rep(paste("lg", lg[i],".", pos[i], sep=""), ncol(tmp.mat))
    full.mat <- cbind(full.mat, tmp.mat)
  }

  ##organizing phenotypes and removing NA individuals
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]
  full.mat <- full.mat[pheno.index,]

  ##fitting lm for all QTLs
  lm.full <- lm(pheno  ~ full.mat)
  if (is.null(addcovar))
    lm.null <- lm(pheno ~ 1)
  else
    lm.null <- lm(pheno ~ addcovar[pheno.index])

  ##correto usar sum of squares??? ou quadrados medios ou dividir por SS/n
  ## for all qtls together
  ss.full <- anova(lm.full)["Residuals", 2]
  ss.null <- anova(lm.null)["Residuals", 2]
  r2.trait <- (ss.null - ss.full) / ss.null
  r2.trait

  ##for each QTL
  if(!is.null(addcovar)){  
    qtls <- unique(colnames(full.mat))
    qtls <- qtls[which(unique(colnames(full.mat)) != "addcovar")]
  }
  else
    qtls <- unique(colnames(full.mat))

  r2.qtls <- rep(NA, length(qtls))
  for (i in 1:length(qtls)){
    ef2rm <- which(qtls[i] == colnames(full.mat))
    ss.qtl <- anova(lm(pheno ~ full.mat[,-ef2rm]))["Residuals", 2]
    r2.qtls[i] <- (ss.qtl - ss.full)/ss.null
  }

  ##print(summary(lm.full))

  if(ls.estimation){
    cat("\nLinear Regression Estimatives using Least Square approach, considering all QTL\n\n")
    print(summary(lm.full)$coefficients)
  }
  
  r2.output <- round(c(r2.trait, r2.qtls) * 100, 4)
  names(r2.output) <- c("R2.trait", paste("R2", qtls, sep="."))
  r2.output
}
