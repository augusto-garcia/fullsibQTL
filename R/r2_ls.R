#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: r2_ls.R                                                       #
#                                                                     # 
# Contains: r2_ls                                                     #
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
# Function: r2_ls                                                     #
#                                                                     #
# It provides an estimation of R2.LS, i.e., part of phenotypic        #
# variation explained by the QTL.                                     #
# Attention: these estimations are done using linear square aproach   #
# and not using MLE, it should be interpredted with some attention    #
#                                                                     #
#######################################################################

## -------------------------
## r2_ls function

#' Calculates the \eqn{R^2} for a group of given QTL.
#' 
#' Considering a group os mapped QTL using either \code{im_scan} or
#' \code{cim_scan}, this function calculates the proportion of all QTL together
#' (and also separetely) explain from total phenotipic variation, using least
#' square estimation.
#' 
#' This function is an alternative way to obtain a estimative of \eqn{R^2} for
#' a model considering simultaneously all QTL mapped on genome scan. This value
#' could be interpreted as an estimation of heritability on the wide sense.
#' 
#' Here, it is not the main goal to develop multiple QTL model using mixture
#' models, so if one really wants this value, we provide an estimation using
#' the least square approach, which is a very good approximation on the results
#' obtained with genome scan.
#' 
#' To obtain these values we first fit a full model with all QTL declared and a
#' second one with no QTL (null model). We obtain the Residual Sum of Squares
#' (\eqn{RSS}) of both models and calculated the \eqn{R^2} as: \deqn{R^2 =
#' \frac{RSS.null - RSS.full}{RSS.null}}
#' 
#' To estimate \eqn{R^2} for the \eqn{i-th} QTL, we obtain the \eqn{RSS} for a
#' model with all QTL dropping the QTL \eqn{i} and estimate: \deqn{R^2_i =
#' \frac{RSS.drop_i - RSS.full}{RSS.null}}
#' 
#' If \code{ls.estimation=TRUE}, the estimation for the model with all QTL
#' adjusted (full model) is printed. This can be done if one would like to
#' check if the estimatives with least square and the ones obtained with genome
#' scan. The values obtained here may differ from the one with genome scan, but
#' not about to change its significance.
#' 
#' @param fullsib An object from class \emph{fullsib}.
#' @param lg Vector of integers indicating the linkage groups for each QTL.
#' @param pos Vector of character (of same length as \code{lg}) indicating the
#' locus that will be taken. This name is found as the row label for the matrix
#' returned by \code{im_scan} or \code{cim_scan}.
#' @param pheno.col Column number in the phenotype matrix (present in
#' \emph{fullsib} object) which should be used as the phenotype.
#' @param addcovar Additive covariates. If it is used, one must indicate the
#' design matrix for those source of variation. It should be noted that
#' additive covariates is included in the model as fixed effects, under
#' ordinary linear regression.
#' @param ls.estimation if \code{TRUE} prints the summary of multiple linear
#' regression for all QTL providing the effects and their significance. If
#' \code{FALSE} no extra details is given.
#' 
#' @return It returns a vector with \eqn{q+1} position (\eqn{q} is the number
#' of QTL) with \eqn{R^2} estimation. In the first position one finds the value
#' for a model fitted with all QTLs simultaneously (identified as
#' \sQuote{R2.trait}), the others positions provides the \eqn{R^2} for each QTL
#' separately (identified as \sQuote{R2.lgX.Y}, where lgX.Y means the lg number
#' X in the Y-th position)
#' 
#' @author Rodrigo Gazaffi \email{rgazaffi@@gmail.com}
#' 
#' @seealso 
#' \code{\link[fullsibQTL]{cim_char}}
#' \code{\link[fullsibQTL]{cim_scan}}
#' \code{\link[fullsibQTL]{im_char}}
#' \code{\link[fullsibQTL]{im_scan}}
#' 
#' @keywords utilities
#' 
#' @examples
#'   data(example_QTLfullsib)
#' 
#'   fullsib <- create_fullsib(example_QTLfullsib,
#'                             list(LG1_final, LG2_final, LG3_final, LG4_final),
#'                             step=0,map.function="kosambi",condIndex=3.5)
#' 
#'   r2_ls(fullsib, pheno.col=1, lg=c(1,1,2,2,2,3,4),
#'      pos=c("M3","M12","M20","M24","M27","M37","M52"))
#' 
#'   covar <- matrix(rep(c(1,-1), each=150), ncol=1)
#'   r2_ls(fullsib, pheno.col=2, addcovar = covar, lg=c(1,1,2,2,2,3,4),
#'      pos=c("M3","M12","M20","M24","M27","M37","M52"),
#'      ls.estimation=TRUE) 
#'

r2_ls <- function(fullsib, lg, pos, pheno.col=1, addcovar = NULL, ls.estimation=FALSE){

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
  r2.output <- data.frame(lg=c(NA,lg),loc=c("Intercept",pos),r2=r2.output)
  class(r2.output) <- c("data.frame","r2ls_out")
  return(r2.output)
}