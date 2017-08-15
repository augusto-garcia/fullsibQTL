#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: cof_selection.R                                               #
#                                                                     # 
# Contains:                                                           #
# cof_selection                                                       #
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
# Function: cof_selection                                             #
#                                                                     #
# Performs the selection of significative markers that is used as     #
# cofactors in CIM analysis. It receives all the information, prepares#
# the formulas and model selection is done by step function. Finally  #
# if deseared  a second selection is done ocnsidering the p-values os #
# selected effects, which are not significant considering a threshold #
# value indicated by user is removed from the model.                  #
#######################################################################

## -------------------------
## cof_selection function

#' Cofactors Selection procedure for CIM analysis
#' 
#' Performs the selection of markers that will be used as cofactors in CIM
#' analysis. The procedure is fully based on \code{step} function from
#' \pkg{stats} package
#' 
#' This function prepares the information contained in the object of class
#' \emph{fullsib} for using on a multiple linear regression. By default the
#' independent variable are molecular markers that are present on the linkages
#' groups (LG), optionally one can also consider the markers that were not
#' placed on the map. Genetic Predictors are obtained by the multipoint
#' probability at the exactly point that the markers is placed on the LG and
#' multiplied by the estimable contrasts.
#' 
#' The genetic predictors are regressed with the phenotype using the
#' \code{step} function. This strategy was done to get advantages on model
#' selection tools avaliable on . In this case, the function performs stepwise
#' procedure using information criteria, such as AIC and BIC (defined by the
#' user). If one would like to use P-value (or another criteria) to select
#' cofactors, one can do by oneself and add these cofactors using
#' \code{cof_definition} function.
#' 
#' As the selection process is based on \code{step}, one can follow the
#' selection just using \code{trace} argument (same present in \code{step})
#' higher than zero.
#' 
#' Remembering that this is an intermediate step for CIM analysis, and
#' sometimes many markers can be selected as cofactors reducing statistical
#' power. To solve this situation, it is considered some arguments that can
#' control the number of included markers, for example \sQuote{n.cofactor}
#' controls the maximum number of included markers established by the user.
#' Other option is \sQuote{stoppage.df}, which the user can define the maximum
#' number of degree of freedom that model can have. In both case, the selection
#' criteria is used to include cofactors, but when model reaches one of these
#' limits FIRST the process is stopped.
#' 
#' An other way to control the selection of markers is to verify the non
#' significative effects. For this, a linear model is fitted \code{lm} function
#' and \code{summary.lm} is used to obtain the p-values of the cofactors
#' effects. So \sQuote{thres.effect} is a threshold probability used to
#' determine which effect are significative. If 1 is used none effect is
#' considered non significative and removed from the analysis.
#' 
#' The method plot was designed for the user evaluated the cofactors saturation
#' and dispersion along the genome. Helping with the decision of window size
#' dimension and best selection options.
#' 
#' @aliases cof_selection
#' @param fullsib an object from class \emph{fullsib}.
#' @param pheno.col Column number in the phenotype matrix (present in
#' \emph{fullsib} object) which should be used as the phenotype.
#' @param addcovar Additive covariates. If it is used, one must indicate the
#' design matrix for those source of variation. It should be noted that
#' additive covariates is included in the model as fixed effects, under
#' ordinary linear regression.
#' @param k value to be used as penalty in the model selection. If \code{k=2}
#' is used corresponds to the Akaike Information Criterion; for \code{k=log(n)}
#' is the Bayesian information criterion (BIC) or sometimes called Schwarz
#' Information criterion.
#' @param n.cofactor maximum number of cofactors that will be present on the
#' model. Default is to consider a model with maximum of 10 cofactors. See
#' details.
#' @param stoppage.df controls the maximum number of markers to be included in
#' the model, using the model's degree of freedom as a limit to stop the
#' process. Default is NULL (not considering this argument), but
#' \eqn{2\times\sqrt(n)} (being \eqn{n} is number of individuals phenotyped)
#' can be used as a limit to control the maximum number of cofactors.
#' @param trace If zero is considered none information about the selection
#' process is showed, if non-negative integer is considered, one can see
#' additional information during the running process. This argument is the same
#' from \code{step}, once the model selection is based on it.
#' @param thres.effect Threshold value to remove non significative effect for
#' selected cofactors. Default is \code{thres.effect = 1}, none cofactors
#' effect is removed from the analysis.
#' @param selection integer indicating which markers should be used for this
#' step. It may assumes values of \sQuote{0}, \sQuote{1} and \sQuote{2}. If
#' \sQuote{0}: only markers that are not included in the map are used on the
#' selection process. If \sQuote{1}: the selection process just consider
#' markers that are on the linkage groups. And, if \sQuote{2}: all markers are
#' used for selection model.
#' 
#' @return An object of class \emph{fullsib_cofactors} returned, which has the
#' same structure of an object of class \emph{fullsib} with the inclusion of an
#' extra component (\sQuote{cofactor}) that is a list with the components:
#' \code{names.cof}, \code{matrix.cof} and \code{trait.cof}.
#' 
#' \code{names.cof} is a data frame showing which are the selected markers and
#' their linkage groups. \code{NA} value is used for markers that are not
#' placed on linkage groups.
#' 
#' \code{matrix.cof} is a matrix with contains all the cofactors effect. During
#' the CIM analysis columns are dropped in function of window size.
#' 
#' \code{trait.cof} is the indication for which pheno.col was considered for
#' cofactor selection step.
#' @author Rodrigo Gazaffi, \email{rgazaffi@@gmail.com}
#' \code{plot_fullsib_cofactors } was based on \code{draw_map} that was written
#' by Marcelo Mollinari \email{mmollina@@usp.br} and modified by Rodrigo
#' Gazaffi, allowing cofactors inclusion.
#' @seealso 
#' \code{\link[fullsibQTL]{create_fullsib}}
#' \code{\link[fullsibQTL]{cof_definition}}
#' 
#' \code{step} from \pkg{stats} package
#' @references Burnham, K.P., Anderson, D.R. (2004) Understanding AIC and BIC
#' in Model Selection. \emph{Sociological Methods & Research}, 33: 261-304
#' 
#' Sakamoto, Y., Ishiguro, M., and Kitagawa, G. (1986), Akaike
#' \emph{Information Criterion Statistics}, Tokyo: KTK Scientific Publishers.
#' @keywords utilities
#' 
#' @examples
#'   data( "example_QTLfullsib" )
#' 
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ),
#'                              step = 0, map.function = "kosambi" , condIndex = 3.5 )
#' 
#' 
#' 
#'   ###############################################
#'   ## cofactor selection using BIC (n.ind = 300)
#' 
#'   ### just using markers that are placed on linkage groups (default)
#'   cofs.fs <- cof_selection( fullsib, pheno.col = 1, k = log( 300 ), selection = 1 )
#'   cofs.fs
#'   plot( cofs.fs )
#' 
#' 
#'   ###just markers that are not placed in LG
#'   ###Attention: just used if there are markers that are not placed on
#'   ###genetic map
#'   \dontrun{
#'   cofs.fs1 <- cof_selection( fullsib, pheno.col = 1, k = log( 300 ), selection = 0 )
#'   }
#'   ###all markers are used
#'   ###Attention: just used if fullsib$unlinked.mkr is not NULL
#'   \dontrun{
#'    cofs.fs2 <- cof_selection( fullsib, pheno.col = 1, k = log( 300 ), selection = 2 )
#'    }
#' 
#' 
#' 
#'   ###############################################
#'   ## cofactor selection using AIC criteria
#' 
#'   \dontrun{  
#'   ### using AIC alowing at maximum of 20 markers
#'   cofs.fs3 <- cof_selection( fullsib, pheno.col = 1, k = 2, n.cofactor =  20 )
#'   plot( cofs.fs3 )
#' 
#'   ### using AIC alowing at maximum of 5 markers
#'   cofs.fs4 <- cof_selection( fullsib, pheno.col = 1, k = 2, n.cofactor = 5 )
#'   plot( cofs.fs4 )
#' 
#'   ### using AIC alowing the model with maximum of 18 D.f.
#'   ### here 6 marker, because each cofactor has 3 effects
#'   cofs.fs5 <- cof_selection( fullsib, pheno.col = 1, k = 2, stoppage.df = 18 )
#'   plot( cofs.fs5 )
#' 
#'   ### it is possible to see the selection process
#'   cofs.fs6 <- cof_selection( fullsib, pheno.col = 1, k = 2, n.cofactor = 5, trace = 1 )
#'   plot( cofs.fs6 )
#' 
#'   ### Same selection, but removing effects that are not signif. at 5%
#'   cofs.fs7 <- cof_selection( fullsib, pheno.col = 1, k = 2, n.cofactor = 5, thres.effect = 0.05 ) 
#'   }
#' 
#'   ##cofactor selection using additive covariate
#'   \dontrun{
#'   covar <- matrix( rep( c( 1, -1 ), each = 150 ), ncol = 1 )
#'   cofs.fs8 <- cof_selection( fullsib, pheno.col = 2, addcovar = covar, k = 2 )
#'   }
#' 

cof_selection <- function(fullsib, pheno.col=1, addcovar=NULL, k=2,
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
      rm.covar <- which(substring(names.cof, 1,5) == "covar")

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
      class(fullsib) <- c("fullsib", "fullsib_cofactors")
    }
  }
  return(fullsib)
} #end function
