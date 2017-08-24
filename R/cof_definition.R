#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: cof_definition.R                                              #
#                                                                     # 
# Contains:                                                           #
# cof_definition                                                      #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# Updated by Rodrigo Amadeu                                           #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
#                                                                     #
# First version: 09/30/2011                                           #
# Last  version: 08/15/2011                                           #
# License: GPL-3                                                      #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: cof_definition.R                                          #
#                                                                     # 
# Allows that one defines their own regions considered as cofactors,  #
# for example, if one wants to considers cofactor by peak  or some    #
# ad-hoc decision                                                     #
#######################################################################

## -------------------------
## cof_definition function

#' Ad-hoc definition of cofactors for CIM analysis
#' 
#' Defines any marker and/or position to be used as a cofactor. It was designed
#' as an alternative way to include cofactor besides using \code{cof_selection}
#' function.
#' 
#' The standard procedure used to perform cofactor selection on this package is
#' based on multiple linear regression using information criteria. However for
#' any reason, one would like to have more flexibility to control the number
#' and location of cofactors to be added in the model.
#' 
#' One possible way of select cofactors is to perform an interval mapping
#' analysis and with the result, one can perform CIM analysis, i.e., first IM
#' is done to detect QTL. After this, the position of mapped QTL can be added
#' for CIM analysis, using \code{cof_definition}. The inclusion and exclusion
#' of cofactor can be performed for some rounds, until one get a final QTL
#' mapping profile.
#' 
#' The method plot was designed for the user evaluated the cofactors saturation
#' and dispersion along the genome. Helping with the decision of window size
#' dimension and best selection options.
#' 
#' Finally, with the development of \code{cof_definition} and
#' \code{cof_selection}. We believe that user has enough flexibility of dealing
#' with selection of markers to be used in CIM.
#' 
#' @param fullsib An object from class \emph{fullsib}.
#' @param pheno.col Column number in the phenotype matrix (present in
#' \emph{fullsib} object) which should be used as the phenotype.
#' @param addcovar Additive covariates. If it is used, one must indicate the
#' design matrix for those source of variation. It should be noted that
#' additive covariates is included in the model as fixed effects, under
#' ordinary linear regression.
#' @param cof.pos character matrix, in each row represents a different cofactor
#' and two columns indicating the linkage group and the name of the position to
#' be considered. To define markers that are not place on markers, one should
#' indicate as \code{NA}. See examples.
#' @param thres.effect Threshold value to remove non significative effect for
#' selected cofactors. Default is \code{thres.effect = 1}, none cofactors
#' effect is removed from the analysis.
#' 
#' @return An object of class \emph{fullsib_cofactors} returned, which has the
#' same structure of an object of class \emph{fullsib} with the inclusion of an
#' extra component (\sQuote{cofactor}) that is a list with the components:
#' 
#' \code{names.cof}, \code{matrix.cof} and \code{trait.cof}.
#' \code{names.cof} is a data frame showing which are the selected markers and
#' their linkage groups. \code{NA} value is used for markers that are not
#' placed on linkage groups.
#' \code{matrix.cof} is a matrix with contains all the cofactors effect. During
#' the CIM analysis columns are dropped in function of window size.
#' \code{trait.cof} is the indication for which pheno.col was considered for
#' cofactor selection step.
#' 
#' @author Rodrigo Gazaffi, \email{rgazaffi@@gmail.com}
#' 
#' \code{plot_fullsib_cofactors } was based on \code{onemap::draw_map} that was
#' written by Marcelo Mollinari \email{mmollina@@usp.br} and modified by
#' Rodrigo Gazaffi, allowing cofactors inclusion.
#' @seealso \code{\link[fullsibQTL]{create_fullsib}},
#' 
#' \code{\link[fullsibQTL]{cof_selection}}
#' 
#' \code{\link[fullsibQTL]{cim_scan}}
#' 
#' \code{\link[fullsibQTL]{cim_char}}
#' @keywords utilities
#' @examples
#'   data( "example_QTLfullsib" )
#' 
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ),
#'                              step = 0, map.function = "kosambi", condIndex = 3.5 )
#' 
#'  \dontrun{
#'   im1 <- im_scan( fullsib, pheno.col = 1, LOD = TRUE )
#'   summary( im1 )
#'   }
#' 
#'   ## using 4 QTL as cofactors (QTL peaks detected using im_scan)
#'   cofs <- matrix( c( "1", "M2",
#'                      "2", "M27",
#'                      "3", "M37",
#'                      "4", "M52"), 4, 2, byrow = TRUE )
#' 
#'   cof_def <- cof_definition( fullsib, pheno.col = 1, cof.pos = cofs )
#'   cof_def
#'   plot( cof_def )
#' 
#'   \dontrun{
#'   covar <- matrix( rep( c( 1, -1 ), each = 150 ), ncol = 1 )
#'  
#'   ##using 4 QTL as cofactors + 1 unlinked marker (just for illustration)
#'   cofs2 <- matrix( c( "1", "M2",
#'                       "2", "M27",
#'                       "3", "M37",
#'                       "4", "M52",
#'                        NA, "M64" ), 5, 2, byrow = TRUE )
#' 
#'   cof2.def <- cof_definition( fullsib, pheno.col = 2, addcovar = covar, 
#'               cof.pos = cofs2, thres.effect = 1 )
#'   cof2.def
#' 
#'   cof3.def <- cof_definition( fullsib, pheno.col = 1, cof.pos = cofs2, thres.effect = 0.05 )
#'   cof3.def
#'   ### realize with thres.effect = 0.05, the 5th cofactor is removed,
#'   ### because it is non significative (never was selected on cof_selection)
#'   }
#' 

cof_definition <- function(fullsib, pheno.col=1, addcovar=NULL,
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
    
      class(fullsib) <- c("fullsib", "fullsib_cofactors")
  }
  return(fullsib)
}
