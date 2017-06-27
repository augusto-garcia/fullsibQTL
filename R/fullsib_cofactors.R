#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: fullsib_cofactors.R                                           #
#                                                                     # 
# Contains:                                                           #
# plot.fullsib_cofactors                                              #
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
# Function: plot.fullsib_cofactors                                    #
#                                                                     #
# This function is based on draw_map presented in onemap pkg. It is   #
# used to plot genetic map with the cofactors on map                  #
#######################################################################

## -------------------------
## plot.fullsib_cofactors function

#' Plot \code{fullsib_cofactors} object
#' 
#' This function is based on \code{\link[onemap]{draw_map}} presented in \pkg{onemap}. It is
#' used to plot genetic map with the cofactors on map    
#' 
#' @aliases plot.cof_selection
#' 
#' @param x An object from class \emph{fullsib_cofactors}.
#' @param horizontal if \code{TRUE}, indicates that the map should be plotted
#' horizontally. Default is \code{FALSE}.
#' @param grid if \code{TRUE}, displays a grid in the background. Default is
#' \code{FALSE}.
#' @param cex.mrk the magnification to be used for markers.
#' @param cex.grp the magnification to be used for group axis annotation.
#' @param ... Further arguments, passed to other methods. Currently ignored.
#' 
#' @return A plot of a map with the \emph{fullsib_cofactors} positions.
#' 
#' @seealso 
#' \code{\link[onemap]{draw_map}},
#' \code{\link[fullsibQTL]{cof_selection}},
#' \code{\link[fullsibQTL]{cof_definition}}
#' 
#' @keywords utilities
#' @examples
#'   data(example_QTLfullsib)
#' 
#'   fullsib <- create_fullsib(example_QTLfullsib,
#'                             list(LG1_final, LG2_final, LG3_final, LG4_final),
#'                             step=0,map.function="kosambi",condIndex=3.5)
#' 
#' 
#' 
#'   ###############################################
#'   ## cofactor selection using BIC (n.ind = 300)
#' 
#'   ### just using markers that are placed on linkage groups (default)
#'   cofs.fs <- cof_selection(fullsib, pheno.col=1, k = log(300), selection =1)
#'   cofs.fs
#'   plot(cofs.fs)
#'   plot(cofs.fs, horizontal = TRUE, grid = TRUE, cex.mrk = 2.5, cex.grp = 1)
#'   

plot.fullsib_cofactors  <- function(x, horizontal = FALSE, grid = FALSE,
                                   cex.mrk = 1, cex.grp = 0.75, ...)
{

  
  ##draw genetic map, as in onemap
  draw_map2(x$map, horizontal, grid, cex.mrk, cex.grp, trait=x$cofactors$trait.cof)
  
  ##find the positions in the genome
  mkr.pos <- rep(NA, length(x$cofactors$names.cof$mkr))
  
  for (i in 1:length(x$cofactors$names.cof$mkr)){
    lg.i <- x$cofactors$names.cof$lg[i]
    mkr.i <- x$cofactors$names.cof$mkr[i]
    tmp <- match(mkr.i, names(x$probs[[ lg.i ]]$newmap$dist))
    if(is.na(tmp))
      mkr.pos[i] <- NA
    else
      mkr.pos[i] <- x$probs[[ lg.i ]]$newmap$dist[tmp]
  }
  
  cofactors <- data.frame(LG = x$cofactors$names.cof$lg,
                          mkr.pos = mkr.pos)
  ##plot the cofactors in the map
  if (horizontal == TRUE)
    points(cofactors$mkr.pos, (length(x$map) + 1 - cofactors$LG),
           pch=17, cex=1.3)
  else
    points(-(length(x$map)+ 1 - cofactors$LG), -cofactors$mkr.pos,
           pch=17, cex=1.3)
}#end function
