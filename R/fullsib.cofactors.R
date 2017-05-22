
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
# Function: plot.fullsib.cofactors                                    #
#                                                                     #
# This function is based on draw.map presented in onemap pkg. It is   #
# used to plot genetic map with the cofactors on map                  #
#######################################################################

plot.fullsib.cofactors <- function(x, horizontal = FALSE, grid = FALSE,
                                   cex.mrk = 1, cex.grp = 0.75,...)
{

  
  ##draw genetic map, as in onemap
  draw.map2(x$map, horizontal, grid, cex.mrk, cex.grp, trait=x$cofactors$trait.cof)
  
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
