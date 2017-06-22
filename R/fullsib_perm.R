
#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: fullsib_perm.R                                                #
#                                                                     # 
# Contains:                                                           #
# summary_fullsib_perm                                                #
# plot_fullsib_perm                                                   #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
# summary_fullsib_perm is based on summary.scanoneperm from qtl pkg   #
# plot_fullsib_perm    is based on plot.scanoneperm from qtl package  #
#                                                                     #
# First version: 09/30/2011 (american date format)                    #
# Last  version: 09/30/2011 (american date format)                    #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#######################################################################
## summary_fullsib_perm                                               #
## function that provides the threshold values, for differents alphas #
## by default it show 0.05, 0.10                                      #
## column one is the classical threshold by Churchill and Doerge, 1994#
## second columns and so, are modification by Chen and Storey, 2006   #
#######################################################################

summary.fullsib_perm <- function(object,
                                 alpha=c(0.05, 0.10), verbose=TRUE, ...)
{
  threshold <-  apply(object, 2, quantile, 1 - alpha, na.rm=T)

  if(length(alpha) == 1 && length(threshold) > 1)
    threshold <- t(threshold)
  
  threshold <- as.data.frame(threshold)
  colnames(threshold) <- colnames(object)
  rownames(threshold) <- round(alpha,2)
  
  if(verbose){
    cat(" Threshold considering", nrow(object), "permutations\n")
    if(ncol(object) == 2){
      cat(" First column indicates", sQuote("classical"), "threshold showed by Churchill and Doerge, 1994\n")
      cat(" Second column means threshold values suggested by Chen and Storey, 2006\n")
    }
  }
  threshold
}

#######################################################################
## plot_fullsib_perm                                                  #
## function that plots the distribution of values obtained by the     #
## permutation test. One can detail wich peak should be plot          #
#######################################################################

plot.fullsib_perm <- function(x, peak=1, ...)
{

  x <- x[,peak]
  n.brk <- ceiling(2*sqrt(length(x)))
  xlim <- c(0,max(as.numeric(x), na.rm=TRUE))
  
  dots <- list(...)
  main <- ""
  xlab <- "maximum"
  
  if("breaks" %in% names(dots)) {
    if("main" %in% names(dots)) {
      if("xlab" %in% names(dots)) {
        hist(x, xlim=xlim, ...)
      }
      else {
        hist(x, xlim=xlim, xlab=xlab, ...)
      }
    }
    else {
      if("xlab" %in% names(dots)) {
        hist(x, xlim=xlim, main=main, ...)
      }
      else {
        hist(x, xlim=xlim, xlab=xlab, main=main, ...)
      }
    }
  }
  else {
    if("main" %in% names(dots)) {
      if("xlab" %in% names(dots)) {
        hist(x, breaks=n.brk, xlim=xlim, ...)
      }
      else {
        hist(x, breaks=n.brk, xlim=xlim, xlab=xlab, ...)
      }
    }
    else {
      if("xlab" %in% names(dots)) {
        hist(x, breaks=n.brk, xlim=xlim, main=main, ...)
      }
      else {
          hist(x, breaks=n.brk, xlim=xlim, xlab=xlab, main=main, ...)
        }
    }
  }
  
  rug(as.numeric(x))
}
