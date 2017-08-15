#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# Contains:                                                           #
# summary.fullsib_scan                                                #
# plot.fullsib_scan                                                   #
# print.fullsib_scan                                                  #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# Reviewed by Rodrigo Amadeu                                          #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
# summary.fullsib_scan is based on summary.scanone from qtl package   #
# plot.fullsib_scan    is fully based on plot.scanone from qtl pkg    #
#                      basically i added more objects to draw and     #
#                      adapted some arguments                         #
#                                                                     #
# First version: 09/30/2011                                           #
# Last  version: 06/22/2017                                           #
# License: GPL-3                                                      #
#                                                                     #
#######################################################################


#######################################################################
# Function: summary_fullsib_scan                                      #
#                                                                     #
# provides a summary of fullsib_scan object (im_scan or cim_scan)     #
# show the highest peak for each group, thr is a threshold value for  #
# printing lg peak.                                                   #
#                                                                     #
#######################################################################

## -------------------------
## summary.fullsib_scan function

#' Summarize QTL mapping scan
#' 
#' Summarizes the QTL mapping scan resulted from the \code{im_scan} and
#' \code{cim_scan}.
#' 
#' This function was designed to get only the maximum peak of each analysed
#' group. For the cases that more than one peak is present should decided the
#' position of the QTL examinating object of class \emph{fullsib_scan} by
#' oneself. We assumed this, because to define a QTL can be very subjective in
#' some situations.
#' 
#' @param object An object of class \emph{fullsib_scan} (output of
#' \code{im_scan} and \code{cim_scan}.
#' @param thr threshold value to declare a QTL, i.e., only prints peaks with
#' value higher than \sQuote{thr}.
#' @param \dots parameters to be passed.
#' 
#' @return It gets the object of class \emph{fullsib_scan} and simplifies in a
#' matrix containing just the maximum peak of each linkage group that exceeds
#' the threshold value defined by the user.
#' 
#' @author Rodrigo Gazaffi \email{rgazaffi@@gmail.com}
#' 
#' @seealso 
#' \code{\link[fullsibQTL]{im_scan}}
#' \code{\link[fullsibQTL]{cim_scan}}
#' \code{\link[fullsibQTL]{im_char}}
#' \code{\link[fullsibQTL]{cim_char}}
#' \code{\link[fullsibQTL]{plot.fullsib_scan}}
#' \code{\link[fullsibQTL]{summary.fullsib_scan}}
#' \code{scanone} from \pkg{qtl} package
#' \code{summary.scanone} from \pkg{qtl} package
#' 
#' @keywords utilities
#' @examples
#' 
#'   data("example_QTLfullsib")
#' 
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ), 
#'                              step = 0, map.function = "kosambi", condIndex = 3.5) 
#' 
#'   im1 <- im_scan( fullsib, lg = "all", pheno.col = 1, LOD = TRUE )
#'   summary( im1 )
#'   summary( im1, thr = 6.5 )
#' 
#' 

summary.fullsib_scan <- function(object, thr=0, ...)
{
  if(!any(class(object) == "fullsib_scan"))
    stop("Input should have class 'fullsib_scan'")

  if(!is.numeric(thr)){
    warning("thr argument need to be a number, equal or higher than 0\n\n")
    thr <- 0
  }
  
  genome.peaks <- NULL
  tmp2 <- NULL

  lg.unique <- unique(object[,1])
  for (i in lg.unique){
    subobject <- subset(object,object[,1]==i)
    max.peak <- which.max(subobject[,3])    
    tmp1 <- subobject[max.peak,]
    tmp2 <- c(tmp2, rownames(subobject)[max.peak])
    genome.peaks <- rbind(genome.peaks,tmp1)
  }
  rownames(genome.peaks) <- tmp2

  if(nrow(genome.peaks) == 1){
    if(genome.peaks[,3] > thr){
      colnames(genome.peaks) <- colnames(object)
      return(genome.peaks)
    }
    else
      cat(" There were no QTL peaks above the threshold.\n")
  }
  else{#more then one peak
    up2thres <- which(genome.peaks[,3] > thr)
    if(length(up2thres) == 0)
      cat(" There were no QTL peaks above the threshold.\n")
    else{ #aqui se tiver um pico consertar
      genome.peaks <- genome.peaks[up2thres,]
      if(is.null(dim(genome.peaks)))
        dim(genome.peaks) <- c(1,4)
      rownames(genome.peaks) <- tmp2[up2thres]
      colnames(genome.peaks) <- colnames(object)
        return(genome.peaks)
      }
    }
}

#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: plot.fullsib_scan.R                                           #
#                                                                     # 
# Contains: plot.fullsib_scan                                         #
#           matchchr                                                  #
#                                                                     #
# Written by Karl Broman,                                             #
# Adapt by Rodrigo Gazaffi: just managed the object class, once the   #
# output structure were fully based on scanone structure              #
# copyright (c) 2008, Karl Broman                                     #
#                                                                     #
# I needed to rename chr to lg to keep the names along the functions. #
# Adapted version: 09/30/2011 (american date format)                  #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

## -------------------------
## plot.fullsib_scan function

#' Plot the QTL mapping scan
#' 
#' Plot the QTL mapping scan resulted from the \code{im_scan} and
#' \code{cim_scan}. This is fully based on \code{plot.scanone} from \pkg{qtl}
#' package.
#' 
#' This function is basically \code{plot.scanone} written by Karl Broman and
#' adapted for this package. We also adapt the functions to receive more object
#' to plot QTL profile.
#' 
#' @param x An object of class \emph{fullsib_scan}, output of \code{im_scan}
#' and \code{cim_scan}.
#' @param x2 Second object of class \emph{fullsib_scan}, optional.
#' @param x3 Third object of class \emph{fullsib_scan}, optional.
#' @param x4 Fourth object of class \emph{fullsib_scan}, optional.
#' @param x5 Fifth object of class \emph{fullsib_scan}, optional.
#' @param lg Indicates which linkage groups will be plotted, default is to plot
#' all groups.
#' @param label.lg for more than one linkage group plotted, a label for each
#' group can be easily added to the plot, the length must be the same of the
#' number of LG.
#' @param cex.axis The same option for argument used in par function for
#' graphical parameters. See \sQuote{?par} for details.
#' @param incl.mkr indicates how the markers will be displayed on the plot. It
#' can assume three options: \code{points}, \code{name}, \code{none}, in which
#' the first indicates the position of the markers with small triangules, the
#' second indicates by its name (sometimes it is useful if only one group is
#' plotted, otherwise it can be difficult to distinguish the names) and the
#' last none additional information is given. Here, the default is to consider
#' \code{incl.mkr = points}.
#' @param cex.incl If incl.mkr is indicated as \code{points} or \code{name}, an
#' integer can be defined to resize the symbols (for points) or text (for
#' names) displayed at the graphic.
#' @param xlim Limits for x-axis (optional).
#' @param ylim Limits for y-axis (optional).
#' @param lty Line types; vector of length 1 until 5, according to the number
#' of object of class \emph{fullsib_scan} plotted.
#' @param lwd Line widths; vector of length 1 until 5, according to the number
#' of object of class \emph{fullsib_scan} plotted.
#' @param add If TRUE, add to a current plot
#' @param gap Gap separating chromosomes (in cM).
#' @param col Line colors; vector of length 1 until 5, according to the number
#' of object of class \emph{fullsib_scan} plotted.
#' @param alternate.lgid If TRUE and more than one linkage group is plotted,
#' alternate the placement of chromosome axis labels, so that they may be more
#' easily distinguished.
#' @param \dots Passed to the function \code{plot} when it is called.
#' 
#' @return None
#' 
#' @author Mainly writted by Karl Broman and present at \code{plot.scanone} and
#' adapted by Rodrigo Gazaffi \email{rgazaffi@@gmail.com}
#' 
#' @seealso 
#' \code{\link[fullsibQTL]{im_scan}}
#' \code{\link[fullsibQTL]{cim_scan}}
#' \code{\link[fullsibQTL]{summary.fullsib_scan}}
#' \code{par}
#' \code{scanone} from \pkg{qtl} package
#' \code{plot.scanone} from \pkg{qtl} package
#' 
#' @keywords plot
#' 
#' @examples
#'   data( "example_QTLfullsib" )
#' 
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ), 
#'                              step = 0, map.function = "kosambi", condIndex = 3.5 ) 
#' 
#' 
#' 
#'   ## running im_scan:
#'   im1 <- im_scan( fullsib, lg = "all", pheno.col = 1, LOD = TRUE )
#'   ## results in LOD Score for all linkage groups
#' 
#'   im2 <- im_scan( fullsib, lg = 1:4, pheno.col = 1, LOD = FALSE )
#'   ## results in -log10(pvalue) for all linkage groups
#' 
#'   plot( im1 ) ##default
#'   ##customizing graphics
#'   plot( im1, lty = 3, lwd = 2, incl.mkr = "points" )
#'   plot( im1, lty = 3, lwd = 2, incl.mkr = "name" )
#'   plot( im1, lty = 3, lwd = 2, incl.mkr = "name", cex.incl = 0.5 )
#'   plot( im1, im2, label.lg = c( "1st LG", "2nd LG","3rd LG","4th LG" ) )
#'   plot( im1, im2, label.lg = c( "1st LG", "2nd LG","3rd LG","4th LG" ),
#'         alternate.lgid=TRUE)
#' 
#'   #add argument is another way of overlapping QTL results
#'   plot( im1, lg = 2, incl.mkr = "name" )
#'   plot( im2, lg = 2, incl.mkr = "name", add = TRUE )
#' 
#' 

plot.fullsib_scan <- function(x, x2, x3, x4, x5, lg, label.lg, cex.axis=1,
                              incl.mkr = c("none", "points", "name"), cex.incl,
                              xlim, ylim, lty = 1, lwd = 2, add = FALSE, gap = 50,
                              col = c("black", "blue", "red", "orange", "limegreen"),
                              alternate.lgid = FALSE, ...)
  {

  if (!any(class(x) == "fullsib_scan") ||
      (!missing(x2) && !any(class(x2) == "fullsib_scan")) ||
      (!missing(x3) && !any(class(x3) == "fullsib_scan")) ||
      (!missing(x4) && !any(class(x4) == "fullsib_scan")) ||
      (!missing(x5) && !any(class(x5) == "fullsib_scan")))
    stop("Input should have class 'fullsib_scan'")

  dots <- list(...)
  
  incl.mkr <- match.arg(incl.mkr)
  switch(EXPR = incl.mkr,
         "none" = {
           CEX.INCL <- NA
         },
         "points" = {
           if(missing(cex.incl))
             cex.incl <- 1.25
           else
             cex.incl <- cex.incl
         },
         "name" = {
           if(missing(cex.incl))
             cex.incl <- 0.80
           else
             cex.incl <- cex.incl          
         })        

           
  if (length(dim(x)) != 2) 
    stop("Argument x must be a matrix or data.frame.")
  if (!missing(x2) && length(dim(x2)) != 2) 
    stop("Argument x2 must be a matrix or data.frame.")
  if (!missing(x3) && length(dim(x3)) != 2)
    stop("Argument x3 must be a matrix or data.frame.")
  if (!missing(x4) && length(dim(x4)) != 2)
    stop("Argument x3 must be a matrix or data.frame.")
  if (!missing(x5) && length(dim(x5)) != 2)
    stop("Argument x3 must be a matrix or data.frame.")


  second <- third <- fourth <- fifth <- TRUE
  
  if (missing(x2) && missing(x3) && missing(x4) && missing(x5)) 
    second <- third <- fourth <- fifth <- FALSE
  if (missing(x5))   fifth <- FALSE
  if (missing(x4))   fourth <- FALSE
  if (missing(x3))   third <- FALSE
  if (missing(x2))   second <- FALSE

  
  matchchr <- function(selection, thelg) {
    if (is.factor(thelg)) 
      thelg <- as.character(thelg)
    if (length(thelg) > length(unique(thelg))) 
      stop("Duplicate group names.")
    if (is.logical(selection)) {
      if (length(selection) != length(thelg)) 
        stop("Logical vector to select groups is the wrong length")
      return(thelg[selection])
    }
    if (is.numeric(selection)) 
      selection <- as.character(selection)
    if (length(selection) > length(unique(selection))) {
      warning("Dropping duplicate groups")
      selection <- unique(selection)
    }
    g <- grep("^-", selection)
    if (length(g) > 0 && length(g) < length(selection)) 
      stop("In selecting groups, all must start with '-' or none should.")
        if (length(g) > 0) {
            selectomit <- TRUE
            selection <- substr(selection, 2, nchar(selection))
          }
        else selectomit <- FALSE
    wh <- match(selection, thelg)
    if (any(is.na(wh))) {
      warning("Lg ", paste(selection[is.na(wh)], collapse = ", "), 
              " not found")
      wh <- wh[!is.na(wh)]
      if (length(wh) == 0) 
        return(thelg)
    }
    if (selectomit) 
      return(thelg[-wh])
    #thelg[sort(wh)]
    thelg[wh]
  }

  out <- x[, 1:4]
  if (second) 
    out2 <- x2[, 1:3]
  if (third) 
    out3 <- x3[, 1:3]
  if (fourth) 
    out4 <- x4[, 1:3]
  if (fifth) 
    out5 <- x5[, 1:3]
  
  
  if (length(lty) == 1) 
    lty <- rep(lty, 5)
  if (length(lwd) == 1) 
    lwd <- rep(lwd, 5)
  if (length(col) == 1) 
    col <- rep(col, 5)
  if (missing(lg) || length(lg) == 0) 
    lg <- unique(as.character(out[, 1]))
  else lg <- matchchr(lg, unique(out[, 1]))

  out <- out[!is.na(match(out[, 1], lg)), ]
  if (second) 
    out2 <- out2[!is.na(match(out2[, 1], lg)), ]
  if (third) 
    out3 <- out3[!is.na(match(out3[, 1], lg)), ]
  if (fourth) 
    out4 <- out4[!is.na(match(out4[, 1], lg)), ]
  if (fifth) 
    out5 <- out5[!is.na(match(out5[, 1], lg)), ]
  
  onelg <- FALSE
  if (length(lg) == 1) {
    gap <- 0
    onelg <- TRUE
  }

  temp <- out

  ##here, it was indicated the range of the chrs
  begend <- matrix(unlist(tapply(temp[, 2], temp[, 1], range)), 
                   ncol = 2, byrow = TRUE)
  rownames(begend) <- unique(out[, 1])
  begend <- begend[as.character(lg), , drop = FALSE]

  len <- begend[, 2] - begend[, 1]
  
  if (!onelg) 
    start <- c(0, cumsum(len + gap)) - c(begend[, 1], 0)
  else start <- 0

  maxx <- sum(len + gap) - gap

  if (all(is.na(out[, 3]))) 
    maxy <- 1
  else maxy <- max(out[, 3], na.rm = TRUE)
 
  if (second) 
    maxy <- max(c(maxy, out2[, 3]), na.rm = TRUE)

  if (third) 
    maxy <- max(c(maxy, out3[, 3]), na.rm = TRUE)

  if (fourth) 
    maxy <- max(c(maxy, out4[, 3]), na.rm = TRUE)
  
  if (fifth) 
    maxy <- max(c(maxy, out5[, 3]), na.rm = TRUE)
  

  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd = FALSE, las = 1)
  on.exit(par(xpd = old.xpd, las = old.las))
  if (missing(ylim)) 
    ylim <- c(0, maxy)
  if (missing(xlim)) {
    if (onelg) 
      xlim <- c(0, max(out[, 2]))
    else xlim <- c(-gap/2, maxx + gap/2)
  }

  if (!add) {
    ##for overlap qtl curves in thge same graphic if TRUE...
    if (onelg) {
      if ("ylab" %in% names(dots)) {
        if ("xlab" %in% names(dots)) {
          plot(0, 0, ylim = ylim, xlim = xlim, type = "n", cex.axis=cex.axis,
               ...)
        }
        else {
          plot(0, 0, ylim = ylim, xlim = xlim, type = "n",  cex.axis=cex.axis,
               xlab = "Map position (cM)", ...)
        }
      }
      else {
        if ("xlab" %in% names(dots)) {
          plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
               ylab = dimnames(out)[[2]][3], cex.axis=cex.axis, ...)
        }
        else {
          plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
               xlab = "Map position (cM)", ylab = dimnames(out)[[2]][3], cex.axis=cex.axis,
               ...)
        }
      }
    }#onelg
    else {
      if ("ylab" %in% names(dots)) {
        if ("xlab" %in% names(dots)) {
          plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
               xaxt = "n", xaxs = "i", cex.axis=cex.axis, ...)
        }
        else {
          plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
               xaxt = "n", xlab = "Group", xaxs = "i",  cex.axis=cex.axis,...)
        }
      }
      else {
        if ("xlab" %in% names(dots)) {
          plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
               xaxt = "n", ylab = dimnames(out)[[2]][3], 
               xaxs = "i",  cex.axis=cex.axis,...)
        }
        else {
          plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
               xaxt = "n", xlab = "Group", ylab = dimnames(out)[[2]][3], 
               xaxs = "i",  cex.axis=cex.axis, ...)
        }
      }
    }
  }

  xtick <- NULL
  xticklabel <- NULL
  for (i in 1:length(lg)) {
    ## here the trick, for plotting many chrs.
    ## X is only one big vector for 1 to k cM
    ## so key point identify were in this vector the chr start and end, and draw it
    x <- out[out[, 1] == lg[i], 2] + start[i] ##axis x
    y <- out[out[, 1] == lg[i], 3] ##axis y
    if (length(x) == 1) {
      g <- max(gap/10, 2)
      x <- c(x - g, x, x + g)
      y <- rep(y, 3)
    }

    ##lines are draw...
    lines(x, y, lwd=lwd[1], lty=lty[1], col=col[1])

    if (!add && !onelg) {
      tloc <- mean(c(min(x), max(x))) ##which is the middle point of the chr...
      xtick <- c(xtick, tloc)

      if (missing(label.lg))
        xticklabel <- c(xticklabel, as.character(lg[i]))
      else xticklabel <- c(xticklabel, label.lg[i])
    }
    if (second) {
      x <- out2[out2[, 1] == lg[i], 2] + start[i]
      y <- out2[out2[, 1] == lg[i], 3]
      if (length(x) == 1) {
        g <- max(gap/10, 2)
        x <- c(x - g, x, x + g)
        y <- rep(y, 3)
      }
      lines(x, y, lty = lty[2], col = col[2], lwd = lwd[2])
    }
    if (third) {
      x <- out3[out3[, 1] == lg[i], 2] + start[i]
      y <- out3[out3[, 1] == lg[i], 3]
      if (length(x) == 1) {
        g <- max(gap/10, 2)
        x <- c(x - g, x, x + g)
        y <- rep(y, 3)
      }
      lines(x, y, lty = lty[3], col = col[3], lwd = lwd[3])
    }

    if (fourth) {
      x <- out4[out4[, 1] == lg[i], 2] + start[i]
      y <- out4[out4[, 1] == lg[i], 3]
      if (length(x) == 1) {
        g <- max(gap/10, 2)
        x <- c(x - g, x, x + g)
        y <- rep(y, 3)
      }
      lines(x, y, lty = lty[4], col = col[4], lwd = lwd[3])
    }
    if (fifth) {
      x <- out5[out5[, 1] == lg[i], 2] + start[i]
      y <- out5[out5[, 1] == lg[i], 3]
      if (length(x) == 1) {
        g <- max(gap/10, 2)
        x <- c(x - g, x, x + g)
        y <- rep(y, 3)
      }
      lines(x, y, lty = lty[5], col = col[5], lwd = lwd[3])
    }

    if (!add) {
      nam <- dimnames(out)[[1]][out[, 1] == lg[i]]
      wh.genoprob <- grep("^c.+\\.loc-*[0-9]+", nam)
      if (length(wh.genoprob) == 0) 
        wh.genoprob <- seq(along = nam)
      else wh.genoprob <- (seq(along = nam))[-wh.genoprob]

      if (incl.mkr == "none"){
        pos <- out[out[, 1] == lg[i], 2][wh.genoprob] + start[i]
        rug(pos, 0.02, quiet = TRUE) ##the "map" on the graphic
      }
      else if (incl.mkr == "points"){
        pos <- out[out[, 1] == lg[i], 2][wh.genoprob] + start[i]
        rug(pos, 0.02, quiet = TRUE) ##the "map" on the graphic
        pos.mkr <- which(substring(names(pos), 1,3) != "loc")
        points(pos[pos.mkr], rep(0,length(pos.mkr)), pch = 17, cex = cex.incl)
      }
      else if (incl.mkr == "name"){
        pos <- out[out[, 1] == lg[i], 2][wh.genoprob] + start[i]
        rug(pos, 0.02, quiet = TRUE) ##the "map" on the graphic
        a <- par("usr")
        mkr.names <- nam[wh.genoprob]
        mkr.names[which(substring(mkr.names, 1,3) == "loc")] <- NA
        text(pos, rep(a[3] + diff(a[3:4]) * 0.03, length(pos)), 
             mkr.names, srt = 90, adj = c(0, 0.5),cex=cex.incl)
      }
    }
  }
  
  if (!add && !onelg) {
    if (!alternate.lgid || length(xtick) < 2) {
      for (i in seq(along = xtick)) {
        axis(side = 1, at = xtick[i], labels = xticklabel[i], cex.axis=cex.axis)      
      }
    }
    else {
      odd <- seq(1, length(xtick), by = 2)
      even <- seq(2, length(xtick), by = 2)
      for (i in odd) {
        axis(side = 1, at = xtick[i], labels = "")
        axis(side = 1, at = xtick[i], labels = xticklabel[i], 
             line = -0.4, tick = FALSE, cex.axis=cex.axis)
        
      }
      for (i in even) {
        axis(side = 1, at = xtick[i], labels = "")      
        axis(side = 1, at = xtick[i], labels = xticklabel[i], 
             line = +0.4, tick = FALSE, cex.axis=cex.axis)
      }
      ##print(xticklabel)
    }
  } 
  invisible()
}

#######################################################################
# Function: print.fullsib_scan                                        #
#                                                                     #
# prints the object of class fullsib_scan (im_scan or cim_scan)       #
# if lg is defined it prints just the intended linkage group          #
# if pos is defined it prints just the markers                        #
# default is to print all linkage groups and all markers              #
#######################################################################

## -------------------------
## print.fullsib_scan function

#' Prints object of class \emph{fullsib_scan}
#' 
#' Prints the object of class \emph{fullsib_scan} created by \code{im_scan} or
#' \code{cim_scan}.
#' 
#' @param x Object of class \emph{fullsib_scan} resulted from \code{im_scan} or
#' \code{cim_scan}.
#' @param lg Linkage group that will be printed. Default is to print all
#' scanned groups.
#' @param pos Position of the QTL, string(s).
#' @param \dots parameters to be passed.
#' 
#' @return None
#' 
#' @author Rodrigo Gazaffi, \email{rgazaffi@@gmail.com}
#' 
#' @seealso 
#' \code{\link[fullsibQTL]{im_scan}}
#' \code{\link[fullsibQTL]{cim_scan}}
#' \code{\link[fullsibQTL]{plot.fullsib_scan}}
#' \code{\link[fullsibQTL]{summary.fullsib_scan}}
#' \code{scanone} from \pkg{qtl} package
#' 
#' @keywords utilities
#' 
#' @examples
#'   data( "example_QTLfullsib" )
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ), 
#'                              step = 0, map.function = "kosambi",
#'                              condIndex = 3.5)
#' 
#'   ## running im_scan:
#'   im1 <- im_scan( fullsib, lg = "all", pheno.col = 1, LOD = TRUE )
#' 
#'   im1
#'   print( im1, lg = 2 )
#'   print( im1, lg = 2, pos = "M27" )
#'   print( im1, lg = c( 1, 4 ) )
#' 

print.fullsib_scan <- function(x, lg, pos, ...){

  if(missing(lg))
    lgs <- unique(x[,"lg"])
  else{
    lgs <- match(lg, unique(x[,"lg"]))
    check.na <- which(is.na(lgs))
    if(length(check.na) != 0){
      lgs <- lgs[-check.na]
    }
    lgs <- unique(x[,"lg"])[lgs]
  }

  if(length(lgs) == 0)
    stop("lg should be integer between ", min(unique(x[,"lg"])), " and ", max(unique(x[,"lg"])))
  else{
    lg.print <- NULL
    lg.print <- function(y) {x[,"lg"] == y}
    if(!missing(pos)){
      temp <- sapply(lgs, lg.print)*(!is.na(match(row.names(x),pos)))
    }else{
      temp <- sapply(lgs, lg.print)
    }
    lines2print <- apply(temp,1, any)
    ##lines2print <- apply(sapply(lgs, function(x) {im.lod[,"lg"] == x}),1, any)   \
    temp <- subset(x, lines2print)
    return(print(temp))
  }
}
