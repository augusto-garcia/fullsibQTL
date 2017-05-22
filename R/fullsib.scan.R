#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# Contains:                                                           #
# summary.fullsib.scan                                                #
# plot.fullsib.scan                                                   #
# print.fullsib.scan                                                  #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
# summary.fullsib.scan is based on summary.scanone from qtl package   #
# plot.fullsib.scan    is fully based on plot.scanone from qtl pkg    #
#                      basically i added more objects to draw and     #
#                      adapted some arguments                         #
#                                                                     #
# First version: 09/30/2011                                           #
# Last  version: 09/30/2011                                           #
# License: GPL-3                                                      #
#                                                                     #
#######################################################################


#######################################################################
# Function: summary.fullsib.scan                                      #
#                                                                     #
# provides a summary of fullsib.scan object (im.scan or cim.scan)     #
# show the highest peak for each group, thr is a threshold value for  #
# printing lg peak.                                                   #
#                                                                     #
#######################################################################

summary.fullsib.scan <- function(object, thr=0,...)
{
  if(!any(class(object) == "fullsib.scan"))
    stop("Input should have class 'fullsib.scan'")

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
# File: plot.fullsib.scan.R                                           #
#                                                                     # 
# Contains: plot.fullsib.scan                                         #
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

plot.fullsib.scan <- function(x, x2, x3, x4, x5, lg, label.lg, cex.axis=1,
                              incl.mkr = c("none", "points", "name"), cex.incl,
                              xlim, ylim, lty = 1, lwd = 2, add = FALSE, gap = 50,
                              col = c("black", "blue", "red", "orange", "limegreen"),
                              alternate.lgid = FALSE, ...)
  {

  if (!any(class(x) == "fullsib.scan") ||
      (!missing(x2) && !any(class(x2) == "fullsib.scan")) ||
      (!missing(x3) && !any(class(x3) == "fullsib.scan")) ||
      (!missing(x4) && !any(class(x4) == "fullsib.scan")) ||
      (!missing(x5) && !any(class(x5) == "fullsib.scan")))
    stop("Input should have class 'fullsib.scan'")

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
# Function: print.fullsib.scan                                        #
#                                                                     #
# prints the object of class fullsib.scan (im.scan or cim.scan)       #
# if lg is defined it prints just the intended linkage group          #
# default is to print all linkage groups                              #
#######################################################################

print.fullsib.scan <- function(x, lg, ...){

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
    temp <- sapply(lgs, lg.print)
    lines2print <- apply(temp,1, any)
    ##lines2print <- apply(sapply(lgs, function(x) {im.lod[,"lg"] == x}),1, any)
    return(print(subset(x, lines2print)))
  }
} 

