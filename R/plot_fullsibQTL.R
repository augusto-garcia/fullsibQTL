#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: plot_fullsibQTL.R                                             #
#                                                                     # 
# Contains:                                                           #
# plot_fullsibQTL                                                     #      
# tidy_fullsibQTL                                                     #
#                                                                     #
# Written by Rodrigo Amadeu                                           #
# rramadeu at gmail dot com                                           #
#                                                                     #
# copyright (c) 2017, Rodrigo Amadeu                                  #
#                                                                     #
# Depends: ggplot2, plotly, htmlwidgets                               #
#                                                                     #
# First version: 06/26/2017                                           #
# Last  version: 08/15/2017                                           #
# License: GPL-3                                                      #
#                                                                     #
#######################################################################

## -------------------------
## plot_fullsibQTL function

#' Multiple Graphics of QTLs Mapping
#' 
#' Given one or more results from 'im_scan' or 'cim_scan' returns a graphic
#' with profiles and map. There is also the possibility to plot QTL positions
#' and effects.
#' 
#' This function first merge the data in a tidy data frame the with
#' \code{\link[fullsibQTL]{tidy_fullsibQTL}}. Then, it makes the respectively
#' plot. The plot is built using \pkg{ggplot2} engine, the interact plot uses
#' \pkg{plotly}. If you are using the interact option with RStudio, it is
#' recommend to plot outside teh RStudio in a browser (e.g. Chrome, Firefox).
#' The returning is a \pkg{ggplot2} object or a \pkg{plotly} object, you can
#' add customizations at them. About how they were made, look at the raw code
#' it is fully commented. You can easily online publish your plot, visit
#' \url{https://plot.ly/r/} to more information. For a complete tutorial, see the vignette:
#' 
#' \code{vignette("Tutorial_fullsibQTL", package = "fullsibQTL")}
#' 
#' The \code{tidy_fullsib} function is an interval (but user available) function to merge and transforms the input data in tidy
#' format.
#' 
#' @aliases plot_fullsibQTL
#' 
#' @param fullsib An object from class \code{fullsib}.
#' @param fullsib.scan An object or a list of objects from class
#' \code{fullsib_scan}.
#' @param r2ls.out An object or a list of objects from class \code{r2ls_out}.
#' @param qtlmapping character vector with the name(s) of \code{fullsib.scan}
#' objects to be ploted.  If \code{qtlmapping=NULL} plot generic names
#' @param thr Numeric vector with threshold values to be plot as lines, similar
#' to \code{abline()}.
#' @param grayscale If \code{TRUE} plots in a gray scale in order to make BW
#' figures.
#' @param lgs Numeric vector with the LGs to be plot. Default is to plots all
#' LGs.
#' @param interact If \code{TRUE} plots a interactive graphic in a \code{html}
#' file.
#' @param file A string. if \code{interact=TRUE}, the \code{html} name to be
#' plot.
#' @param folder A string. If \code{interact=TRUE}, the folder where the created
#' \code{html} will be located. \code{default=getwd()}.
#' @param browser If \code{browser=TRUE}, open the generated \code{html}.
#' @param width numeric argument to be passed for \pkg{plotly}::\code{\link{ggplotly}} internal. Default is
#' (250 x number of LGs) pixels.
#' @param height numeric argument to be passed for \pkg{plotly}::\code{\link{ggplotly}} internal. Default is
#' 1000 pixels.
#' @param lg.colors string vector with the color names for each LG, argument will be passed for 
#' \pkg{ggplot2}::\code{\link{scale_color_manual}}.
#' 
#' @return The resulting graphic is divided per linkage group, it has the map
#' representantion at the bottom, a curve for each one of the \code{fullsib.scan} 
#' arguments. It can also have the QTL representation for each curve. If the user
#' sets interact argument for \code{TRUE}, it will have the features of \pkg{plotly} 
#' package as select, zoom, pan, enable/disable curves, it is also include in the 
#' tooltip information of LOD, QTL estimated effect, position, etc.
#' 
#' @author Rodrigo Amadeu, \email{rramadeu@@gmail.com}
#' 
#' @seealso 
#' \code{\link[fullsibQTL]{create_fullsib}}
#' \code{\link[fullsibQTL]{im_scan}}
#' \code{\link[fullsibQTL]{cim_scan}}
#' \code{\link[fullsibQTL]{r2_ls}}
#' \code{\link[fullsibQTL]{tidy_fullsibQTL}}
#' \code{\link[ggplot2]{ggplot}}
#' \code{\link[plotly]{ggplotly}}
#' 
#' @keywords utilities
#' 
#' @examples
#' 
#'   library( fullsibQTL )
#'   data( "example_QTLfullsib" )
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ),
#'                              step = 1, map.function = "kosambi", condIndex = 3.5 )
#' 
#'   ## IM
#'   im <- im_scan( fullsib )
#'   qtls.im <- r2_ls( fullsib, pheno.col = 1, lg = c( 1, 1, 2, 2, 2, 3, 3, 4 ),
#'                     pos = c( "loc13", "M12", "loc45", "loc77", "loc108", "M33", "loc65", "M52" ) )
#'   ## CIM
#'   cofs.fs <- cof_selection( fullsib, pheno.col = 1, k = log( 300 ), n.cofactor = 10 )
#'   cim <- cim_scan( cofs.fs )
#'   qtls.cim <- r2_ls( fullsib, pheno.col = 1, lg = c( 1, 1, 2, 2, 2, 3, 3, 4 ),
#'                     pos = c( "loc13", "M12", "loc45", "loc77", "loc108", "M33", "loc65", "M52" ) )
#' 
#'   ## Minimal plot
#'   plot_fullsibQTL( fullsib = fullsib,
#'                    fullsib.scan = list( cim, im ) )
#'   
#'   ## QTL plot
#'   plot_fullsibQTL( fullsib = fullsib,
#'                    fullsib.scan = list( cim, im ),
#'                    r2ls.out = list( qtls.cim, qtls.im ) )
#'   
#'   ## Customizing
#'   plot_fullsibQTL( fullsib = fullsib,
#'                    fullsib.scan = list( cim, im ),
#'                    r2ls.out = list( qtls.cim, qtls.im ),
#'                    qtlmapping = c( "cim", "im" ) )
#'                   
#'   plot_fullsibQTL( fullsib = fullsib,
#'                    fullsib.scan = list( cim, im ),
#'                    r2ls.out = list( qtls.cim, qtls.im ),
#'                    qtlmapping = c( "cim", "im" ),
#'                    grayscale = TRUE,
#'                    thr = c( 3.5, 4 ) )
#'                   
#'   plot_fullsibQTL( fullsib = fullsib,
#'                    fullsib.scan = list( cim, im ),
#'                    r2ls.out = list( qtls.cim, qtls.im ),
#'                    qtlmapping = c( "cim", "im" ),
#'                    interact = TRUE )
#' 

plot_fullsibQTL <- function(fullsib = NULL, fullsib.scan = NULL,
                            r2ls.out = NULL, qtlmapping = NULL, lgs=NULL,
                            thr=NULL, grayscale=FALSE,
                            interact = FALSE, file = "fullsibQTL.html", folder=getwd(), browser=TRUE,
                            height = NULL, width = NULL, lg.colors=NULL){
  
  ## Defining global variable to avoid conflicts during R CMD CHECK: ggplot2 issue
  pos.cM = LOD = dummy = loc = r2.qtl = NULL
  
  if(is.null(fullsib))
    stop(deparse("fullsib object must be present"))
  
  if(is.null(fullsib.scan))
    stop(deparse("at list one fullsib.scan object must be present"))
  
  if(any(class(fullsib.scan) == "fullsib_scan"))
    fullsib.scan <- list(fullsib.scan)
  
  if(!is.null(r2ls.out))
    if(any(class(r2ls.out) == "r2ls_out"))
      r2ls.out <- list(r2ls.out)
  
  df <- tidy_fullsibQTL(fullsib, fullsib.scan, r2ls.out, qtlmapping)
  if(!is.null(lgs)){
    df <- df[df$lg %in% lgs, ]
  }else{
    lgs <- unique(df$lg)
  }
  df <- unique(df)
  
  if(!is.null(thr))
    if(class(thr) != "numeric" | !is.vector(thr))
      stop(deparse("thr object must be a numeric vector"))
    
  ## Setting y variable for QTL position
  df$dummy <- scales::rescale(as.numeric(df$qtlmapping),to=c(-0.05,-0.95))

  ## Plotting
  suppressWarnings( ## aes error 'Ignoring unknown aesthetics', used for plotly tooltip graphic
    p <- ggplot() +
      
      ## Add LOD curves 
      geom_line(data=df[df$plot=="lod",],
                aes(x=pos.cM,y=LOD,color=qtlmapping)) +
      
      ## Add QTL position
      geom_point(data=df[(df$r2.qtl > 0 & df$plot=="lod"),], 
                 aes(x=pos.cM,y=dummy,color=qtlmapping,
                     label1=loc,label2=r2.qtl,label3=LOD),
                 shape=17,size=2)+ ## Point char
      
      ## Add Map
      geom_line(data=df[df$plot=="map",],
                aes(x=pos.cM,y=-1)) +
      
      geom_point(data=df[df$plot=="map",],
                 aes(x=pos.cM,y=-1,label1=loc),
                 shape=3,size=1) +
      
      ## Facetting
      facet_grid(.~ lg, scales = "free_x")+
      
      ## Labs
      labs(x = "cM",
           y = "LOD Score",
           title = "Linkage Group") +
      theme(strip.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5)))
    
    ## Setting theme
    if(grayscale){
      p <- p + scale_color_grey(" ") + theme_classic()
      #print("grayscale")
    }else{
      if(is.null(lg.colors)){
      p <- p + scale_color_discrete(" ")
      }else{
        p <- p + scale_color_manual(" ",values = lg.colors)
      }
    }
    
    ## Adding intercepts
    if(!is.null(thr))
      for(i in 1:length(thr))
        p <- p + geom_hline(yintercept=thr[i],linetype=i+1,color="dark gray")
  
  if(interact){
    if(is.null(height))
      height = 500
  
    if(is.null(width))
      width = length(lgs)*250
    
    p <- plotly::ggplotly(p, autosize=T,
                          tooltip = c("LOD", "pos.cM", "loc", "r2.qtl", "qtlmapping"), 
                          width = width, height = height)
    
    htmlwidgets::saveWidget(p, file=file.path(folder, file))
    if(browser){
      openHTML <- function(x,y) browseURL(paste0('file://', file.path(y, x)))
      openHTML(file,folder)
    }
  }else{
    return(p)
  }
}

## -------------------------
## tidy_fullsibQTL function

#' Merges fullsibQTL data in a single and tidy data table
#' 
#' Merges in a single and tidy \code{data.table} objects from one or more of 
#' following classes: \code{fullsib}, \code{fullsib_scan}, and \code{r2ls_out}.
#' This function was developed to be used as part of the 
#' \code{\link[fullsibQTL]{plot_fullsibQTL}}.
#' 
#' @aliases tidy_fullsibQTL
#' 
#' @param fullsib An object from class \code{fullsib}.
#' @param fullsib.scan A list of object(s) from class \code{fullsib_scan}.
#' @param r2ls.out A list of object(s) from class \code{r2ls_out}.
#' @param qtlmapping character vector with the name(s) of \code{fullsib.scan}
#' objects to be ploted.  If \code{qtlmapping=NULL} plot generic names.
#' 
#' @return A merged \code{data.table} with the input data
#' 
#' @author Rodrigo Amadeu, \email{rramadeu@@gmail.com}
#' 
#' @seealso
#' \code{\link[fullsibQTL]{create_fullsib}}
#' \code{\link[fullsibQTL]{im_scan}}
#' \code{\link[fullsibQTL]{cim_scan}}
#' \code{\link[fullsibQTL]{r2_ls}}
#' \code{\link[fullsibQTL]{tidy_fullsibQTL}}
#' \code{\link[fullsibQTL]{plot_fullsibQTL}}
#' \code{\link[ggplot2]{ggplot}}
#' 
#' @keywords utilities
#' 
#' @examples
#' 
#'   library(fullsibQTL)
#'   data("example_QTLfullsib")
#'   fullsib <- create_fullsib(example_QTLfullsib,
#'                              list(LG1_final, LG2_final, LG3_final, LG4_final),
#'                              step=1,map.function="kosambi",condIndex=3.5)
#' 
#'   ## IM
#'   im <- im_scan(fullsib)
#'   qtls.im <- r2_ls(fullsib, pheno.col=1, lg=c(1,1,2,2,2,3,3,4),
#'                  pos=c("loc13","M12","loc45","loc77","loc108","M33","loc65","M52"))
#'   ## CIM
#'   cofs.fs <- cof_selection(fullsib, pheno.col=1, k=log(300), n.cofactor = 10)
#'   cim<-cim_scan(cofs.fs)
#'   qtls.cim <- r2_ls(fullsib, pheno.col=1, lg=c(1,1,2,2,2,3,3,4),
#'                     pos=c("loc13","M12","loc45","loc77","loc108","M33","loc65","M52"))
#' 
#'   ## Merging and transforming to tidy data
#'   tidy_fullsibQTL(fullsib, fullsib.scan=list(im,cim),
#'                   r2ls.out= list(qtls.im,qtls.cim), qtlmapping=c("im","cim"))
#' 

tidy_fullsibQTL <- function(fullsib, fullsib.scan, r2ls.out=NULL, qtlmapping=NULL){
  
  ## Checking fullsib.scan data input
  if(class(fullsib.scan)!="list")
    stop(deparse("fullsib.scan object must be a list"))
  if(any(unlist(lapply(lapply(fullsib.scan,class),head,n=1))!="fullsib_scan"))
    stop(deparse("All fullsib.scan objects in the list must be from 'fullsib_scan' class"))
  
  n <- length(fullsib.scan)
  
  ## Checking qtlmapping
  if(is.null(qtlmapping)){
    qtlmapping <- paste0("qtlmap",1:n)
  }else{
    if(class(qtlmapping)!="character" || !is.vector(qtlmapping))
      stop(deparse("qtlmapping object must be a character vector"))
  }
  
  ## Checking data length
  if(length(fullsib.scan) != length(qtlmapping))
    stop(deparse("all objects must have the same length"))
  
  ## Checking r2ls.out data input
  if(!is.null(r2ls.out)){
    if(class(r2ls.out)!="list")
      stop(deparse("r2ls.out object must be a list"))
    if(any(unlist(lapply(lapply(r2ls.out,class),tail,n=1))!="r2ls_out"))
      stop(deparse("All r2ls.out objects in the list must be from 'r2ls_out' class"))
    if(length(fullsib.scan) != length(r2ls.out))
      stop(deparse("all objects must have the same length"))
  }
  
  ## Building data.frame for plot
  for(i in 1:n){
    markers_temp <- rownames(fullsib.scan[[i]])
    rownames(fullsib.scan[[i]]) <- NULL
    fullsib.scan[[i]] <- data.frame(fullsib.scan[[i]],
                                    loc=markers_temp,
                                    qtlmapping=qtlmapping[i])
    if(!is.null(r2ls.out)){
      r2ls.out[[i]] <- data.frame(r2ls.out[[i]],
                                  qtlmapping=qtlmapping[i])
    }
  }
  
  if(n>1){
    for(i in 2:n){
      if(i==2){
        temp.fs <- rbind(fullsib.scan[[i-1]],fullsib.scan[[i]])
        if(!is.null(r2ls.out))
          temp.r2ls <- rbind(r2ls.out[[i-1]][-1,],r2ls.out[[i]][-1,])
      }else{
        temp.fs <- rbind(temp.fs,fullsib.scan[[i]])
        if(!is.null(r2ls.out))
          temp.r2ls <- rbind(temp.r2ls,r2ls.out[[i]][-1,])
      }
    }
  }else{
    temp.fs <- fullsib.scan[[1]]
    if(!is.null(r2ls.out))
      temp.r2ls <- r2ls.out[[1]][-1,]
  }
  
  df.out <- cbind(index=1:nrow(temp.fs),temp.fs,r2.qtl=0)
  if(!is.null(r2ls.out))
    df.out$r2.qtl[merge(df.out,temp.r2ls,by=c("lg","loc","qtlmapping"))$index] <- temp.r2ls$r2
  df.out$plot <- !is.na(match(as.character(df.out$loc),fullsib$mapped.mkrs))
  class(df.out) <- c("data.frame","tidy.fullsibQTL")
  df.out <- df.out[order(df.out$qtlmapping,df.out$lg,df.out$pos.cM),]
  
  temp <- df.out[df.out$plot,]
  temp$plot <- "map"
  temp$LOD <- -1
  temp$qtlmapping <- "map"
  df.out$plot <- "lod"
  df.out <- rbind(df.out,temp)
  return(df.out)
}
