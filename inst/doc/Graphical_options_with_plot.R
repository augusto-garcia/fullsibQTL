## ----knitr_init, echo=FALSE, cache=FALSE---------------------------------
library(knitr)

## From ggplot2 rmd vignette options
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>", 
                      fig.width = 6, 
                      fig.height = 6, 
                      fig.align = "center",
                      dev = 'png',
                      dpi = 36)

## ---- eval=TRUE, echo=FALSE,results='hide', cache=FALSE------------------
  library(onemap)

## ---- eval=TRUE, echo=TRUE, results='hide'-------------------------------
  library( "fullsibQTL" )
  data( "example_QTLfullsib" )
  fullsib <- create_fullsib( example_QTLfullsib,
                             list( LG1_final, LG2_final, LG3_final, LG4_final ),
                             step = 1,
                             map.function = "kosambi",
                             condIndex = 3.5 )

  ## IM
  im <- im_scan( fullsib )
  qtls.im <- r2_ls( fullsib, 
                    pheno.col = 1, 
                    lg = c( 1, 1, 2, 2, 2, 3, 3, 4 ),
                    pos = c( "loc13", "M12", "loc45", "loc77", "loc108", "M33", "loc65", "M52" ) )
  
  ## CIM
  cofs.fs <- cof_selection( fullsib, 
                            pheno.col = 1, 
                            k=log(300), 
                            n.cofactor = 10)
  cim <- cim_scan( cofs.fs )
  qtls.cim <- r2_ls( fullsib, 
                     pheno.col = 1, 
                     lg = c( 1, 1, 2, 2, 2, 3, 3, 4 ),
                     pos = c( "loc13", "M12", "loc45", "loc77", "loc108", "M33", "loc65", "M52" ) )
  
  ## Getting im and cim peaks for alpha 0.05.
  im.peak1 <- summary( im_perm, alpha = 0.05 )[ 1, 1 ]
  im.peak2 <- summary( im_perm, alpha = 0.05 )[ 1, 2 ]
  
  cim.peak1 <- summary( cim_perm, alpha = 0.05 )[ 1, 1 ]
  cim.peak2 <- summary( cim_perm, alpha = 0.05 )[ 1, 2 ]

## ------------------------------------------------------------------------
plot( im )

## ------------------------------------------------------------------------
## argument lty
## Changes the line type
## There are 6 different line types as follow:
plot( im )
abline( h = c( 1, 2, 3, 4, 5, 6 ), lty = c ( 1, 2, 3, 4, 5, 6 ) )

## argument lwd
### changes the line width. It changes in a relative proportion to 
### the default. 2 is twice the default width, 3 is triple of the 
### default, and so on.

## You can set any of them with the plot:
plot( im, lty = 3, lwd = 2 )

## argument lg
### selects which LG to plot
## Plotting just LG1
plot( im, lty = 3, lwd = 2, lg = 1 )

## Plotting LG1 and LG2
plot( im, lty=3, lwd = 2, lg = c ( 1 , 2 ) )

## argument label.lg
### labels the LGs
plot( im, lty = 3, lwd = 2, lg = c( 1, 2 ), 
      label.lg = c( "Chr 1", "Chr 2" ) )

## argument cex.axis  
### changes the magnification of axis annotation relative to cex
plot( im, lty = 3, lwd = 2, lg = c( 1, 2 ), 
      label.lg = c( "Chr 1", "Chr 2" ), cex.axis = 1 )
plot( im, lty = 3, lwd = 2, lg = c( 1, 2 ), 
      label.lg = c("Chr 1","Chr 2"), cex.axis = 1.5 )

## argument incl.mkr
### includes the marker position there are 3 options: 
### `points`, `name`, and, `none`. `none` is default.
plot( im, lty = 3, lwd = 2, lg = c( 1, 2 ), 
      label.lg = c( "Chr 1", "Chr 2" ), cex.axis = 1.2, 
      incl.mkr = "points" )

## argument cex.incl
### changes the magnification of incl.mrk annotation 
### relative to cex.
plot( im, lty = 3, lwd = 2, lg = c( 1, 2 ), 
      label.lg = c( "Chr 1", "Chr 2" ), cex.axis = 1.2, 
      incl.mkr = "name", cex.incl = .75 )

## argument gap
### enlarges the gap between chromosomes (in CM)
plot( im, lty = 3, lwd = 2, lg = c( 1, 2 ), 
      label.lg = c( "Chr 1", "Chr 2" ), cex.axis = 1.2, 
      incl.mkr = "name", cex.incl = .75, gap = 50 )

## argument col
### changes the line colors
plot( im, lty = 3, lwd = 2, lg = c( 1, 2 ), 
      label.lg = c( "Chr 1", "Chr 2" ), cex.axis = 1.2, 
      incl.mkr = "name", cex.incl = .75, gap = 50, col="red")

## ------------------------------------------------------------------------
plot( im, lty=1, lwd=2, lg=c(1,2), cex.axis=1.2, 
      incl.mkr="name", cex.incl = .9, gap = 50, col="red", 
      alternate.lgid = TRUE,
      ## Configuring with R plot options:
      main = "IM - Trait - Crop", 
      xlab = "Linkage Group", 
      ylab = "LOD score", 
      yaxt = "n",
      ylim = c(0,10))
axis( side = 2, at = seq( 0, 10, 1 ), lwd = 1 )
abline( h = im.peak1, lty = 2 )
abline( h = im.peak2, lty = 3 )
legend( "topleft", lty = c( 1, 2, 3, 0 ), cex = 0.75, 
        col = c( "red", "black", "black" ), 
        legend = c( "IM", "1st Peak Threshold", "2nd Peak Threshold" ) )

## ------------------------------------------------------------------------
plot( im,cim, lty = c(2,1), lwd = 2, incl.mkr = "name", 
      cex.incl = 0.7, gap = 50, cex.axis = 0.8, 
      col = c("blue","red"),
      main = "IM - Trait - Crop", 
      xlab = "Linkage Group", 
      ylab = "LOD score", 
      yaxt = "n",
      ylim = c( 0, 15 ) )
axis( side = 2, at = seq( 0, 15, 1 ), lwd = 1 )
abline( h = im.peak1, lty = 2, lwd = 0.75 )
abline( h = cim.peak1, lty = 2, lwd = 0.75 )
legend( "topleft", lty = c( 2, 1, 2, 2 ), cex = 0.5, 
        col = c( "red", "blue", "black", "black" ), 
        legend = c( "IM", "CIM", "1st Threshold","2nd Threshold") )

## ------------------------------------------------------------------------
print( cim, lg = 1, pos = c( "loc13", "M12" ) )
print( cim, lg = 2, pos = c( "loc45", "loc77" ,"loc108" ) )
print( cim, lg = 3, pos = c( "M33", "loc65" ) )
print( cim, lg = 4, pos = "M52" )

## ------------------------------------------------------------------------
points.loc <- c( 13, 110, 45, 77, 108, 20, 65, 60 )
points.loc <- points.loc + c( 0, 0,
                              140 + 50, 140 + 50, 140 + 50,
                              140 + 50 + 140 + 50, 140 + 50 + 140 + 50,
                              140 + 50 + 140 + 50 + 140 + 50 )

## ------------------------------------------------------------------------
plot( im, cim, lty = c( 2, 1 ), lwd = 2, 
      incl.mkr = "name", cex.incl = 0.5, gap = 50,
      cex.axis = 0.8, col = c( "blue", "red" ),
      main = "IM - Trait - Crop", 
      xlab = "Linkage Group", 
      ylab = "LOD score", 
      yaxt = "n",
      ylim = c( -2, 15 ) )
axis( side=2, at = seq( 0, 15, 1 ), lwd = 1 )
abline( h = im.peak1, lty = 2, lwd = 0.75 )
abline( h = cim.peak1, lty = 2, lwd = 0.75 )
points( y = rep( 0, 8 ), x = points.loc, pch = 17, col = "dark green" )
legend( "topleft", lty = c( 2, 1, 2, 2, 0 ), cex = 0.75, 
        pch = c( NA, NA, NA, NA, 17 ), 
        col = c( "red", "blue", "black", "black", "dark green" ), 
        legend = c( "IM", "CIM", "1st Threshold", "2nd Threshold", "QTL"),
        bty = 'n', y.intersp = .75 )

## ------------------------------------------------------------------------
sessionInfo()

