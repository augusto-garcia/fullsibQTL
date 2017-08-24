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

## ---- eval=FALSE---------------------------------------------------------
#  vignette( "QTL_mapping_with_partially_informative_markers", package = "fullsibQTL" )

## ---- eval=FALSE---------------------------------------------------------
#  vignette( "Graphical_options_with_plot_fullsibQTL", package = "fullsibQTL" )

## ---- eval=FALSE---------------------------------------------------------
#  vignette( "Graphical_options_with_plot", package = "fullsibQTL" )

## ---- eval=FALSE---------------------------------------------------------
#  install.packages( "fullsibQTL" )

## ---- eval=FALSE---------------------------------------------------------
#  ## Installing and loading devtools
#  install.packages( "devtools" )
#  library( "devtools" )
#  
#  ## Installing from github the fullsibQTL
#  install_github( "augusto-garcia/fullsibQTL" )

## ---- eval=FALSE---------------------------------------------------------
#  library(fullsibQTL)
#  fs_data <- read_onemap( dir = "C:/workingdirectory", inputfile = "filename.raw" )
#  
#  ## To read the example data from this same package:
#  fs_data <- read_onemap( dir = system.file("extdata", package = "fullsibQTL"),
#                          inputfile = "example_QTLfullsib.raw" )
#  
#  ## Within the Linkage Groups done:
#  fsib <- create_fullsib( fs_data,
#                          map.list = list( LG1_final, LG2_final, LG3_final, LG4_final ),
#                          step = 1, map.function = "kosambi" )

## ---- eval=FALSE---------------------------------------------------------
#  fs_data

## ---- eval=FALSE---------------------------------------------------------
#  fsib

## ---- eval=FALSE---------------------------------------------------------
#  system.file( "extdata", "example_QTLfullsib.raw", package = "fullsibQTL" )

## ---- eval=TRUE, echo=FALSE----------------------------------------------
rm(list=ls())
library(fullsibQTL)

## ---- eval=TRUE,cache=FALSE----------------------------------------------
## Loading the example data
data( "example_QTLfullsib" ) 

## Checking if the objects were loaded
ls() 

## Printing some of the loaded objects
example_QTLfullsib
LG1_final

## ---- eval=TRUE, cache=FALSE---------------------------------------------
fsib <- create_fullsib( example_QTLfullsib,
                        map.list = list( LG1_final, LG2_final, LG3_final, LG4_final ),
                       step = 1, map.function = "kosambi" )
## Printing it
fsib

## ---- eval=TRUE, cache=FALSE---------------------------------------------
im <- im_scan( fullsib = fsib, lg = "all", pheno.col = 1, LOD = TRUE )

## ---- eval=TRUE, cache=FALSE, R.options=list(max.print=50)---------------
im

## ---- eval=TRUE, cache=FALSE, R.options=list(max.print=50)---------------
print( im, lg = 1 )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
summary( im )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
plot( im )

## ---- eval=TRUE, cache=TRUE----------------------------------------------
plot_fullsibQTL( fullsib = fsib, fullsib.scan = im )

## ---- eval=FALSE, cache=FALSE--------------------------------------------
#  plot_fullsibQTL( fullsib = fsib, fullsib.scan = im, interact = TRUE )

## ---- eval=FALSE, cache=FALSE--------------------------------------------
#  vignette("Graphical_options_with_plot_fullsibQTL" ,package="fullsibQTL")

## ---- eval=FALSE, cache=FALSE--------------------------------------------
#  vignette("Graphical_options_with_plot" ,package="fullsibQTL")

## ---- eval=TRUE, cache=FALSE---------------------------------------------
covar <- matrix( rep( c( 1 ,-1 ), each = 150 ), ncol = 1 )
im_covar <- im_scan( fsib, pheno.col = 1, addcovar = covar )

## ---- eval=FALSE, cache=FALSE, collapse=TRUE-----------------------------
#  ## Perfoming the permutation
#  set.seed( 1234 )
#  im_perm <- im_scan( fsib, pheno.col = 1, n.perm = 1000,
#                      write.perm = "im_permutations.txt" )

## ---- eval=FALSE, cache=FALSE, collapse=TRUE-----------------------------
#  save( im_perm, file = "im_permutations.Rdata" )
#  load( "im_permutations.Rdata" )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
## Summary of the permutation test
summary( im_perm )
summary( im_perm, alpha = 0.05 )

## Getting the peaks for alpha 0.05.
im.peak1 <- summary( im_perm, alpha = 0.05 )[ 1, 1 ]
im.peak1
im.peak2 <- summary( im_perm,alpha = 0.05 )[ 1, 2 ]
im.peak2

## Plotting the second peak in the graphic
plot( im )
abline( h = im.peak2, col = "red" )

## Plotting both thresholds
plot( im, lty = 1, lwd = 2, incl.mkr = NULL, cex.incl = 0.7, 
      cex.axis = 0.8, col = "blue", ylab = "LOD Score", xlab = "Linkage Group", 
      main = "IM - Trait 1" )
abline( h = im.peak1, lty = 2, lwd = 1, col = "black" )
abline( h = im.peak2, lty = 2, lwd = 2, col = "black" )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
plot( im_perm, peak = 1 )
abline( v = im.peak1 )
abline( v = im.peak1, col = "red" )

plot( im_perm, peak = 2 )
abline( v = 3.419819, col = "blue" )
abline( v = 2.497431, col = "red" )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
summary( im )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
qtl1 <- im_char( fsib, pheno.col = 1, lg = 1, pos = "loc4" )
qtl2 <- im_char( fsib, pheno.col = 1, lg = 2, pos = "M27" )
qtl3 <- im_char( fsib, pheno.col = 1, lg = 3, pos = "loc64" )
qtl4 <- im_char( fsib, pheno.col = 1, lg = 4, pos = "loc56" )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
qtl1

## ---- eval=TRUE, cache=FALSE---------------------------------------------
cbind( qtl1, qtl2, qtl3, qtl4 )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
draw_phase( fsib, qtl1, 0.05 )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
get_segr( qtl1 )
get_segr( qtl2 )
get_segr( qtl3 )
get_segr( qtl4 )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
cofs_fs <- cof_selection( fsib, pheno.col = 1, k = log( 300 ), n.cofactor = 10 )
plot( cofs_fs )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
covar <- matrix( rep( c( 1, -1 ), each = 150 ), ncol = 1 )
cofs_fs2 <- cof_selection( fsib, pheno.col = 1, addcovar = covar, 
                           k = 2, thres.effect = 0.05 )

## ------------------------------------------------------------------------
cofs <- matrix(c("1","M2",
                 "2","M27",
                 "3","M33",
                 "4","M52"), 4, 2, byrow = TRUE )
cofs_fsdef <- cof_definition( fsib, pheno.col = 1, 
                              thres.effect = 0.05, cof.pos = cofs)

## ---- eval=TRUE, cache=FALSE---------------------------------------------
cim <- cim_scan( fullsib = cofs_fs, lg = "all", ws = 10, pheno.col = 1, LOD = TRUE )

## ---- eval=FALSE, cache=FALSE--------------------------------------------
#  set.seed( 4321 )
#  cim_perm <- cim_scan( fullsib = cofs_fs, lg = "all", pheno.col = 1,
#                        LOD = TRUE, n.perm = 1000,
#                        write.perm = "cim_permutations.txt" )

## ---- eval=FALSE, cache=FALSE, collapse=TRUE-----------------------------
#  save( cim_perm, file = "cim_permutations.Rdata" )
#  load( "cim_permutations.Rdata" )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
## Summary of the permutation test
summary( cim_perm )
## just with alpha=0.05 
summary( cim_perm, alpha = 0.05 )

## Getting the peaks for alpha 0.05. 
cim.peak1 <- summary(cim_perm, alpha=0.05)[1,1]
cim.peak1
cim.peak2 <- summary(cim_perm, alpha=0.05)[1,2]
cim.peak2

## Ploting the thrshold with the second peak
plot( cim )
abline( h = cim.peak2, col = "red" )

## Plotting both thresholds
plot( cim, lty = 1, lwd = 2, incl.mkr = NULL, cex.incl = 0.7, 
      cex.axis = 0.8, col = "red", ylab = "LOD Score", xlab = "Linkage Group", 
      main = "CIM - Trait 1" )
abline( h = cim.peak1, lty = 2, lwd = 1, col = "black" )
abline( h = cim.peak2, lty = 2, lwd = 2, col = "black" )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
plot( cim_perm, peak = 1 )
abline( v = cim.peak1, col = "red" )

plot( cim_perm, peak = 2 )
abline( v = cim.peak2, col = "red" )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
summary( cim )

## ---- eval=TRUE, cache=FALSE---------------------------------------------
## Plotting the CIM results
plot( cim, lty = 1, lwd = 2, incl.mkr = NULL, cex.incl = 0.7, 
      cex.axis = 0.8, col = "red", ylab = "LOD Score", xlab = "Linkage Group", 
      main = "CIM - Trait 1" )
abline( h = cim.peak1, lty = 2, lwd = 1, col = "black" )
abline( h = cim.peak2, lty = 2, lwd = 2, col = "black" )

## ---- eval=FALSE, cache=FALSE--------------------------------------------
#  plot_fullsibQTL( fullsib = fsib, fullsib.scan = cim, interact = TRUE )
#  # plot not shown

## ---- eval=TRUE, cache=FALSE---------------------------------------------
## Visual Analysis of every location of LG1
## QTLs peak at loc13 and M12
print( cim,lg = 1 )

## ---- eval=FALSE, cache=FALSE--------------------------------------------
#  ## Visual Analysis LG2
#  ## QTLs peaks at loc45, loc77, and loc 108
#  print( cim, lg = 2 )
#  
#  ## Visual Analysis LG3
#  ## QTLs peaks at M33 and loc65
#  print( cim,lg = 3 )
#  
#  ## Visual Analysis LG4
#  ## QTL peak at M52
#  print( cim,lg = 4 )

## ------------------------------------------------------------------------
print( cim, lg = 1, pos = c( "loc13", "M12", "loc45" ) )

## ---- eval=FALSE, cache=FALSE--------------------------------------------
#  print( cim, lg=2, pos = c( "loc45", "loc77", "loc108" ) )
#  print( cim, lg=3, pos = c( "M33", "loc65" ) )
#  print( cim, lg=4, pos = "M52" )

## ------------------------------------------------------------------------
cim_QTL1_lg1 <- cim_char( fullsib = cofs_fs, lg = 1, pheno.col = 1, pos = "loc13" )
cim_QTL2_lg1 <- cim_char( fullsib = cofs_fs, lg = 1, pheno.col = 1, pos = "M12" )
cim_QTL1_lg2 <- cim_char( fullsib = cofs_fs, lg = 2, pheno.col = 1, pos = "loc45" )
cim_QTL2_lg2 <- cim_char( fullsib = cofs_fs, lg = 2, pheno.col = 1, pos = "loc77" )
cim_QTL3_lg2 <- cim_char( fullsib = cofs_fs, lg = 2, pheno.col = 1, pos = "loc108" )
cim_QTL1_lg3 <- cim_char( fullsib = cofs_fs, lg = 3, pheno.col = 1, pos = "M33" )
cim_QTL2_lg3 <- cim_char( fullsib = cofs_fs, lg = 3, pheno.col = 1, pos = "loc65")
cim_QTL1_lg4 <- cim_char( fullsib = cofs_fs, lg = 4, pheno.col = 1 , pos = "M52" )

## Printing the object
cim_QTL1_lg1

## ------------------------------------------------------------------------
## Printing in a single table
cbind( cim_QTL1_lg1, cim_QTL2_lg1, cim_QTL1_lg2, cim_QTL2_lg2, 
       cim_QTL3_lg2, cim_QTL1_lg3, cim_QTL2_lg3, cim_QTL1_lg4)

## ------------------------------------------------------------------------
get_segr( cim_QTL1_lg1, probs1 = 0.05, probs2 = 0.05 )
get_segr( cim_QTL2_lg1, probs1 = 0.05, probs2 = 0.05 )
get_segr( cim_QTL1_lg2, probs1 = 0.05, probs2 = 0.05 )
get_segr( cim_QTL2_lg2, probs1 = 0.05, probs2 = 0.05 )
get_segr( cim_QTL3_lg2, probs1 = 0.05, probs2 = 0.05 )
get_segr( cim_QTL1_lg3, probs1 = 0.05, probs2 = 0.05 )
get_segr( cim_QTL2_lg3, probs1 = 0.05, probs2 = 0.05 )
get_segr( cim_QTL1_lg4, probs1 = 0.05, probs2 = 0.05 )

## ------------------------------------------------------------------------
draw_phase( cofs_fs, cim_QTL1_lg1 )
draw_phase( cofs_fs, cim_QTL2_lg1 )
draw_phase( cofs_fs, cim_QTL1_lg2 )
draw_phase( cofs_fs, cim_QTL2_lg2 )
draw_phase( cofs_fs, cim_QTL3_lg2 )
draw_phase( cofs_fs, cim_QTL1_lg3 )
draw_phase( cofs_fs, cim_QTL2_lg3 )
draw_phase( cofs_fs, cim_QTL1_lg4 )

## ------------------------------------------------------------------------
## Estimating QTL effects without covariate
qtls.cim <- r2_ls( fsib, pheno.col = 1, lg = c( 1, 1, 2, 2, 2, 3, 3, 4 ),
                   pos = c( "loc13", "M12", "loc45", "loc77", "loc108", "M33", "loc65", "M52" ) )

## ------------------------------------------------------------------------
## Estimating QTL effects with covariates
covar <- matrix( rep( c( 1, -1 ), each = 150 ), ncol = 1 )
r2_ls( fsib, pheno.col = 2, addcovar = covar, lg = c( 1, 1, 2, 2, 2, 3, 3, 4 ),
       pos = c( "loc13", "M12", "loc45", "loc77", "loc108", "M33", "loc65", "M52" ) )

## ------------------------------------------------------------------------
plot_fullsibQTL( fullsib = fsib,
                 fullsib.scan = cim,
                 r2ls.out = qtls.cim,
                 qtlmapping = "CIM",
                 thr = c( 4.023724, 2.900906 ) )

## ------------------------------------------------------------------------
sessionInfo()

