## ----knitr_init, echo=FALSE, cache=FALSE---------------------------------
library(knitr)

## From ggplot2 rmd vignette options
knitr::opts_chunk$set(collapse = FALSE,
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
plot_fullsibQTL( fullsib = fullsib,
                 fullsib.scan = cim )

## ------------------------------------------------------------------------
plot_fullsibQTL( fullsib = fullsib,
                 fullsib.scan = cim,
                 qtlmapping = "CIM" )

## ------------------------------------------------------------------------
plot_fullsibQTL( fullsib=fullsib,
                 fullsib.scan=list(cim,im),
                 qtlmapping = c( "CIM", "IM" ) )

## ------------------------------------------------------------------------
plot_fullsibQTL( fullsib = fullsib,
                 fullsib.scan = cim,
                 r2ls.out = qtls.cim,
                 qtlmapping = "CIM" )

## ------------------------------------------------------------------------
plot_fullsibQTL( fullsib = fullsib,
                 fullsib.scan = list( cim, im ),
                 r2ls.out = list( qtls.cim, qtls.im ),
                 qtlmapping = c( "CIM", "IM" ) )

## ------------------------------------------------------------------------
summary( cim_perm )

plot_fullsibQTL( fullsib = fullsib,
                 fullsib.scan = list ( cim, im ),
                 qtlmapping = c( "CIM", "IM" ),
                 thr = c(4.023724, 2.900906 ) )

## ------------------------------------------------------------------------
plot_fullsibQTL( fullsib = fullsib,
                 fullsib.scan = list( cim, im ),
                 r2ls.out = list( qtls.cim, qtls.im ),
                 qtlmapping = c( "CIM", "IM" ),
                 thr = c( 4.023724, 2.900906 ),
                 grayscale = TRUE )

## ------------------------------------------------------------------------
plot_fullsibQTL( fullsib = fullsib,
                 fullsib.scan = list( cim, im ),
                 r2ls.out = list( qtls.cim, qtls.im ),
                 qtlmapping = c( "CIM", "IM" ),
                 thr = c( 4.023724, 2.900906 ),
                 lgs = 2 )

## ------------------------------------------------------------------------
plot_fullsibQTL( fullsib = fullsib,
                 fullsib.scan = list( cim, im ),
                 r2ls.out = list( qtls.cim, qtls.im ),
                 qtlmapping = c( "CIM", "IM" ),
                 thr = c( 4.023724, 2.900906 ),
                 lgs = c( 2, 4 ) )

## ---- eval=FALSE---------------------------------------------------------
#  ## Setting the file name
#  plot_fullsibQTL( fullsib = fullsib,
#                   fullsib.scan = list( cim, im ),
#                   r2ls.out = list( qtls.cim, qtls.im ),
#                   qtlmapping = c ( "CIM", "IM" ),
#                   thr = c( 4.023724, 2.900906 ),
#                   interact = TRUE,
#                   file = "cim_im.html" )
#  
#  ## Setting the folder name
#  plot_fullsibQTL( fullsib = fullsib,
#                   fullsib.scan = list( cim, im ),
#                   r2ls.out = list( qtls.cim, qtls.im ),
#                   qtlmapping = c( "CIM", "IM" ),
#                   thr = c( 4.023724, 2.900906 ),
#                   interact = TRUE,
#                   file = "cim_im.html",
#                   folder = "/home/your_user_name/graphics" )
#  
#  ## Setting to not open the browser
#  plot_fullsibQTL( fullsib = fullsib,
#                   fullsib.scan = list( cim, im ),
#                   r2ls.out = list( qtls.cim, qtls.im ),
#                   qtlmapping = c( "CIM", "IM" ),
#                   thr = c( 4.023724, 2.900906 ),
#                   interact = TRUE,
#                   file = "cim_im.html",
#                   folder = "/home/your_user_name/Pictures",
#                   browser = FALSE )

## ---- eval=FALSE---------------------------------------------------------
#  plot_fullsibQTL( fullsib = fullsib,
#                   fullsib.scan = list( cim, im ),
#                   r2ls.out = list( qtls.cim, qtls.im ),
#                   qtlmapping = c( "CIM", "IM" ),
#                   thr = c( 4.023724, 2.900906 ),
#                   lgs = c( 2, 4 ),
#                   interact = TRUE,
#                   lg.colors = c( "darkgreen", "green" ),
#                   height = 600,
#                   width = 800,
#                   file = "cim_im_greens.html",
#                   folder="/home/your_user_name/graphics" )

## ------------------------------------------------------------------------
sessionInfo()

