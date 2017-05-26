 
#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: im_char.R                                                     #
#                                                                     # 
# Contains:                                                           #
# im_char                                                             #
#                                                                     #
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
# Function: im_char                                                   #
#                                                                     #
# This function receives the position that qtl caracterization is done#
# i.e., one estimate the genetic effects of QTL and test if they are  #
# significative, and segregation. Linkage phase is obtained by the    #
# intepretation of QTL effects signals and better view with function  #
# draw_phase that used the result of this function...                 #
# all the MLE and tests are done in C code                            #
#######################################################################

im_char <- function(fullsib, pheno.col=1, addcovar=NULL, lg, pos,
                    maxit=1000, tol=1e-10, verbose=FALSE){
  
  ##checking input arguments
  if (!any(class(fullsib) == "fullsib"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib'")
  
  if (is.na(match(pheno.col, 1:ncol(fullsib$pheno))) )
    stop("pheno.col must be a number between 1 and",ncol(fullsib$pheno))
  
  if(missing(lg))
    stop("lg argument is missing")
  if(!(lg %in% seq(1,length(fullsib$map),1)))
    stop("lg must be an integer between 1 and ", length(fullsib$map))
  
  if(missing(pos))
    stop("pos argument is missing")
  ##it can be a number from 1 to last position of the group, good to loops...
  if(is.numeric(pos) && !(pos %in% seq(1,length(fullsib$probs[[lg]]$newmap$dist),1)))
    stop("pos argument may be an integer between 1 and ",
         length(fullsib$probs[[lg]]$newmap$dist))
  
  if(is.character(pos)){
    pos.arg <- match.arg(pos, names(fullsib$probs[[lg]]$newmap$dist))
    if(!(pos.arg %in% names(fullsib$probs[[lg]]$newmap$dist)))     
      stop("pos argument may be an string that contain the name of the loci: please check it")
    pos <- which(names(fullsib$probs[[lg]]$newmap$dist) == pos.arg)
  }


  if(is.null(addcovar)==TRUE){
    n.ac <- 0
    ac <- NULL
  }
  else{
    ac <- as.matrix(addcovar)
    n.ac <- ncol(ac)
  }
  
  if (maxit < 1) 
    stop("maxit must be >= 1")
  
  if (tol <= 0) 
    stop("tol must be > 0")
  
  
  ##keeping non NA elements
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]
  
  if(n.ac > 0) ac <- ac[pheno.index,]
  
  if(n.ac > 0)
    lm1 <- lm(pheno ~ ac)
  else
    lm1 <- lm(pheno ~ 1)

  ac <- model.matrix(lm1)
  n.ac <- ncol(ac)
  Zls <- t(solve(crossprod(ac), t(ac)))
  resid0 <- lm1$resid
  gamma <- coef(lm1)

  ##getting the flanking markers or closest markers
  name.pos <- names(fullsib$probs[[lg]]$newmap$dist[pos])
  if(substring(name.pos, 1, 3) == "loc"){
    mkrs <- which(substring(names(fullsib$probs[[lg]]$newmap$dist),1,3) != "loc")
    flank <- which(cumsum(pos < mkrs) == 1)
    mk.rgt <- names(fullsib$probs[[lg]]$newmap$dist)[ mkrs[flank] ]
    mk.lft <- names(fullsib$probs[[lg]]$newmap$dist)[ mkrs[flank - 1] ]
    closest.mkr <-  paste(mk.lft,mk.rgt, sep="-")
  }
  else
    closest.mkr <- name.pos

  z <- .C("R_char_qtl",
          as.integer(length(pheno.index)),
          as.integer(fullsib$probs[[lg]]$colin$geno.class[pos]),
          as.integer(fullsib$probs[[lg]]$colin$type[pos]),      
          as.double(fullsib$probs[[lg]]$cond.prob[pheno.index,pos,]),
          as.double(ac),
          as.integer(n.ac),
          as.double(gamma),
          as.double(Zls),
          as.double(pheno),        
          result=as.double(rep(0, 12)), #this is LOD
          as.integer(maxit),
          as.double(tol),
          as.integer(verbose), PACKAGE="fullsibQTL")$result
  ##as.integer(verbose))$result 
  z <- as.matrix(z)

  pval <- -log10(pchisq(z[1]*4.61,
                        fullsib$probs[[lg]]$colin$geno.class[pos] - 1,
                        lower.tail = FALSE)) 

  ##convert LOD to pvalue: H04, H05 and H06
  z[9:11] <- pchisq(z[9:11]*4.61, 1, lower.tail = FALSE)
  z <- rbind(lg, fullsib$probs[[lg]]$newmap$dist[pos], pval, z)
 
  rownames(z) <- c("LG", "pos", "-log10(pval)","LOD_Ha", "mu",
                   "alpha_p", "LOD_H1","alpha_q","LOD_H2", "delta_pq",
                   "LOD_H3", "H4_pvalue", "H5_pvalue", "H6_pvalue", "model")
  
  colnames(z) <- closest.mkr
  structure(z, class = c("fullsib_char", "matrix"))
}#end function
