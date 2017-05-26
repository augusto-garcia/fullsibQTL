#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: cim_char.R                                                    #
#                                                                     # 
# Contains: cim_char                                                  #
#                                                                     #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
#                                                                     #
# First version: 2011/30/09                                           #
# Last  version: 2011/30/09                                           #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: cim_char                                                  #
#                                                                     #
# This function receives the position that qtl caracterization is done#
# i.e., one estimate the genetic effects of QTL and test if they are  #
# significative, and segregation. Linkage phase is obtained by the    #
# intepretation of QTL effects signals and better view with function  #
# draw_phase that used the result of this function.                   #
# all the MLE and tests are done in C code                            #
#######################################################################

cim_char <- function(fullsib, pheno.col=1, ws = 10, lg, pos,
                     maxit=1000,tol=1e-10, icim=FALSE, verbose=FALSE)
{
  ##checking input arguments
  if (!any(class(fullsib) == "fullsib_cofactors"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib_cofactors': execute 'cof.selection' or 'cof_definition' first")
  
  if(ws < 0)
    stop("window size must be non negative number")
  ws <- ws / 2 #similar to qtl
  
  if (is.na(match(pheno.col, 1:ncol(fullsib$pheno))) )
    stop("pheno.col must be a number between 1 and", ncol(fullsib$pheno))
  
  if(missing(lg))
    stop("lg argument is missing")
  if(!(lg %in% seq(1,length(fullsib$map),1)))
    stop("lg must be an integer between 1 and ", length(fullsib$map))
  
  if(missing(pos))
    stop("pos argument is missing")
  if(is.character(pos)){
    pos.arg <- match.arg(pos, names(fullsib$probs[[lg]]$newmap$dist))
    if(!(pos.arg %in% names(fullsib$probs[[lg]]$newmap$dist)))  
      stop("pos argument may be an string that contain the name of the loci: please check this name")
    pos <- which( names(fullsib$probs[[lg]]$newmap$dist) == pos.arg)
  } 
  
  if(pheno.col != fullsib$cofactors$trait.cof)
    stop("pheno.col argument is different than pheno.col used to obtain cofactors")
  
  ##keeping non NA elements
  pheno <- fullsib$pheno[,pheno.col]
  pheno.index <- which(!is.na(pheno))
  pheno <- pheno[pheno.index]
  if(icim == TRUE)
    pheno.icim <- pheno
  
  ##obtaining the model that will be used for that position
  cim.models <- model_select_lg(fullsib, ws = ws, LG=lg)

  for(mods in 1:length(cim.models)){
    find.pos <- which(pos %in% cim.models[[mods]]$position)
    if(length(find.pos) > 0)
      flag.mkr <- mods
  }

  if(cim.models[[flag.mkr]]$mkr2rm[1] == "none")
    cim.cof.matrix <- fullsib$cofactors$matrix.cof[pheno.index,]
  else{
    mkr2rm <-do.call("c",lapply(cim.models[[flag.mkr]]$mkr2rm,
           function(x) which(colnames(fullsib$cofactors$matrix.cof) == x)))

    cim.cof.matrix <- fullsib$cofactors$matrix.cof[pheno.index,-mkr2rm]
    
  }

  
  unar <- matrix(1,nrow(cim.cof.matrix))
  colnames(unar) <- "intercept"
  cim.cof.matrix <- cbind(unar,cim.cof.matrix)
  
  if(icim == FALSE){
    ##traditional analysis
    Zls <- t(solve(crossprod(cim.cof.matrix), t(cim.cof.matrix)))
    gamma <- coef(lm(pheno ~ cim.cof.matrix -1))
  }
  else{
    ##icim proposed -  remove cofactors from phenotype
    pheno.rm <- lm(pheno.icim ~ cim.cof.matrix - 1)
    pheno <- pheno.rm$residuals
    gamma <- coef(pheno.rm)[1]#for intercept
    Zls <- t(solve(crossprod(unar), t(unar)))
    cim.cof.matrix <- unar
  }
  
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
  
  dim_probs <- dim(fullsib$probs[[lg]]$cond.prob)
  z <- .C("R_char_qtl",
          as.integer(length(pheno.index)),        
          as.integer(fullsib$probs[[lg]]$colin$geno.class[pos]),
          as.integer(fullsib$probs[[lg]]$colin$type[pos]),      
          as.double(fullsib$probs[[lg]]$cond.prob[pheno.index,pos,]),
          as.double(cim.cof.matrix),
          as.integer(ncol(cim.cof.matrix)),
          as.double(gamma),#
          as.double(Zls), #
          as.double(pheno),       
          result=as.double(rep(0, 12)),
          as.integer(maxit),
          as.double(tol),
          as.integer(verbose),PACKAGE="fullsibQTL")$result
          ##as.integer(verbose))$result 
  z <- as.matrix(z)


  pval <-  -log10(pchisq(z[1]*4.61,  fullsib$probs[[lg]]$colin$geno.class[pos] - 1, lower.tail = FALSE))
  #print(z[9:11])
  z[9:11] <- pchisq(z[9:11]*4.61, 1, lower.tail = FALSE)
  #print(z[9:11])
  
  z <- rbind(lg, fullsib$probs[[lg]]$newmap$dist[pos], pval, z)

  rownames(z) <- c("LG", "pos", "-log10(pval)", "LOD_Ha", "mu",
                   "alpha_p", "LOD_H1", "alpha_q", "LOD_H2", "delta_pq",
                   "LOD_H3", "H4_pvalue", "H5_pvalue", "H6_pvalue", "model")
  colnames(z) <- closest.mkr

  structure(z, class = c("fullsib_char", "matrix"))

}


