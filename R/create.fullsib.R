
#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: create.fullsib.R                                              #
#                                                                     # 
# Contains:                                                           #
# create.fullsib                                                      #
# print.fullsib                                                       #
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
# Function: create.fullsib                                            #
#                                                                     #
# Receives all the necessary information to perform QTL mapping       #
# One inserts the map, phenotype(s) and unlinked markers (if there is)#
# It is also required to define the step to obtain the multipoint prob#
# as well error.function and map function (can choice between 3)      #
# It returns a object from class 'fullsib'                            #
#######################################################################

create.fullsib <- function(input.obj, map.list, step=0, error.prob=1e-4,
                           map.function=c("kosambi", "haldane","morgan"), condIndex=3.5) {

  ## checking input.obj
  if (!any(class(input.obj) == "outcross"))
    stop(deparse(substitute(input.obj)), " is not an object of class 'outcross'")
  if(is.null(input.obj$pheno))
    stop("There is no phenotypes avaliable in ", deparse(substitute(input.obj)))

  ## checking map.list
  if (!any(class(map.list) == "list" | class(map.list) == "sequence")) 
    stop(deparse(substitute(map.list)),
         " is not an object of class 'list' or 'sequence'")
  if (class(map.list) == "sequence") 
    map.list <- list(map.list)

  ## checking object names
  ##pegue os nome do input.obj e coloque em x$data.name
  ##verifique se tem rf.2pts
  
  ##checking map.function
  map.function <- match.arg(map.function)

  ##obtaining unliked markers if any
  mapped.mkrs <- unlist(lapply(map.list, function(x) x$seq.num))
  not.mapped <- which(is.na(match(1:ncol(input.obj$geno), mapped.mkrs)))

  fullsib <- list(map = map.list,
                  unlinked.mkr = not.mapped, 
                  pheno = input.obj$pheno)

  ##obtaining multipoint probabilities used for im and/or cim analysis
  fullsib <- calc.probs(fullsib, step=step, error.prob=error.prob, map.function=map.function, condIndex)
  class(fullsib) <- "fullsib"
  fullsib

}

#######################################################################
#                                                                     #
# Function: print.fullsib                                             #
#                                                                     #
# Resume the main information on the object from fullsib class        #
#######################################################################

print.fullsib <- function (x,...) {
  LGtotal <- length(which(sapply(x$map, class) == "sequence"))
  #LG.length <- sapply(x$map, function(z) max(c(0, cumsum(get(.map.fun)(z$seq.rf)))))
  LG.length <- sapply(x$map, function(z) max(c(0, cumsum(kosambi(z$seq.rf)))))
  LG.mkrs <- sapply(x$map, function(z) length(z$seq.num))
  LG.dens <- LG.length/LG.mkrs
  cat("  This is an object of class 'fullsib'\n\n")
  ##cat("  The linkage map has", LGtotal, "groups, with", formatC(sum(LG.length, na.rm = TRUE)), "cM and", sum(LG.mkrs, na.rm = TRUE), "markers\n")
  cat("  The linkage map has", LGtotal, "groups, with", round(sum(LG.length, na.rm = TRUE),2), "cM and", sum(LG.mkrs, na.rm = TRUE), "markers\n")
  cat("  No. individuals genotyped: ",  get(x$map[[1]]$data.name)$n.ind, "\n")

  for (i in 1:length(x$map)){
    type.mkrs <- names(table(get(x$map[[i]]$data.name)$segr.type.num[x$map[[i]]$seq.num]))
    type.mkrs <- sapply(as.character(type.mkrs),
                        function(x){
                          switch(EXPR=x, "1" = "A", "2" = "B1", "3" = "B2",
                                 "4" = "B3", "5" = "C", "6" = "D1", "7" = "D2")})
    names(type.mkrs) <- NULL
    cat("  Group", formatC(i), ":", formatC(round(LG.length[i],2)), "cM, ", LG.mkrs[i], "markers",
        paste("(", paste(type.mkrs, collapse=", "), ")", sep=""), "\n")
  }
  if (!is.null(x$unlinked.mkr)) 
    cat("  And", length(x$unlinked.mkr), "unlinked markers\n")
  if (ncol(x$pheno) == 1) 
    cat("\n\n  One phenotype is avaliable for QTL mapping\n")
  else cat("\n\n ", ncol(x$pheno), "phenotypes are avaliable for QTL mapping\n")

  if(x$probs[[1]]$step == 0)
    cat("  Multipoint probability for QTL genotype was obtained at markers postions\n")
  else
    cat("  Multipoint probability for QTL genotype was obtained for each",
         x$probs[[1]]$step, "cM\n\n")

  if (any(class(x) == "fullsib.cofactors")){
    cat("\n  Cofactor selection was done for pheno.col =", x$cofactors$trait.cof, "\n")
    cat("  Markers selected for CIM analysis:", length(x$cofactors$names.cof$lg), "\n")
    
    entrance <- rownames(x$cofactors$names.cof)
    entrance <- cbind(entrance, x$cofactors$names.cof)
    
    lg.width <- nchar(as.character(as.character(entrance[,2])))
    mkr.width <- nchar(as.character(as.character(entrance[,3])))
    
    cat("  ", formatC("LG", width=max(lg.width)), rep("", 4),
        formatC("Marker", width=max(mkr.width)), "\n")  
    
    for(i in 1:nrow(entrance)){
      cat("  ", formatC(entrance[i,2], width=max(lg.width)),
          rep("", 4), as.character(entrance[i,3]), "\n")
    }
  }
}#end print.fullsib

