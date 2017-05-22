#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: draw.phase.R                                                  #
#                                                                     # 
# Contains: draw.phase                                                #
#                                                                     #
#                                                                     #
# An adaption by Rodrigo Gazaffi, for the print.sequence function     #
# present in onemap package                                           #
#                                                                     #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
#                                                                     #
# First version: 09/30/2011 (american date format)                    #
# Last  version: 09/30/2011 (american date format)                    #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: draw.phase.R                                              #
#                                                                     # 
# Identifies the linkage phase between QTL and markers                #
# It represents the linkage group, with the mapped QTL.               #
# the linkage group information is in fullsib object and              #
# QTL is in fschar, probs is used to identify if the allele is        #
# significative or not                                                #
#######################################################################


draw.phase <- function(fullsib, fschar, probs=0.05){
  
  ##checking argument
  if (!any(class(fullsib) == "fullsib"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib'")

  ##checking argument
  if (!any(class(fschar) == "fullsib.char"))
    stop(sQuote(deparse(substitute(fschar))),
         " is not an object of class 'fullsib.char'")

  ##checking argument
  if(probs > 1 || probs < 0){
    break("probs argument need to be a probability between 0 and 1")
    probs <- 0
  }

  ##extracting argument from fullsib.char object
  LG <- fschar["LG",]
  qtl.cM <- fschar["pos", ]
  used.model <- fschar["model",]
  alpha.p <- fschar["alpha_p",]
  lod.ap <- fschar["LOD_H1",]
  alpha.q <- fschar["alpha_q",]
  lod.aq <- fschar["LOD_H2",]
  
  ## returning the linkage phase for QTL
  QTL.phase <- get.phase(used.model, alpha.p, lod.ap, alpha.q, lod.aq, probs)

  ## finding where qtl will be draw...
  mkrs <- which(substring(names(fullsib$probs[[LG]]$newmap$dist),1,3) != "loc")
  #distances <- c(0, cumsum(get(.map.fun)(fullsib$map[[LG]]$seq.rf)))
  distances <- c(0, cumsum(kosambi(fullsib$map[[LG]]$seq.rf)))
  mk.rgt <- which(cumsum(qtl.cM < distances)==1)

  ## print.seq
  marnames <- colnames(get(fullsib$map[[LG]]$data.name)$geno)[fullsib$map[[LG]]$seq.num]
  
  link.phases <- matrix(NA,length(fullsib$map[[LG]]$seq.num),2)
  link.phases[1,] <- rep(1,2)
  for (i in 1:length(fullsib$map[[LG]]$seq.phases)) {
    switch(EXPR=fullsib$map[[LG]]$seq.phases[i],
           link.phases[i+1,] <- link.phases[i,]*c(1,1),
           link.phases[i+1,] <- link.phases[i,]*c(1,-1),
           link.phases[i+1,] <- link.phases[i,]*c(-1,1),
           link.phases[i+1,] <- link.phases[i,]*c(-1,-1),
           )
  }
  link.phases <- apply(link.phases,1,function(x) paste(as.character(x),collapse="."))
  parents <- matrix("",length(fullsib$map[[LG]]$seq.num),4)
  
  for (i in 1:length(fullsib$map[[LG]]$seq.num))
    parents[i,] <- return.geno(get(fullsib$map[[LG]]$data.name)$segr.type[fullsib$map[[LG]]$seq.num[i]],link.phases[i])

  ##adapting for QTL printing
  longest.name <- max(nchar(c("QTL", marnames)))                          
  marnames.qtl <- formatC(c("QTL", marnames),flag="-")
  qtl.label <- marnames.qtl[1]
  marnames <- marnames.qtl[-1]
  distances <-formatC(distances,format="f",digits=2,width=7)
  distances.qtl <- formatC(qtl.cM, format="f",digits=2,width=7)
  
  cat("\nPrinting QTL and its linkage phase between markers across the LG:\n\n")
  cat("Markers",rep("",max(longest.name-7,0)+9),"Position",rep("",10),"Parent 1","     ","Parent 2\n\n")

  
  if(length(mk.rgt) == 1){
    ##before the QTL
    for (i in 1: (mk.rgt-1)){
      cat(marnames[i],rep("",max(7-longest.name,0)+10),distances[i],rep("",10),parents[i,1],"|  |",parents[i,2],"     ",parents[i,3],"|  |",parents[i,4],"\n")
    }
    
    ##QTL information - keypoint...
    cat(qtl.label,rep("",max(7-longest.name,0)+10), distances.qtl, rep("",9), QTL.phase[1],"|  |", QTL.phase[2],"   ", QTL.phase[3],"|  |", QTL.phase[4],"\n")

    ##after the QTL
    for (i in mk.rgt: length(fullsib$map[[LG]]$seq.num))
      cat(marnames[i],rep("",max(7-longest.name,0)+10),distances[i],rep("",10),parents[i,1],"|  |",parents[i,2],"     ",parents[i,3],"|  |",parents[i,4],"\n")
  }
  else {#QTL mapped in the last position
    ##before QTL
    for (i in 1: length(fullsib$map[[LG]]$seq.num))
      cat(marnames[i],rep("",max(7-longest.name,0)+10),distances[i],rep("",10),parents[i,1],"|  |",parents[i,2],"     ",parents[i,3],"|  |",parents[i,4],"\n")

    ##QTL information
    cat(qtl.label,rep("",max(7-longest.name,0)+10), distances.qtl, rep("",9), QTL.phase[1],"|  |", QTL.phase[2],"   ", QTL.phase[3],"|  |", QTL.phase[4],"\n")
  }
  
  cat("\nP1 and Q1 has positive effect (increase phenotypic value)\n")
  cat("P2 and Q2 has negative effect (reduce phenotypic value) \n")
  if(length(which(substring(QTL.phase,2,2) == "0")) > 1)
    cat("P0 and Q0 has neutral effect (non signif.)\n")
}


