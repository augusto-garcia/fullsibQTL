#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: draw_phase.R                                                  #
#                                                                     # 
# Contains: draw_phase                                                #
#                                                                     #
#                                                                     #
# An adaption by Rodrigo Gazaffi, for the print.sequence function     #
# present in onemap package                                           #
#                                                                     #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
# Updated by Rodrigo Amadeu                                           #
#                                                                     #
# First version: 09/30/2011 (american date format)                    #
# Last  version: 08/15/2011 (american date format)                    #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: draw_phase.R                                              #
#                                                                     # 
# Identifies the linkage phase between QTL and markers                #
# It represents the linkage group, with the mapped QTL.               #
# the linkage group information is in fullsib object and              #
# QTL is in fschar, probs is used to identify if the allele is        #
# significative or not                                                #
#######################################################################

## -------------------------
## draw_phase function

#' Represents the linkage phase between markers and QTL
#' 
#' Locate the QTL in a linkage group, representing the significative parental
#' effects from the result of \code{im_char} and \code{cim_char}
#' 
#' @param fullsib An object from class \emph{fullsib} or
#' \emph{fullsib_cofactors}.
#' @param fschar An object from class \emph{fullsib_char}.
#' @param probs probability used to identify significative parental effect.
#' 
#' @return It returns text output indicating in a linkage group linkage phase
#' between markers and QTL.
#' 
#' @author Based on \code{print.sequence} from \pkg{onemap} and it was written
#' by Marcelo Mollinari \email{mmollina@@usp.br}. It was modified allowing the
#' inclusion of QTL, by Rodrigo Gazaffi \email{rgazaffi@@gmail.com}
#' 
#' @seealso
#' \code{\link[fullsibQTL]{im_char}}
#' \code{\link[fullsibQTL]{cim_char}}
#' \code{\link[fullsibQTL]{get_segr}}
#' @keywords utilities
#' @examples
#' 
#' 
#'   data( "example_QTLfullsib" )
#' 
#'   fullsib <- create_fullsib( example_QTLfullsib,
#'                              list( LG1_final, LG2_final, LG3_final, LG4_final ),
#'                              step = 0, map.function = "kosambi",condIndex = 3.5 )
#' 
#' 
#'   ###############################################
#'   ## cofactor selection using BIC (n.ind = 300)
#'   cofs.fs <- cof_selection( fullsib, pheno.col = 1, k = log( 300 ),
#'                             selection = 1 ) 
#' 
#'   \dontrun{
#'   cim1 <- cim_scan( cofs.fs, pheno.col = 1, ws = 22, LOD = TRUE, icim = FALSE )
#'   summary( cim1 )
#'   plot( cim1 )
#'   }
#' 
#'   qtl.lg3 <- cim_char( cofs.fs, pheno.col = 1, ws = 22, lg = 3, pos = "M38" )
#'   draw_phase( cofs.fs, qtl.lg3, probs = 0.05 )
#' 
#'   qtl.lg4 <- cim_char( cofs.fs, pheno.col = 1, ws = 22, lg = 4, pos = "M52" )
#'   draw_phase( cofs.fs, qtl.lg4, probs = 0.05 )
#' 
 
draw_phase <- function(fullsib, fschar, probs=0.05){
  
  ##checking argument
  if (!any(class(fullsib) == "fullsib"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib'")

  ##checking argument
  if (!any(class(fschar) == "fullsib_char"))
    stop(sQuote(deparse(substitute(fschar))),
         " is not an object of class 'fullsib_char'")

  ##checking argument
  if(probs > 1 || probs < 0){
    break("probs argument need to be a probability between 0 and 1")
    probs <- 0
  }

  ##extracting argument from fullsib_char object
  LG <- fschar["LG",]
  qtl.cM <- fschar["pos", ]
  used.model <- fschar["model",]
  alpha.p <- fschar["alpha_p",]
  lod.ap <- fschar["LOD_H1",]
  alpha.q <- fschar["alpha_q",]
  lod.aq <- fschar["LOD_H2",]
  
  ## returning the linkage phase for QTL
  QTL.phase <- get_phase(used.model, alpha.p, lod.ap, alpha.q, lod.aq, probs)

  ## finding where qtl will be draw...
  mkrs <- which(substring(names(fullsib$probs[[LG]]$newmap$dist),1,3) != "loc")
  #distances <- c(0, cumsum(get(.map.fun)(fullsib$map[[LG]]$seq.rf)))
  distances <- c(0, cumsum(kosambi(fullsib$map[[LG]]$seq.rf)))
  mk.rgt <- which(cumsum(qtl.cM < distances)==1)

  ## print.seq
  marnames <- colnames(fullsib$map[[LG]]$data.name$geno)[fullsib$map[[LG]]$seq.num]
  
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
    parents[i,] <- return_geno(fullsib$map[[LG]]$data.name$segr.type[fullsib$map[[LG]]$seq.num[i]],link.phases[i])

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