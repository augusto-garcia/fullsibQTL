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
# Updated by Rodrigo Amadeu and Cristiane Taniguti                    #
#                                                                     #
# First version: 09/30/2011 (american date format)                    #
# Last  version: 11/08/2020 (american date format)                    #
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
# QTL is in fullsib.char, probs is used to identify if the allele is        #
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
#' @param fullsib.char An object from class \emph{fullsib_char}.
#' @param fullsib.scan An object from class \emph{fullsib_scan}.
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
#' @export 
draw_phase <- function(fullsib, fullsib.char, probs=0.05, fullsib.scan = NULL){
  
  ##checking argument
  if (!any(class(fullsib) == "fullsib"))
    stop(sQuote(deparse(substitute(fullsib))),
         " is not an object of class 'fullsib'")
  
  ##checking argument
  if (!any(class(fullsib.char) == "fullsib_char"))
    stop(sQuote(deparse(substitute(fullsib.char))),
         " is not an object of class 'fullsib_char'")
  
  ##checking argument
  if (!(is(fullsib.scan, "fullsib_scan") | is.null(fullsib.scan)))
    stop(sQuote(deparse(substitute(fullsib.scan))),
         " is not an object of class 'fullsib_scan'")
  
  ##checking argument
  if(probs > 1 || probs < 0){
    stop("probs argument need to be a probability between 0 and 1")
    probs <- 0
  }
  
  ##extracting argument from fullsib_char object
  LG <- fullsib.char["LG",]
  qtl.cM <- fullsib.char["pos", ]
  used.model <- fullsib.char["model",]
  alpha.p <- fullsib.char["alpha_p",]
  lod.ap <- fullsib.char["LOD_H1",]
  alpha.q <- fullsib.char["alpha_q",]
  lod.aq <- fullsib.char["LOD_H2",]
  
  ## returning the linkage phase for QTL
  QTL.phase <- get_phase(used.model, alpha.p, lod.ap, alpha.q, lod.aq, probs)
  
  ## finding where qtl will be draw...
  mkrs <- which(substring(names(fullsib$probs[[LG]]$newmap$dist),1,3) != "loc")
  #distances <- c(0, cumsum(get(.map.fun)(fullsib$map[[LG]]$seq.rf)))
  distances <- c(0, cumsum(kosambi(fullsib$map[[LG]]$seq.rf)))
  
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
  marnames.qtl <- formatC(c("QTL", marnames),flag="-")
  qtl.label <- marnames.qtl[1]
  marnames <- marnames.qtl[-1]
  
  all.dist <- unique(c(qtl.cM, distances))
  mk.rm <- length(c(qtl.cM, distances)) != length(all.dist)
  all.dist.sort <- sort(all.dist)
  pos.qtl <- which(match(all.dist.sort, all.dist) == 1)
  
  obj.phases <- matrix(NA, nrow= length(all.dist), ncol=6)
  
  if(mk.rm)  obj.phases[,1] <- marnames else obj.phases[-pos.qtl,1] <- marnames 
  obj.phases[pos.qtl,1] <- qtl.label
  obj.phases[,2] <- all.dist.sort
  if(mk.rm) obj.phases[,3:6] <- parents else obj.phases[-pos.qtl,3:6]  <- parents
  obj.phases[pos.qtl,3:6] <- QTL.phase
  
  # Confidence interval
  if(!is.null(fullsib.scan)){
    LOD1 <- fullsib.char[4,1] - 1
    LOD2 <- fullsib.char[4,1] - 2
    pos.qtl <- which(fullsib.scan[,1] == fullsib.char[1,1] & fullsib.scan[,2] == fullsib.char[2,1])
    lg <- fullsib.scan[fullsib.scan[,1] == fullsib.char[1,1],3]
    
    LOD1.right <- LOD2.right <- vector()
    i <- pos.qtl 
    while(i <= length(lg)){
      if(lg[i] > LOD1) LOD1.right <- c(LOD1.right, i) else break
      i <- i + 1
    }
    
    i <- pos.qtl
    while(i <= length(lg)){
      if(lg[i] > LOD2) LOD2.right <- c(LOD2.right, i) else break
      i <- i + 1
    }
    
    LOD1.left <- LOD2.left <- vector()
    i <- pos.qtl -1
    while(i != 0){
      if(lg[i] > LOD1) LOD1.left <- c(LOD1.left, i) else break
      i <- i - 1
    }
    
    i <- pos.qtl -1 
    while(i != 0){
      if(lg[i] > LOD2) LOD2.left <- c(LOD2.left, i) else break
      i <- i - 1
    }
    
    LOD1 <- fullsib.scan[sort(c(LOD1.left, LOD1.right)),2]
    LOD2 <- fullsib.scan[sort(c(LOD2.left, LOD2.right)),2]
  } else {
    LOD1 <- NA
    LOD2 <- NA
  }
  
  colnames(obj.phases) <- c("markers", "positions", "P1_1", "P1_2", "P2_1", "P2_2")
  obj.phases <- as.data.frame(obj.phases)
  obj.phases[,2] <- as.numeric(as.character(obj.phases[,2]))
  obj.phases <- cbind(obj.phases, LG=fullsib.char[1,1])
  obj.phases <- list(obj.phases, LOD1 = c(min(LOD1), max(LOD1)), LOD2 = c(min(LOD2), max(LOD2)))
  class(obj.phases) <- "fullsib_phases"
  return(obj.phases)
}

#' Print markers and QTL phases
#' 
#' @param x object of class fullsib_phases
#' @param ... currently ignored
#' 
#' @method print fullsib_phases
#' @export
print.fullsib_phases <- function(x, ...){
  x <- x[[1]]
  longest.name <- max(nchar(as.character(x[,1])))                          
  x[,c(1,3:6)] <- apply(x[,c(1,3:6)], 2, function(x) formatC(x,flag="-"))
  x[,2] <- formatC(round(x[,2],2), format="f",digits=2,width=7)
  
  cat("\nPrinting QTL and its linkage phase between markers across the LG:\n\n")
  cat("Markers",rep("",max(longest.name-7,0)+9),"Position",rep("",10),"Parent 1","     ","Parent 2\n\n")
  
  for (i in 1:dim(x)[1]){
    cat(x[i,1],rep("",max(7-longest.name,0)+10),x[i,2],rep("",10),x[i,3],"|  |",x[i,4],"     ",x[i,5],"|  |",x[i,6],"\n")
  }
  
  cat("\nP1 and Q1 have positive effect (increase phenotypic value)\n")
  cat("P2 and Q2 have negative effect (reduce phenotypic value) \n")
  cat("P0 and Q0 have neutral effect (non signif.)\n")
}

#' Plot markers and QTL phases
#' 
#' @param x object or list of objects of class fullsib_phases
#' @param CI plots confidence interval using 1 LOD (LOD1) or 2 LODs (LOD2)
#' @param ... currently ignored
#' 
#' @rawNamespace import(ggplot2, except = last_plot)
#' @importFrom reshape2 melt
#' 
#' @export
plot_fullsib_phases <- function(x, CI = NULL, ...){
  if(!is(x, "fullsib_phases")){
    if(is(x, "list")){
      if(!all(sapply(x, function(y) is(y, "fullsib_phases")))) 
        stop("All objects in list must be of class fullsib_phases.\n")
      if(!is.null(CI)) 
        stop("This function don't support plot CI for more than one fullsib_phases object.\n")
      group <- unique(unlist(as.vector(sapply(x, function(y) y[[1]][,7]))))
      if(length(group) != 1) 
        stop("This function only plots multiple QTLs in the sample linkage group. 
              If you want to plot multiple groups, check ggarrange function from ggpubr package. \n")
      qtl.list <- lapply(x, function(y) y[[1]][which(y[[1]][,1] == "QTL"),])
      qtl.list <- do.call(rbind, qtl.list)
      qtl.list <- qtl.list[,-7]
      x <- x[[1]]
    }
  }
  
  if(!is.null(CI)){
    if(!(CI %in% c("LOD1", "LOD2"))){
      stop("CI argument can receive only 'NULL', 'LOD1' and 'LOD2'\n")
    }
    if(is.na(x[[2]][1])) stop("To plot the confidence interval you must
                              add the fullsib_scan object to run the draw_phases.\n")
  }
  
  df.tot <- x[[1]]
  if(!exists("qtl.list")) qtl.list <- df.tot[which(df.tot[,1] == "QTL"),-7]
  qtl.list[,3:6] <- apply(qtl.list[,3:6], 2, as.character)
  qtl.pos <- melt(qtl.list, measure.vars = colnames(qtl.list)[-c(1,2)], variable.name = "haplo")
  qtl.pos <- cbind(qtl.pos, parent = sapply(strsplit(as.character(qtl.pos[,3]), "_"), "[[",1))
  
  df <- df.tot[-which(df.tot[,1] == "QTL"),-7]
  df[,2] <- c(0,df[,2][-1]-diff(df[,2])/2)
  
  df[,3:6] <- apply(df[,3:6], 2, as.character)
  df_melt <- melt(df[,-1], measure.vars = colnames(df[,-c(1,2)]))
  df_melt$variable <- factor(df_melt$variable, levels = sort(levels(df_melt$variable)))
  df_melt$value <- as.numeric(as.factor(df_melt$value))
  
  # Creating alpha
  alpha <- vector()
  for(i in 1:length(df_melt$value)){
    if(df_melt$value[i] == 1){
      alpha <- c(alpha, c(1,0,0,0,0))
    } else if(df_melt$value[i] == 2){
      alpha <- c(alpha, c(0,1,0,0,0))
    } else if(df_melt$value[i] == 3){
      alpha <- c(alpha, c(0,0,1,0,0))
    }else if(df_melt$value[i] == 4){
      alpha <- c(alpha, c(0,0,0,1,0))
    }else if(df_melt$value[i] == 5){
      alpha <- c(alpha, c(0,0,0,0,1))
    }
  }
  
  df_melt <- data.frame(pos = rep(df_melt[,1], each=5), 
                        haplo = rep(df_melt[,2], each=5),
                        allele = c("a", "b", "c", "d", "o"), 
                        alpha,
                        parent = rep(sapply(strsplit(as.character(df_melt[,2]), "_"), "[[",1), each=5))
  
  shapes <- c("\u25B2","\u25BC","\u25C6", "\u25B3", "\u25BD", "\u25C7")
  names(shapes) <- c("P1", "P2", "P0", "Q1", "Q2", "Q0")
  
  p <- ggplot(df_melt, aes(x = pos))  +
    geom_line(aes(y = haplo, col= allele, alpha = alpha), size = 5) +
    scale_alpha_continuous(range = c(0,1)) +
    geom_point(data = qtl.pos, aes(x = positions, 
                                   y = haplo, shape = factor(value)), 
               size = 7, na.rm = T) + 
    scale_color_viridis_d() +
    scale_shape_manual(values= shapes, drop=F) +
    labs(col = "Alleles", x = "position (cM)", y = element_blank(), shape="QTL effects") +
         # caption = "\n\nP1 and Q1 have positive effect (increase phenotypic value)\n 
         # P2 and Q2 have negative effect (reduce phenotypic value) \n
         # P0 and Q0 have neutral effect (non signif.)\n") + 
    theme(panel.background = element_blank(), 
          legend.position="bottom", 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    guides(alpha=FALSE) 
  
  if(!is.null(CI)){
    if(CI == "LOD1") {
      LOD.pos <- c(x[[2]][1], x[[2]][2])
    } else {
      LOD.pos <- c(x[[3]][1], x[[3]][2])
    }
    
    LOD <- df[,2] >= LOD.pos[1] & df[,2] <= LOD.pos[2]
    
    qtls <- matrix(NA, nrow=4, ncol=dim(df)[1])
    qtls[1:4, which(LOD)] <- as.matrix(df.tot[which(df.tot[,1] == "QTL"),3:6])
    qtls <- t(qtls)
    colnames(qtls) <- paste0(colnames(df)[3:6])
    positions <- rep(NA, dim(qtls)[1])
    positions[which(LOD)] <- seq(LOD.pos[1], LOD.pos[2], length.out = sum(LOD))
    qtls <- data.frame(positions, qtls)
    qtls[,2:5] <- apply(qtls[,2:5], 2, as.character)
    qtls_melt <- melt(qtls, measure.vars = colnames(qtls)[-1])

    p <- p + geom_line(data = qtls_melt, aes(x = as.numeric(positions), 
                                             y = variable), na.rm = T, size = 1.5)
  }
  p
}


