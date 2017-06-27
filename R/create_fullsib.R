#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: create_fullsib.R                                              #
#                                                                     # 
# Contains:                                                           #
# create_fullsib                                                      #
# print_fullsib                                                       #
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
# Function: create_fullsib                                            #
#                                                                     #
# Receives all the necessary information to perform QTL mapping       #
# One inserts the map, phenotype(s) and unlinked markers (if there is)#
# It is also required to define the step to obtain the multipoint prob#
# as well error.function and map function (can choice between 3)      #
# It returns a object from class 'fullsib'                            #
#######################################################################

## -------------------------
## create_fullsib function

#' Creates the main object to perform QTL mapping on \pkg{fullsibQTL}
#' 
#' First function that will be used on this package. It receives the data file
#' containing all genotypes and phenotypes for all individuals on the mapping
#' population. Also, requires the genetic map obtained from \pkg{onemap}. On
#' this step, probabilities for QTL genotypes are obtained using Hidden Markov
#' Models technology.
#' 
#' In a fullsib cross, a priori we have four genetic classes segregating on the
#' progeny. However it is not always possible to identify four different
#' classes. This happens because some regions on the genetic map may be
#' saturated with molecular markers which are not fully informative. For
#' example, if a given linkage group only have markers segregating as 1:1
#' fashion, we are able to identify two differents genotipic classes, resulting
#' in only one genetic effect significative.  This can also happens with
#' markers segregating in 1:2:1 (only two of three effects are estimable) and
#' 3:1 fashion, as well some this markers. When this situation is verified,
#' colinearity on the QTL genotype probabilities is found and may lead some
#' problems on the analysis.
#' 
#' To avoid this kind of problems, we investigated for each position on the
#' genome, the number of different genetic classes and if the any relationship
#' between then. If four classes are avaliable, the model is runned as usual.
#' On the other hand, we infer the estimable effects and use just them on the
#' analysis. For details of how model works under colinearity, check
#' \pkg{fullsibQTL} Vignette.
#' 
#' The method print returns a text output summarizing the object of class
#' \emph{fullsib} or \emph{fullsib_cofactors}. Briefly, it indicates the total
#' number of linkage groups, markers and map total size. For each group, it is
#' indicated, its size, number of marker on the group and the markers
#' segregation found on the group. It is also showed the grid used to obtain
#' multipoint probabilities on the genome. If the object is also from
#' \emph{fullsib_cofactors}, it is indicated the trait that selection was done,
#' as well the markers location. If a selected marker is not in the map,
#' \code{NA} will be indicated on its corresponding linkage groups.
#' 
#' @aliases create_fullsib print.fullsib
#' 
#' @param input.obj an object of class \emph{outcross}, obtained with
#' \code{read.outcross} function from \pkg{onemap} package (>2.0-3), where it
#' is also possible to include phenotypic values.
#' @param map.list A map obtained with \pkg{onemap}, it is either an object of
#' class \emph{sequence} or list of objects from the same class.
#' @param step Defines the distance in centiMorgans (cM) between positions on
#' the linkage groups that genotype probabilities are obtained. By default, the
#' probabilities is calculated only at the marker locus \sQuote{step=0}. If
#' \sQuote{step=1}, for example, at every 1 cM a virtual marker will be created
#' to evaluate if there is assocition with phenotipic values.
#' @param error.prob Probability of genotyping error.
#' @param map.function Mapping function that will be used to obtain the
#' probabilities. One can choose between \emph{kosambi}, \emph{haldane} or
#' \emph{morgan}.
#' @param condIndex threshold used to evaluate the level of informativeness for
#' QTL genotypes probabilities. See details.
#' 
#' @return An object of class \emph{fullsib} is returned, which has the
#' following components:
#' 
#' \item{map}{List with objects of class \emph{sequence} obtained with
#' \pkg{onemap}. Each position of the list contains a different linkage group.}
#' 
#' \item{unlinked.mkr}{Vector contaning markers that are not placed in the
#' genetic map. The markers are represented by the indices as done with
#' \pkg{onemap} package.}
#' 
#' \item{pheno}{Matrix containing the phenotypes to be used for QTL mapping.
#' Rows represents individuals and columns, the traits.}
#' 
#' \item{probs}{List (of size equal to \code{map}) containing the multipoint
#' probabilities and additional required information. Each position of the list
#' has four components:
#' 
#' \sQuote{cond.prob}: an array containing the multipoint probabilities of
#' dimension [n.ind, n.locus, 4], where n.ind is the number of individuals and
#' n.locus the total of position (on a linkage group) where the probabilities
#' were considered
#' 
#' \sQuote{newmap}: list containing three components: \code{dist}, the cumulate
#' distance between all markers (including virtual markers); \code{phase},
#' linkage phases between them, \code{type}, markers segregation type.
#' 
#' \sQuote{step} is the walking grid used to estimate the probabilities along
#' the linkage group;
#' 
#' \sQuote{colin} has information concern on the locus informativeness, i.e,
#' \code{geno.class} indicates how many different genotypes means are estimable
#' (2, 3 or 4) and \code{type} indicates what kind the colinearity one may have
#' on these probs (0, 1, \ldots{}, 9). For further details for colinearity, see
#' the \pkg{fullsibQTL} Vignette.}
#' @author Rodrigo Gazaffi, \email{rgazaffi@@gmail.com}
#' @seealso 
#' \code{\link[onemap]{read_onemap}}
#' \code{\link[fullsibQTL]{example_QTLfullsib}}
#' 
#' @references
#' Wu, R., Ma, C.-X., Painter, I., Zeng, Z.-B. (2002a) Simultaneous maximum
#' likelihood estimation of linkage and linkage phases in outcrossing species.
#' \emph{Theoretical Population Biology} 61: 349-363
#' 
#' Wu, R., Ma, C.-X., Wu, S. S., Zeng, Z.-B. (2002b). Linkage mapping of
#' sex-specific differences. \emph{Genetical Research} 79: 85-96
#' @keywords IO
#' @examples
#'   data(example_QTLfullsib)
#' 
#'   ##fullsib is main object to perform QTL maping
#' 
#'   fsib <- create_fullsib(example_QTLfullsib,
#'                          list(LG1_final, LG2_final, LG3_final, LG4_final),
#'                          step=1,map.function="kosambi",condIndex=3.5)
#' 

create_fullsib <- function(input.obj, map.list, step=0, error.prob=1e-4,
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
  fullsib <- calc_probs(fullsib, step=step, error.prob=error.prob, map.function=map.function, condIndex)
  fullsib$mapped.mkrs <- colnames(input.obj$geno)[mapped.mkrs]
  class(fullsib) <- "fullsib"
  
  return(fullsib)
}

## -------------------------
## print.fullsib function

#######################################################################
#                                                                     #
# Function: print.fullsib                                             #
#                                                                     #
# Resume the main information on the object from fullsib class        #
#######################################################################

print.fullsib <- function(x, ...) {
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

  if (any(class(x) == "fullsib_cofactors")){
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
}#end print_fullsib