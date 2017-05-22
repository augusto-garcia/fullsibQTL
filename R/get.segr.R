#######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: get.segr.R                                                    #
#                                                                     # 
# Contains: get.segr                                                  #
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
# Function:   get.segr.R                                              #
#                                                                     # 
# Receives the object of class fullsib.char with QTL information and  #
# returns the QTL segregation                                         #
#                                                                     #
#######################################################################


get.segr <- function(fschar, probs1=.05, probs2=.05){

  if (!any(class(fschar) == "fullsib.char"))
    stop(sQuote(deparse(substitute(fschar))),
         " is not an object of class 'fullsib.char'")

  if((probs1 < 0 ) || (probs1 > 1))
    stop("probs1 argument should be a probability between 0 and 1")
  
  if((probs2 < 0 ) || (probs2 > 1))
    stop("probs2 argument should be a probability between 0 and 1")
  
  ##iddentify performed tests
  used.model <- fschar["model",]
  
  ##step1 in lod
  step1 <- fschar[c("LOD_H1", "LOD_H2", "LOD_H3"),]
  ##threshold probs1 em thr.lod
  thr.lod <- qchisq(probs1, 1, lower.tail=FALSE)/4.61
  
  ##step2 in pvalue
  step2 <- fschar[c("H4_pvalue", "H5_pvalue", "H6_pvalue"),]

  #######################################################
  ## model used when three genetics effects are estimable
  #######################################################

  if(used.model == 0){
    signif.lod <- step1 > thr.lod

    if (sum(signif.lod) == 0) # no signif ef.
      return(cat("No QTL is considered with probs1 = ", probs1, "\n"))
    
    if (sum(signif.lod) == 1) # only one signif ef.
      return(cat("QTL segregation is 1:1\n"))
    
    if (sum(signif.lod) == 2){ # two signf ef.
      signif.lod <- paste(which(signif.lod),collapse=":")
      switch(EXPR = signif.lod,
             '1:2' = {
               if(step2[1] < probs2)
                 return(cat("QTL segregation is 1:1:1:1\n"))
               else
                 return(cat("QTL segregation is 1:2:1\n"))
             },
             '1:3' = {
               if(step2[2] < probs2)
                 return(cat("QTL segregation is 1:1:1:1\n"))
               else
                 return(cat("QTL segregation is 1:2:1\n"))
             },
             '2:3' = {
               if(step2[3] < probs2)
                 return(cat("QTL segregation is 1:1:1:1\n"))
               else
                 return(cat("QTL segregation is 1:2:1\n"))
               
             })#end-switch
    }

    if (sum(signif.lod) == 3){ # three signf ef.
      signif.pvalue <- as.character(sum(step2 > probs2))
      switch(EXPR = signif.pvalue,
             '3' = {return(cat("QTL segregation is 3:1\n"))},
             '1' = {return(cat("QTL segregation is 1:2:1\n"))},
             {return(cat("QTL segregation is 1:1:1:1\n"))})
    }
  } #end - used.model==0

  #######################################################
  ## model used when just one additive effect is estimable
  #######################################################
  
  if(used.model == 1){ ## just ap is estimable
    if(step1[1] > thr.lod)
      return(cat("QTL segregation is 1:1\n"))
    else
      return(cat("No QTL is considered with probs1 = ", probs1, "\n"))
  }

  if(used.model == 2){ ## just aq is estimable
    if(step1[2] > thr.lod)
      return(cat("QTL segregation is 1:1\n"))
    else
      return(cat("No QTL is considered with probs1 = ", probs1, "\n"))
  }
  
  #######################################################
  ## model used when "two" additive effect is estimable
  #######################################################

  if((used.model == 3) || (used.model == 4)){ ##1 and (4 or 5)
    signif.lod <- step1[1:2] > thr.lod
    if(sum(signif.lod, na.rm=T) == 0) return(cat("No QTL is considered with probs1 =", probs1, "\n"))
    if(sum(signif.lod, na.rm=T) == 1){
      if(signif.lod[1] == TRUE) return(cat("QTL segregation is 1:1\n"))
      if(signif.lod[2] == TRUE) return(cat("QTL segregation is 1:2:1\n"))
    }
    if(sum(signif.lod, na.rm=T) == 2){
      if(step2[1] < probs2)
        return(cat("QTL segregation is 1:2:1\n"))
      else
        return(cat("QTL segregation is 3:1\n"))
    }
  }

    if((used.model == 5) || (used.model == 6)){ ## (6 or 7) and 2
    signif.lod <- step1[1:2] > thr.lod
    if(sum(signif.lod, na.rm=T) == 0) return(cat("No QTL is considered with probs1 =", probs1, "\n"))
    if(sum(signif.lod, na.rm=T) == 1){
      if(signif.lod[1] == TRUE) return(cat("QTL segregation is 1:2:1\n"))
      if(signif.lod[2] == TRUE) return(cat("QTL segregation is 1:1\n"))
    }
    if(sum(signif.lod, na.rm=T) == 2){
      if(step2[1] < probs2)
        return(cat("QTL segregation is 1:2:1\n"))
      else
        return(cat("QTL segregation is 3:1\n"))
    }
  }
  
  if((used.model == 7) || (used.model == 8)){ ##(8 or 9)  and 3
    signif.lod <- step1[c(1,3)] > thr.lod
    if(sum(signif.lod, na.rm=T) == 0) return(cat("No QTL is considered with probs1 =", probs1, "\n"))
    if(sum(signif.lod, na.rm=T) == 1){
      if(signif.lod[1] == TRUE) return(cat("QTL segregation is 1:2:1\n"))
      if(signif.lod[2] == TRUE) return(cat("QTL segregation is 1:1\n"))
    }
    if(sum(signif.lod, na.rm=T) == 2){
      if(step2[2] < probs2)
        return(cat("QTL segregation is 1:2:1\n"))
      else
        return(cat("QTL segregation is 3:1\n"))
    }              
  }
  
  #######################################################
  ## model used when LG have mkr with 3:1 fashion in coupling
  #######################################################

  if(used.model > 8) {
    if(step1[1] > thr.lod)
      return(cat("QTL segregation is 3:1\n"))
    else
      return(cat("No QTL is considered with probs1 =", probs1, "\n"))
  }
  
}
