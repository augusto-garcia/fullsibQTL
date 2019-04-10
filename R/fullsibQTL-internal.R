######################################################################
#                                                                     #
# Package: fullsibQTL                                                 #
#                                                                     #
# File: fullsibQTL-internal.R                                         #
#                                                                     #
# Contains:                                                           #
# calc_probs*          used in create.fullsib                         #
# calc_genoprobs:*     used in calc_probs                             #
# colin_check**        used in calc_probs                             #
# geno_class_identify  used in calc_probs                             #
# model_select_lg:     used in cim.char                               #
# model_select:        used in cim.scan                               #
# prob2cof:            used in cof.def..., cof.sel...  and r2.model   #
# unlink2cof:          used in cof_definition and cof_selection       #
# get_phase            used in draw.phase                             #
# viewer               used in plot_fullsibQTL                        #
#                                                                     #
# Written by Rodrigo Gazaffi                                          #
# copyright (c) 2011, Rodrigo Gazaffi                                 #
# *written by Marcelo Mollinari and Rodrigo Gazaffi                   #
# **written by Gabriel R A Margarido and Rodrigo Gazaffi              #
#                                                                     #
# First version: 09/30/2011                                           #
# Last  version: 09/30/2011                                           #
# License: GPL-3                                                      #
#                                                                     #
#######################################################################

#######################################################################
#                                                                     #
# Function: calc_probs                                                #
#                                                                     #
# Function to obtain multipoint QTL probability for all chromossome   #
# it is used in create.fullsib                                        #
#######################################################################

calc_probs <- function(fullsib, step, error.prob, map.function, condIndex)
{
  ##show the time during calculation of probs...
  pb <- txtProgressBar(style=3)
  setTxtProgressBar(pb, 0)

  for (i in 1:length(fullsib$map)){

    setTxtProgressBar(pb, i/length(fullsib$map))
    
    seq.mark <- fullsib$map[[i]]$seq.num
    dist.mark <- cumsum(c(0,       kosambi(fullsib$map[[i]]$seq.rf)))
    obj.name <- fullsib$map[[1]]$data.name
    #print(obj.name)
    #print(get(obj.name))
    names(dist.mark) <- colnames(get(obj.name)$geno)[seq.mark]

    fullsib$probs[[i]] <- calc_genoprobs(dist=dist.mark,
                     geno=get(obj.name)$geno[,seq.mark],
                     type=get(obj.name)$segr.type.num[seq.mark],
                     phase=fullsib$map[[i]]$seq.phases,
                     step = step, error.prob=error.prob,
                     map.function=map.function)     
    fullsib$probs[[i]]$step <- step
    
    if(is.null(fullsib$probs[[i]]$cond.prob))
      next
    else{
      tmp.colin <- tmp.geno <-
        rep(NA, length(fullsib$probs[[i]]$cond.prob[1,,1]))
      for (j in 1:length(fullsib$probs[[i]]$cond.prob[1,,1])){
        probs2check <- fullsib$probs[[i]]$cond.prob[,j,]
        tmp.colin[j] <- colin_check(probs2check, condIndex)
        tmp.geno[j] <- geno_class_identify(tmp.colin[j])
      }
      fullsib$probs[[i]]$colin$type <- tmp.colin
      fullsib$probs[[i]]$colin$geno.class <- tmp.geno
    }#end if
  }#end for-i
  close(pb)
  fullsib
}


#######################################################################
#                                                                     #
# Function: calc_genoprobs                                            #
#                                                                     #
# Function to obtain multipoint QTL probability for each chromossome  #
# Probabilites are obtained by HMM approach similiar Wu et al. 2002   # 
# This is done using C codes - this is already inside Onemap          #
#######################################################################

create_map <-function(dist,phase,type,step)    {
    ## function that creates a map to be used in calc_genoprobs
    nam<-names(dist)
    if(step<0) stop("step must be >= 0.")
    minloc <- min(dist)
    dist <- dist-minloc
    if(step==0) return(list(dist=dist+minloc, type=type, phase=phase))
    ##no 'loc' markers
    else{
        a <- seq(floor(min(dist)),max(dist),by = step)
        #print(a)
        #print(dist)
        #distn <- names(dist)
        #print(as.integer(dist))
        
        check.na <- match(as.integer(a),as.integer(dist))
        ##to make sure it works well  - da problema com decimais
        
        #print(check.na)
        #print(is.na(check.na))
        #print(which(is.na(a)))
        #print(a)
        #print(dist)
        #print(check.na)       
        #print(is.na(check.na))
        #print(a[is.na(check.na)])
        
        #a <- a[is.na(match(a,dist))]
        a <- a[is.na(check.na)]
        #print(a)
        if(length(a) == 0){
            ##huge step and no 'loc' was obtained - similar step==0 
            return(list(dist=dist+minloc, type=type, phase=phase))
        }
        else{
            ## step used to create 'loc' markers
            names(a) <- paste("loc",a,sep="")
            dist.temp<-sort(c(a,dist))+minloc
            index<-match(nam,names(dist.temp))
            type.temp<-rep(1,length(dist.temp))
            type.temp[index]<-type
            phase.temp<-rep(1,length(dist.temp))
            phase.temp[index]<-c(phase,1)
            #print(dist.temp)
            return(list(dist=dist.temp,
                        type=type.temp,
                        phase=phase.temp))
        }
    }
}

calc_genoprobs <- function(dist, geno, type, phase, step,
                           error.prob, map.function)
{


    ## map function
    mf.k <- function(d) 0.5*tanh(d/50)
    mf.h <- function(d) 0.5*(1-exp(-d/50))
    mf.m <- function(d) sapply(d,function(a) min(a/100,0.5))
    
    if(map.function=="kosambi") mf <- mf.k
    else if(map.function=="morgan") mf <- mf.m
    else mf <- mf.h
    
    
    ## don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 1) {
        error.prob <- 1-1e-50
        warning("error.prob shouldn't be > 1!")
    }
    
    n.ind <- nrow(geno)
    n.mar <- ncol(geno)
    
    ## new map and recombination fractions
    map <- create_map(dist,phase,type,step)
    rf <- mf(diff(map$dist))
    
    ## new genotype matrix with pseudomarkers filled in
    newgen <- matrix(ncol=length(map$dist),nrow=nrow(geno))
    dimnames(newgen) <- list(NULL,names(map$dist))
    newgen[,colnames(geno)] <- geno 
    newgen[is.na(newgen)] <- 0
    n.pos <- ncol(newgen)
    marnames <- names(map$dist)
    
    
    ##this is already inside Onemap
    #z<-.C("onemap:::calc_genoprob_outbred",
    ## z<-.C("calc_genoprob_outbred",
    z<-.C("onemap_calc_genoprob_outbred",
        as.integer(n.ind),
        as.integer(n.pos),
        as.integer(map$type),
        as.integer(map$phase), 
        as.integer(newgen),
        as.double(rf),   
        as.double(error.prob),
        genoprob=as.double(rep(0,4*n.ind*n.pos)))
        #genoprob=as.double(rep(0,4*n.ind*n.pos)))
        #PACKAGE="onemap")
    
    ## re-arrange marginal probabilites
    z2 <- array(z$genoprob,dim=c(n.ind,n.pos,4))
    dimnames(z2) <- list(NULL, marnames, c("P1Q1","P1Q2","P2Q1","P2Q2"))
    
    if(length(head(map$phase,-1)) == 0){
        return(list(cond.prob=z2, newmap=map))
    }#two position in this group
    
    else{
        ##more then 2 positions
        link.phases <- matrix(NA,length(map$phase),2)
        link.phases[1,] <- rep(1,2) 
        for (i in 1:(length(map$phase)-1)) {
            switch(EXPR=map$phase[i],
                   link.phases[i+1,] <- link.phases[i,]*c(1,1), 
                   link.phases[i+1,] <- link.phases[i,]*c(1,-1),
                   link.phases[i+1,] <- link.phases[i,]*c(-1,1),
                   link.phases[i+1,] <- link.phases[i,]*c(-1,-1))
        }
        link.phases <- apply(link.phases,1,
                             function(x) paste(as.character(x),collapse="."))
        ord.phase.num<-rep(NA, length(link.phases))
        ord.phase.num[which(link.phases=="1.1")]<-1
        ord.phase.num[which(link.phases=="1.-1")]<-2
        ord.phase.num[which(link.phases=="-1.1")]<-3
        ord.phase.num[which(link.phases=="-1.-1")]<-4
        ord.phase.mat<-matrix(c(1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1),4,4)
        
        for(i in 1:n.ind){
            for(j in 1:length(ord.phase.num)){
                z2[i,j,]<-z2[i,j,ord.phase.mat[,ord.phase.num[j]]]     
            }
        }
        return(list(cond.prob=z2, newmap=map))
    }
}#end function

#######################################################################
#                                                                     #
# Function: calc_genoprobs                                            #
#                                                                     #
# Linkage groups with differents segregation markers may differ about #
# the information on the probabilities. So it is necessary to identify#
# if this happens and then check what kind of QTL model we can search #
#                                                                     #
# All possible contrasts are indicated here                           #
# contrast  <- matrix(c(1, 1,-1,-1,         # 1                       #
#                       1,-1, 1,-1,         # 2                       #
#                       1,-1,-1, 1,         # 3                       #
#                       0, 0, 1,-1,         # 4                       #
#                       1,-1, 0, 0,         # 5                       #
#                       0, 1, 0,-1,         # 6                       #
#                       1, 0,-1, 0,         # 7                       #
#                       0, 1,-1, 0,         # 8                       #
#                       1, 0, 0,-1,         # 9                       #
#                       1/3, 1/3, 1/3,  -1, #10                       #
#                       1/3, 1/3,  -1, 1/3, #11                       #
#                       1/3,  -1, 1/3, 1/3, #12                       #
#                        -1, 1/3, 1/3, 1/3) #13                       #
#                                                                     #
# If the result is 0 to 12 is the model we need to used               #
# colin check output for contrasts showed previously                  #
# 0: contr 1, 2 e 3  (original model)                                 #
# 1: contr 1         (similar a backcross)                            #
# 2: contr 2         (similar a backcross)                            #
# 3: contr 1 e 4                                                      #
# 4: contr 1 e 5                                                      #
# 5: contr 6 e 2                                                      #
# 6: contr 7 e 2                                                      #
# 7: contr 8 e 3                                                      #
# 8: contr 9 e 3    (similar a F2)                                    #
# 9: contr 10                                                         #
#10: contr 11                                                         #
#11: contr 12                                                         #
#12: contr 13                                                         #
#######################################################################
                                                  
colin_check <- function (probs_1pos, condIndex)
{
  decomp <- svd(probs_1pos)
  conditionIndex <- max(decomp$d)/decomp$d

  ## For our datasets 3.5 is working very well
  if (sum(conditionIndex > condIndex) == 1) {
    V <- decomp$v
    phi <- t(apply(V,1,function(y) y^2/decomp$d^2))
    phi_k <- apply(phi,1,sum)
    
    Pi <- t(apply(phi,2,function(y) y/phi_k))
    probs.similar <- sort(tail(order(Pi[4,]),2));
    switch(EXPR=paste(probs.similar,collapse=":"),
           "1:2" = {contr.type <- 3},
           "3:4" = {contr.type <- 4},
           "1:3" = {contr.type <- 5},
           "2:4" = {contr.type <- 6},
           "1:4" = {contr.type <- 7},
           "2:3" = {contr.type <- 8},
           )
  }
  ## For our datasets 3.5 is working very well
  else if (sum(conditionIndex > condIndex) == 2) {
    V <- decomp$v;
    phi <- t(apply(V,1,function(y) y^2/decomp$d^2));
    phi_k <- apply(phi,1,sum);
    Pi <- t(apply(phi,2,function(y) y/phi_k));
    probs.similar1 <- sort(tail(order(Pi[3,]),2))
    probs.similar2 <- sort(tail(order(Pi[4,]),2))
    switch(EXPR=paste(c(probs.similar1,probs.similar2),collapse=":"),
           "1:2:3:4" = ,
           "3:4:1:2" = {contr.type <- 1},
           "1:3:2:4" = ,
           "2:4:1:3" = {contr.type <- 2},
           {contr.type <- 10}
           )
}
  else
    contr.type <- 0
  contr.type
}#end function

#######################################################################
#                                                                     #
# Function: geno_class_identify                                       #
#                                                                     #
# This function returns the number of estimable means in a single cM  #
# under the colinearity relationship that is obtained with colin_check#
# This information is passed as argument for EM algorithm done in C   #
#######################################################################

geno_class_identify <- function(colin_number)
{  
  if (colin_number == 0) 
    geno.class <- 4
  ##no colinearity - probs fully informative
  else if (colin_number > 2 && colin_number < 9)
    geno.class <- 3
  ## 3,4,5,6,7,8 - has 3 estimable means like "B" mkr
  else 
    geno.class <- 2
  ##1,2 (D's) e 9,10,11,12 (C's)
  geno.class
}#end function


#######################################################################
#                                                                     #
# Function: model_select_lg                                           #
#                                                                     #
# For CIM analysis cofactors are included on the model, obtained with #
# either cof_selection or cof_definition function. However, depending #
# on the position (defined the window size) some cofactors are removed#
# from the analysis. So, this function indicated which cofactors are  #
# dropped from the model and for which position this is done          #
#######################################################################

model_select_lg <- function(fullsib, ws, LG){

  ##is there any cof in this group?
  if.cof <- which(fullsib$cofactors$names.cof$lg == LG)
  
  if(length(if.cof) != 0){ #if there is cofactor
    ##which are the positions that REAL MARKERS are placed (location)
    mkr_lg <- which(substring(names(fullsib$probs[[LG]]$newmap$dist),1,3) != "loc") 
    mkr_lg <- fullsib$probs[[LG]]$newmap$dist[mkr_lg]

    ##which are the positions that COFACTORS are placed (location)
    cof.num <- which(names(fullsib$probs[[LG]]$newmap$dist) %in%
                     as.character(fullsib$cofactors$names.cof$mkr))
    cof.cm <- fullsib$probs[[LG]]$newmap$dist[cof.num]

    ##For cof.cm i will get the longest distance: ws_1 or ws_2
    
    ##finding adjacent markers: part1
    wsl1 <- which(mkr_lg %in%cof.cm == TRUE) - 1
    wsl1 <- ifelse(wsl1 < 1, yes = 1, no = wsl1)
    
    wsr1 <- which(mkr_lg %in%cof.cm == TRUE) + 1
    wsr1 <-  ifelse(wsr1 > length(mkr_lg),
                    yes =  length(mkr_lg), no = wsr1)
    
    ##finding the distance limits for cofactor part2
    wsl2 <- cof.cm - ws
    wsr2 <- cof.cm + ws
    
    wsl2 <- sapply(wsl2, function(x) ifelse(sum(cumsum(mkr_lg > x)) == 0, yes = length(mkr_lg), no = which(cumsum(mkr_lg > x) == 1) - 1))
    wsl2 <- ifelse(wsl2 < 1, yes = 1, no = wsl2)
    wsr2 <- sapply(wsr2, function(x) ifelse(sum(cumsum(mkr_lg > x)) == 0, yes = length(mkr_lg), no = which(cumsum(mkr_lg > x) == 1)))
    wsr2 <- ifelse(wsr2 > length(mkr_lg), yes =  length(mkr_lg), no = wsr2)
    
    ##which parts gives a more extense interval
    wslfinal <- apply(cbind(wsl1, wsl2), 1, min)
    wsrfinal <- apply(cbind(wsr1, wsr2), 1, max)

    ##define limits for window size: 
    wsfinal <- cbind(wslfinal, wsrfinal)   
    cim.window.limits <- function(x) {
      xres <- which(names(fullsib$probs[[LG]]$newmap$dist) %in%
                    names(mkr_lg[x]) );
      xres[2] <- ifelse (xres[2] == length(fullsib$probs[[LG]]$newmap$dist), yes = xres[2], no = xres[2] - 1)
      xres
    }    
    wsfinal.npos <- apply(wsfinal, 1, cim.window.limits)
    ##now wsfinal.npos have the positions that one marker should be dropped

    ##next step is to check for intersections in those positions,
    ##or if there are regions that should drop more than one marker
    n2drop <- matrix(0, ncol=length(fullsib$probs[[LG]]$newmap$dist), nrow=1)
    n2drop <- apply(wsfinal.npos, 2, function(x) {n2drop[ seq(x[1], x[2], 1)] <- 1;n2drop})
    all.models <- apply(n2drop, 1, paste, collapse="")
    tab.models <- names(table(all.models))

    ##positions that share the same cofactor to be removed
    out.cof <- sapply(tab.models, function(x) which(all.models == x),
                      simplify = FALSE)
    out.cof.end <- list()   
    for(j in 1:length(tab.models)){
      out.cof.end[[j]] <- vector("list",2)
      names(out.cof.end[[j]]) <- c("mkr2rm", "position")
      if(tab.models[j] == paste(rep("0", length(colnames(n2drop))),collapse="")){
        out.cof.end[[j]]$mkr2rm <- "none" 
        out.cof.end[[j]]$position <-  out.cof[[as.character(tab.models[j])]]
      }
      else{
        out.cof.end[[j]]$mkr2rm <- colnames(n2drop)[as.logical(as.integer(strsplit(tab.models[j], split="", useBytes=TRUE)[[1]]))]
        out.cof.end[[j]]$position <- out.cof[[as.character(tab.models[j]) ]]
      }     
    }#end for
  }#end if
  else{#if there is no cofactor
    out.cof.end <- list()   
    out.cof.end[[1]] <- vector("list",2)
    names(out.cof.end[[1]]) <- c("mkr2rm", "position")
    out.cof.end[[1]]$mkr2rm <- "none"
    out.cof.end[[1]]$position <- seq(1,length(fullsib$probs[[LG]]$newmap$dist), 1)
  }
  out.cof.end
}#end function


#######################################################################
#                                                                     #
# Function: model_select                                              #
#                                                                     #
# This function returns a list, indication for each group (chr or LG) #
# what cofactors should be used for CIM analysis and in which position#
# this is just model_select_lg, used for all groups                   #
#######################################################################

model_select <- function(fullsib, ws){
  cof.output <- vector("list", length(fullsib$map))
  names(cof.output) <- paste("LG", 1: length(fullsib$map), sep="-")
  for (LG in 1:length(fullsib$map))
    cof.output[[LG]] <- model_select_lg(fullsib, ws, LG)
  cof.output
}#end function


#######################################################################
#                                                                     #
# Function: prob2cof                                                  #
#                                                                     #
# This function considers a probability for one position and transform#
# in genetic predictor (based on colinearity evidence). The predictors#
# are used are cofactor for CIM analysis.                             #
#######################################################################

prob2cof <- function(prob, colin_type){
  colin_type <- as.character(colin_type)
  switch(EXPR=colin_type,
         "0" = prob %*% matrix(c(1, 1,-1,-1,  1,-1, 1,-1,  1,-1,-1,1),ncol=3),
         "1" = prob %*% matrix(c(1, 1,-1,-1), ncol=1),
         "2" = prob %*% matrix(c(1,-1, 1,-1), ncol=1),
         "3" = prob %*% matrix(c(1, 1,-1,-1,  0, 0, 1,-1),ncol=2),
         "4" = prob %*% matrix(c(1, 1,-1,-1,  1,-1, 0, 0),ncol=2),
         "5" = prob %*% matrix(c(0, 1, 0,-1,  1,-1, 1,-1),ncol=2),
         "6" = prob %*% matrix(c(1, 0,-1, 0,  1,-1, 1,-1),ncol=2),
         "7" = prob %*% matrix(c(0, 1,-1, 0,  1,-1,-1, 1),ncol=2),
         "8" = prob %*% matrix(c(1, 0, 0,-1,  1,-1,-1, 1),ncol=2),
         "9" = prob %*% matrix(c(1/3, 1/3, 1/3,  -1),ncol=1),
         "10"= prob %*% matrix(c(1/3, 1/3,  -1, 1/3),ncol=1),
         "11"= prob %*% matrix(c(1/3,  -1, 1/3, 1/3),ncol=1),
         "12"= prob %*% matrix(c( -1, 1/3, 1/3, 1/3),ncol=1)
         ) #end-switch
}#end function

#######################################################################
#                                                                     #
# Function: unlink2cof                                                #
#                                                                     #
# This function allows using unlinked markers as cofactors            #
#######################################################################
unlink2cof <- function(mk.unmap, obj.name, pheno.index){ 
  cof.unmap <- vector("list", length(mk.unmap))
  for (j in 1:length(mk.unmap)){
    #for each markers: see the segregation
    segr <- get(obj.name)$segr.type.num[ mk.unmap[j] ]
    if(segr == 1){
      ##if markers A (1:1:1:1)
      cof.unmap[[j]] <- matrix(0, nrow(get(obj.name)$geno),3)    
      geno1 <- which(get(obj.name)$geno[,mk.unmap[j]] ==1)
      geno2 <- which(get(obj.name)$geno[,mk.unmap[j]] ==2)
      geno3 <- which(get(obj.name)$geno[,mk.unmap[j]] ==3)
      geno4 <- which(get(obj.name)$geno[,mk.unmap[j]] ==4)    
      if(length(geno1) > 0)
        cof.unmap[[j]][geno1, ] <- t(apply(cof.unmap[[j]][geno1, ], 1, function(x) x <- c( 1, 1, 1)))
      if(length(geno2) > 0)
        cof.unmap[[j]][geno2, ] <- t(apply(cof.unmap[[j]][geno2, ], 1, function(x) x <- c( 1,-1,-1)))
      if(length(geno3) > 0)
        cof.unmap[[j]][geno3, ] <- t(apply(cof.unmap[[j]][geno3, ], 1, function(x) x <- c(-1, 1,-1)))
      if(length(geno4) > 0)
        cof.unmap[[j]][geno4, ] <- t(apply(cof.unmap[[j]][geno4, ], 1, function(x) x <- c(-1,-1, 1)))
    }
    if(segr==2){
      ##if markers B1 (1:2:1)
      cof.unmap[[j]] <- matrix(0, nrow(get(obj.name)$geno),2)       
      geno1 <- which(get(obj.name)$geno[,mk.unmap[j]] ==1)
      geno2 <- which(get(obj.name)$geno[,mk.unmap[j]] ==2)
      geno3 <- which(get(obj.name)$geno[,mk.unmap[j]] ==3)
      if(length(geno1) > 0)
        cof.unmap[[j]][geno1, ] <- t(apply(cof.unmap[[j]][geno1, ], 1, function(x) x <- c( 1, 0)))            
      if(length(geno2) > 0)
        cof.unmap[[j]][geno2, ] <- t(apply(cof.unmap[[j]][geno2, ], 1, function(x) x <- c(-1, 1)))
      if(length(geno3) > 0)
        cof.unmap[[j]][geno3, ] <- t(apply(cof.unmap[[j]][geno3, ], 1, function(x) x <- c(-1,-1)))
    }
    if(segr == 3){
      ##if markers B2 (1:2:1)
      cof.unmap[[j]] <- matrix(0, nrow(get(obj.name)$geno),2)     
      geno1 <- which(get(obj.name)$geno[,mk.unmap[j]] ==1)
      geno2 <- which(get(obj.name)$geno[,mk.unmap[j]] ==2)
      geno3 <- which(get(obj.name)$geno[,mk.unmap[j]] ==3)
      if(length(geno1) > 0)
        cof.unmap[[j]][geno1, ] <- t(apply(cof.unmap[[j]][geno1, ], 1, function(x) x <- c( 0, 1)))
      if(length(geno2) > 0)
        cof.unmap[[j]][geno2, ] <- t(apply(cof.unmap[[j]][geno2, ], 1, function(x) x <- c( 1,-1)))
      if(length(geno3) > 0)
        cof.unmap[[j]][geno3, ] <- t(apply(cof.unmap[[j]][geno3, ], 1, function(x) x <- c(-1,-1)))
    }
    if(segr == 4){
      ##if markers B3 (1:2:1)
      cof.unmap[[j]] <- matrix(0, nrow(get(obj.name)$geno),2)   
      geno1 <- which(get(obj.name)$geno[,mk.unmap[j]] ==1)
      geno2 <- which(get(obj.name)$geno[,mk.unmap[j]] ==2)
      geno3 <- which(get(obj.name)$geno[,mk.unmap[j]] ==3)
      if(length(geno1) > 0)
        cof.unmap[[j]][geno1, ] <-t( apply(cof.unmap[[j]][geno1, ], 1, function(x) x <- c( 1, 1)))
      if(length(geno2) > 0)
        cof.unmap[[j]][geno2, ] <-t( apply(cof.unmap[[j]][geno2, ], 1, function(x) x <- c( 0,-1)))
      if(length(geno3) > 0)
        cof.unmap[[j]][geno3, ] <-t( apply(cof.unmap[[j]][geno3, ], 1, function(x) x <- c(-1, 1)))
    }
    if(segr > 4){
      ##if segr C (5), D1 (6) or D2 (7)
      cof.unmap[[j]] <- as.matrix(get(obj.name)$geno[,mk.unmap[j] ])
      cof.unmap[[j]][which( cof.unmap[[j]] == 2)] <- -1
    }
    cof.unmap[[j]] <- as.matrix(cof.unmap[[j]][pheno.index, ])
  }
  names(cof.unmap) <- paste("mkr", mk.unmap, "cof", sep=".")
  cof.unmap
}#end function



#######################################################################
#                                                                     #
# Function: get_phase                                                 #
#                                                                     #
# This function draw the QTL configuration in the linkage group       #
#######################################################################

####################
get_phase <- function(model, alpha.p, lod.ap, alpha.q, lod.aq, probs){
  result <- rep(NA,4)
  model <- as.character(model)


  thr <- qchisq(probs,1, lower.tail=FALSE)   ##p.value to LRT
  thr <- thr / 4.61  ## LRT to LOD

  switch(EXPR=model,
         '1' = { #D1
           if(lod.ap > thr){#if ap is **
             if(alpha.p > 0)
               result[1:2] <- c("P1", "P2")
             else
               result[1:2] <- c("P2", "P1")
           }
           else
             result[1:2] <- c("P0", "P0")
           result[3:4] <- c("Q0", "Q0")
         },
         
         '2' = { #D2
           result[1:2] <- c("P0", "P0")
           if(lod.aq > thr){#if aq is **
             if(alpha.q > 0)
               result[3:4] <- c("Q1", "Q2")
             else
               result[3:4] <- c("Q2", "Q1")
           }
           else
             result[3:4] <- c("Q0", "Q0")

         },
         
         '7' = {
           ## it means the contrasts are:
           ## 0  1 -1  0 (8) 'it looks like imprinting  contrasts'
           ## 1 -1 -1  1 (3) original contrast for dominance
           ## for me this situation we can not infer linkage phase
           result[1:4] <- c("P0", "P0", "Q0","Q0")          
         },
         
         '8' =, '9'= { #like F2

           if(lod.ap > thr){#if ap is **
             if(alpha.p > 0)
               result[1:4] <- c("P1", "P2", "Q1", "Q2")
             else
               result[1:4] <- c("P2", "P1", "Q2", "Q1")
           }
           else
             result[1:4] <- c("P0","P0", "Q0","Q0")
         },
         
         '0' =, '3'=, '4'=,'5'=, '6'= {

           ##phase for parent P
           if(lod.ap > thr){#if ap is **
             if(alpha.p > 0)
               result[1:2] <- c("P1", "P2")
             else
               result[1:2] <- c("P2", "P1")
           }
           else
             result[1:2] <- c("P0","P0")
           ##phase for parent Q
           if(lod.aq > thr){#if aq is **
             if(alpha.q > 0)
               result[3:4] <- c("Q1", "Q2")
             else
               result[3:4] <- c("Q2", "Q1")
           }
           else
             result[3:4] <- c("Q0","Q0")
         },
         '10' = {
             if(lod.ap > thr){#if ap is **
                 if(alpha.p > 0)
                     result[1:4] <- c("P1", "P2", "Q1", "Q2")
                 else
                     result[1:4] <- c("P2", "P1", "Q2", "Q1")
             }
             else
                 result[1:4] <- c("P0", "P0", "Q0", "Q0")
         },
                  
         {
           warning("\n Should not be here! default options from switch in get_phase was used\n")
         })#end-switch
  result
}


draw_map2 <- function(map.list, horizontal = FALSE, names = FALSE, grid = FALSE,
                       cex.mrk = 1, cex.grp = 0.75, trait) {
  
  if (!any(class(map.list) == "list" | class(map.list) == "sequence")) 
    stop(deparse(substitute(map.list)),
         " is not an object of class 'list' or 'sequence'")
  
  if (class(map.list) == "sequence") 
    map.list <- list(map.list)
  j <- 1
  out <- data.frame()
  pos <- NULL
  marker <- NULL
  for (i in length(map.list):1) {
    if (!any(class(map.list[[i]]) == "sequence")) 
      stop("Object ", i, " in map.list is not an object of class 'sequence'")
    if (is.null(map.list[[i]]$seq.like)) 
      stop("Parameters are not estimated for object ", 
           i, " in map.list")
    #map <- cumsum(c(0, get(.map.fun)(map.list[[i]]$seq.rf)))
    map <- cumsum(c(0, kosambi(map.list[[i]]$seq.rf)))
    marnames <- colnames(get(map.list[[i]]$data.name, pos = 1)$geno)[map.list[[i]]$seq.num]
    out <- rbind(out, data.frame(dist = map, pos = j, marker = marnames))
    j <- j + 1
  }
  x <- tapply(out$dist, out$pos, max)
  y <- unlist(unique(out[2]))
  out.fake <- data.frame(dist = rep(c(0, max(out$dist)), max(y) + 
                           2), pos = c(0:(max(y) + 1)))
  if (horizontal == TRUE) {
    plot(out.fake, axes = FALSE, col = 0, xlab = "Distance (cM)",  ylab = "",
         main = paste("Genetic Map with cofactors for CIM analysis, \n considering pheno.col = ", trait, sep=" "), cex.main=.9)

    points(out[1:2], pch = "|", cex = cex.mrk, xlab = "Distance (cM)", 
           ylab = "", main = "Genetic Map")
    axis(1, at = seq(from = 0, to = 10 * round(max(x)/10), 
              by = 10), labels = seq(from = 0, to = 10 * round(max(x)/10), 
                          by = 10), cex.axis = 0.75)
    axis(2, y, paste("Group", rev(y)), lwd = 0, las = 2, 
         cex.axis = cex.grp)
    if (grid == TRUE) 
      abline(v = seq(from = 0, to = 10 * round(max(x)/10), 
               by = 10), lty = 2, lwd = 0.5, col = 2)
    for (i in y) {
      if (names == TRUE) 
        text(x = unlist(subset(out, pos == i, dist)), 
             y = i + max(y)/80,
             labels = unlist(subset(out, pos == i, marker)),
             srt = 90, cex = cex.mrk * 0.75, adj = c(0, 0.5))
      lines(x = c(0, x[i]), y = c(y[i], y[i]))
    }
  }
  else {
    plot(-out.fake[2:1], axes = FALSE, col = 0, ylab = "Distance (cM)", xlab = "",
         main = paste("Genetic Map with cofactors for CIM analysis, \n considering pheno.col = ", trait, sep=" "), cex.main=.9)

    points(-out[2:1], pch = 95, cex = cex.mrk)
    axis(2, cex.axis = 0.75,
         at = -seq(from = 0, to = 10 * round(max(x)/10), by = 10),
         labels = seq(from = 0, to = 10 * round(max(x)/10), by = 10), las = 2)
    
    axis(1, -y, paste("Group", rev(y)), lwd = 0, las = 2, 
         cex.axis = cex.grp)
    if (grid == TRUE) 
      abline(h = -seq(from = 0, to = 10 * round(max(x)/10), 
                      by = 10), lty = 2, lwd = 0.5, col = 2)
    for (i in y) {
      if (names == TRUE) 
        text(x = -i + max(y)/100,
             y = -unlist(subset(out, pos == i, dist)),
             labels = unlist(subset(out, pos == i, marker)),
             cex = cex.mrk * 0.75, adj = c(0, 0.5))
      lines(y = c(-0.2, -x[i] + 0.2), x = c(-y[i], -y[i]))
    }
  }
}

#######################################################################
#                                                                     #
# Function: viewer                                                    #
#                                                                     #
# Gets the use viewer to plot interactive graphic (plot_fulssibQTL.R) #
# from rstudioapi::viewer?? help                                      #
#######################################################################

#viewer <- function(...){
#  getOption("viewer")
#}
  
                                           
                                           
## return_geno function originally written for onemap package                                           
 #######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: return_geno.R                                                 #
# Contains: return_geno                                               #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2009, Gabriel R A Margarido                           #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# Function to create diplotypes based on segregation type and linkage phase
return_geno <-
function(segr.type, link.phases) {
  switch(EXPR=segr.type,
         'A.1' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","c","d")),
                  '1.-1'  = return(c("a","b","d","c")),
                  '-1.1'  = return(c("b","a","c","d")),
                  '-1.-1' = return(c("b","a","d","c"))
                  )
         },
         'A.2' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","a","c")),
                  '1.-1'  = return(c("a","b","c","a")),
                  '-1.1'  = return(c("b","a","a","c")),
                  '-1.-1' = return(c("b","a","c","a"))
                  )
         },
         'A.3' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","c","o")),
                  '1.-1'  = return(c("a","b","o","c")),
                  '-1.1'  = return(c("b","a","c","o")),
                  '-1.-1' = return(c("b","a","o","c"))
                  )
         },
         'A.4' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","o","b","o")),
                  '1.-1'  = return(c("a","o","o","b")),
                  '-1.1'  = return(c("o","a","b","o")),
                  '-1.-1' = return(c("o","a","o","b"))
                  )
         },
         'B1.5' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","a","o")),
                  '1.-1'  = return(c("a","b","o","a")),
                  '-1.1'  = return(c("b","a","a","o")),
                  '-1.-1' = return(c("b","a","o","a"))
                  )
         },
         'B2.6' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","o","a","b")),
                  '1.-1'  = return(c("a","o","b","a")),
                  '-1.1'  = return(c("o","a","a","b")),
                  '-1.-1' = return(c("o","a","b","a"))
                  )
         },
         'B3.7' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","a","b")),
                  '1.-1'  = return(c("a","b","b","a")),
                  '-1.1'  = return(c("b","a","a","b")),
                  '-1.-1' = return(c("b","a","b","a"))
                  )
         },
         'C.8' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","o","a","o")),
                  '1.-1'  = return(c("a","o","o","a")),
                  '-1.1'  = return(c("o","a","a","o")),
                  '-1.-1' = return(c("o","a","o","a"))
                  )
         },
         'D1.9' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","c","c")),
                  '1.-1'  = return(c("a","b","c","c")),
                  '-1.1'  = return(c("b","a","c","c")),
                  '-1.-1' = return(c("b","a","c","c"))
                  )
         },
         'D1.10' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","a","a")),
                  '1.-1'  = return(c("a","b","a","a")),
                  '-1.1'  = return(c("b","a","a","a")),
                  '-1.-1' = return(c("b","a","a","a"))
                  )
         },
         'D1.11' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","o","o")),
                  '1.-1'  = return(c("a","b","o","o")),
                  '-1.1'  = return(c("b","a","o","o")),
                  '-1.-1' = return(c("b","a","o","o"))
                  )
         },
         'D1.12' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("b","o","a","a")),
                  '1.-1'  = return(c("b","o","a","a")),
                  '-1.1'  = return(c("o","b","a","a")),
                  '-1.-1' = return(c("o","b","a","a"))
                  )
         },
         'D1.13' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","o","o","o")),
                  '1.-1'  = return(c("a","o","o","o")),
                  '-1.1'  = return(c("o","a","o","o")),
                  '-1.-1' = return(c("o","a","o","o"))
                  )
         },
         'D2.14' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("c","c","a","b")),
                  '1.-1'  = return(c("c","c","b","a")),
                  '-1.1'  = return(c("c","c","a","b")),
                  '-1.-1' = return(c("c","c","b","a"))
                  )
         },
         'D2.15' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","a","a","b")),
                  '1.-1'  = return(c("a","a","b","a")),
                  '-1.1'  = return(c("a","a","a","b")),
                  '-1.-1' = return(c("a","a","b","a"))
                  )
         },
         'D2.16' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("o","o","a","b")),
                  '1.-1'  = return(c("o","o","b","a")),
                  '-1.1'  = return(c("o","o","a","b")),
                  '-1.-1' = return(c("o","o","b","a"))
                  )
         },
         'D2.17' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","a","b","o")),
                  '1.-1'  = return(c("a","a","o","b")),
                  '-1.1'  = return(c("a","a","b","o")),
                  '-1.-1' = return(c("a","a","o","b"))
                  )
         },
         'D2.18' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("o","o","a","o")),
                  '1.-1'  = return(c("o","o","o","a")),
                  '-1.1'  = return(c("o","o","a","o")),
                  '-1.-1' = return(c("o","o","o","a"))
                  )
         }
         )
}
