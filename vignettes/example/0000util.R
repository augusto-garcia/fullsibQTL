
######################################################################
#
# create.map
#
# create a new map with inserted inter-marker locations
#
# Note: map is a vector or a matrix with 2 rows
# 
# stepwidth = "fixed" is what R/qtl uses; stepwidth="variable" is for
#     Brian Yandell and the bmqtl package
######################################################################
create.map <-
    function(map, step, off.end, stepwidth = c("fixed", "variable"))
{
    stepwidth <- match.arg(stepwidth)
    if(step<0 || off.end<0) stop("step and off.end must be > 0.")
    
    if(!is.matrix(map)) { # sex-ave map
        if(stepwidth == "variable") {
            if(off.end > 0) {
                tmp <- names(map)
                ## Append and prepend by off.end value (exact here, no roundoff).
                map <- c(map[1] - off.end, map, map[length(map)] + off.end)
                names(map) <- c("loc000", tmp, "loc999")
            }
            if(step == 0)
                return(unclass(map))
        
            ## Determine differences and expansion vector.
            dif <- diff(map)
            expand <- pmax(1, floor(dif / step))
            
            ## Create pseudomarker map.
            a <- min(map) + cumsum(c(0, rep(dif / expand, expand)))

            ## Names are marker names or locNNN.
            namesa <- paste("loc", seq(length(a)), sep = "")
            namesa[cumsum(c(1, expand))] <- names(map)
            names(a) <- namesa
            
            return(unclass(a))
        }
        if(length(map) == 1) { # just one marker!
            if(off.end==0) {
                if(step == 0) step <- 1
                nam <- names(map)
                map <- c(map,map+step)
                names(map) <- c(nam,paste("loc",step,sep=""))
            }
            else {
                if(step==0) m <- c(-off.end,off.end)
                else m <- seq(-off.end,off.end,by=step)
                m <- m[m!=0]
                names(m) <- paste("loc",m,sep="")
                map <- sort(c(m+map,map))
            }
            return(map)
        }
        
        minloc <- min(map)
        map <- map-minloc
        
        if(step==0 && off.end==0) return(map+minloc)
        else if(step==0 && off.end > 0) {
            a <- c(floor(min(map)-off.end),ceiling(max(map)+off.end))
            names(a) <- paste("loc", a, sep="")
            return(sort(c(a,map))+minloc)
        }
        ##############################################################################
        else if(step>0 && off.end == 0) {
            print("oi")

            print(map)
            print(min(map))
            print(floor(min(map)))
            a <- seq(floor(min(map)),max(map),
                     by = step)
            if(any(is.na(match(a, map)))) {
                print(a)
                a <- a[is.na(match(a,map))]
                print(a)
                names(a) <- paste("loc",a,sep="")
                print(a) 
                return(sort(c(a,map))+minloc)
            }
            else return(map+minloc)
        }
        ##############################################################################
        else {
            a <- seq(floor(min(map)-off.end),ceiling(max(map)+off.end+step),
                     by = step)
            a <- a[is.na(match(a,map))]
            
                                        # no more than one point above max(map)+off.end
            z <- (seq(along=a))[a >= max(map)+off.end]
            if(length(z) > 1) a <- a[-z[-1]]
            
            names(a) <- paste("loc",a,sep="")
            return(sort(c(a,map))+minloc)
        }
    } # end sex-ave map
    else { # sex-specific map
        if(stepwidth == "variable") {
            if(off.end > 0) {
                tmp <- dimnames(map)[[2]]
                map <- cbind(map[, 1] - off.end, map, map[, length(map)] + off.end)
                dimnames(map) <- list(NULL, c("loc000", tmp, "loc999"))
            }
            if(step == 0)
                return(unclass(map))
            
            ## Determine differences and expansion vector.
            dif <- diff(map[1, ])
            expand <- pmax(1, floor(dif / step))
            
            ## Create pseudomarker map.
            a <- min(map[1, ]) + cumsum(c(0, rep(dif / expand, expand)))
            b <- min(map[2, ]) + cumsum(c(0, rep(diff(map[2, ]) / expand, expand)))
            
            namesa <- paste("loc", seq(length(a)), sep = "")
            namesa[cumsum(c(1, expand))] <- dimnames(map)[[2]]
            map <- rbind(a,b)
            dimnames(map) <- list(NULL, namesa)
            
            return(unclass(map))
      }
        minloc <- c(min(map[1,]),min(map[2,]))
        map <- map-minloc
        markernames <- colnames(map)
        
        if(step==0 && off.end==0) return(map+minloc)
        else if(step==0 && off.end > 0) {
            if(ncol(map)==1) { # only one marker; assume equal recomb in sexes
                L1 <- L2 <- 1
            }
            else {
                L1 <- diff(range(map[1,]))
                L2 <- diff(range(map[2,]))
            }
            
            a <- c(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end))
            names(a) <- paste("loc", a, sep="")
            b <- c(floor(min(map[2,])-off.end)*L2/L1,
                   ceiling(max(map[2,])+off.end)*L2/L1)
            n <- c(names(a)[1],markernames,names(a)[2])
            map <- cbind(c(a[1],b[1]),map,c(a[2],b[2]))
            dimnames(map) <- list(NULL,n)
            return(map+minloc)
        }
        else if(step>0 && off.end == 0) {
            if(ncol(map)==1) return(map+minloc)
            
            a <- seq(floor(min(map[1,])),max(map[1,]),
                     by = step)
            a <- a[is.na(match(a,map[1,]))]
            
            b <- sapply(a,function(x,y,z) {
                ZZ <- min((seq(along=y))[y > x])
                (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1] }, map[1,],map[2,])
            m1 <- c(a,map[1,])
            m2 <- c(b,map[2,])
            names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
            return(rbind(sort(m1),sort(m2))+minloc)
        }
        else {
            a <- seq(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end+step),
                     by = step)
            a <- a[is.na(match(a,map[1,]))]
                                        # no more than one point above max(map)+off.end
            z <- (seq(along=a))[a >= max(map[1,])+off.end]
            if(length(z) > 1) a <- a[-z[-1]]
            
            b <- sapply(a,function(x,y,z,ml) {
              if(x < min(y)) {
                  return(min(z) - (min(y)-x)/diff(range(y))*diff(range(z)) - ml)
              }
              else if(x > max(y)) {
                  return(max(z) + (x - max(y))/diff(range(y))*diff(range(z)) - ml)
              }
              else {
                  ZZ <- min((seq(along=y))[y > x])
                  (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1]
              }
          }, map[1,],map[2,], minloc[2])
            m1 <- c(a,map[1,])
            m2 <- c(b,map[2,])
            names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
            return(rbind(sort(m1),sort(m2))+minloc)
        }
    }
}
