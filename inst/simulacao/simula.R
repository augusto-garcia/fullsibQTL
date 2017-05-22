require(fullsibQTL)
set.map.fun <- "kosambi"
fs.data <- read.outcross.pheno(file="QTLfullsib.txt")
##qchisq(0.05/(65*64/2),1,lower.tail=F)/4.61
twopts <- rf.2pts(fs.data, max.rf=.4, LOD=3.9)
all <- make.seq(twopts, "all")
LGs <- group(all)
##
LG1 <- make.seq(LGs,1)
LG1.end <- make.seq(order.seq(LG1,n.init=6), "force")
LG1.end$seq.rf <- rep(0.09868766, 14)
##
LG2 <- make.seq(LGs,2)
LG2.end <- make.seq(order.seq(LG2,n.init=6), "force")
LG2.end$seq.rf <- rep(0.09868766, 14)
##
LG3 <- make.seq(LGs,3)
LG3.end <- make.seq(order.seq(LG3,n.init=6), "force")
LG3.end$seq.rf <- rep(0.09868766, 14)
##
LG4 <- make.seq(LGs,4)
LG4.end <- make.seq(order.seq(LG4,n.init=6), "force")
LG4.end$seq.rf <- rep(0.09868766, 14)

#mf.k <- function(d) 0.5*tanh(d/50)
#mf.k(10)
save("QTLexample.RData", list(fs.data, LG1.end, LG2.end, LG3.end, LG4.end))

fsib <- create.fullsib(fs.data, list(LG1.end, LG2.end, LG4.end, LG4.end))
