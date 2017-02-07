source("model.R")

G$retention = "immed"
G$task = "free"

G$ll <- 20

G$nPrevLists = 5;
G$prevLL = rep(20,G$nPrevLists)

G$recTime <- 60

G$openSet = TRUE

G$durs <- function(ll){ c(rep(.75,ll), .75)} # last element is endDur

P$Gfun <- function(x){ sample.int(6, x, TRUE, c(5,4,6,5,2,1))}
P$Gendfun <- function(x){ sample.int(6, 1, TRUE, c(5,4,6,5,2,1))}

P$vSize <- 30;
P$intruderAct <- .002

P$lastProb <- 0.9
P$onlyLastProb <-  0.0
P$firstCue <- 0.3

P$iRT <- 1.5
P$totOmmT <- 7

P$dursScale <- 0.3

#G$x = c(.005, .3, .5, .003,
#        .45, 6)

P$sigma_v =  .005#  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .5 # iPhi
P$theta = .003 # T0
P$etaO = .45 # eta^O outIntG: output interference
P$T_G = 6 # giveUpG: omission criterion to give up on group

res <- model(G,P)

plot(getAccFree(res[res$serpos<21,])$pcor)

plot(getFRP(res, 20)$prob)

plot(getlagCRP(res[res$recalled!=0,], 20)$lagrec)

plot(getCondRec(res, 20)$condRec[1:8])