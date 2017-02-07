source("model.R")

G$retention = "immed"
G$task = "serialUG"

G$ll <- 6

G$nPrevLists = 5;
G$prevLL = rep(6,G$nPrevLists)

G$recTime <- 10000

G$openSet = TRUE

G$durs <- function(ll){ c(rep(.5,ll), .5)} # last element is endDur

P$Gfun <- function(x){ 6}
P$Gendfun <- function(x){ 6}

P$vSize <- 6;
P$intruderAct <- .126

P$lastProb <- 0

P$totOmmT <- 5

P$dursScale <- 0.5

P$sigma_v = .005 #  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .5 # iPhi
P$theta = -1000 # T0
P$etaO = .6 # eta^O outIntG: output interference
P$T_G = 1000 # giveUpG: omission criterion to give up on group

res <- model(G,P)
spcAcc <- ddply(res, .(outpos), summarise, pcor = mean((serpos==outpos) & (recalled==1)))
plot(spcAcc$outpos, spcAcc$pcor, ylim=c(0,1))

par(new=T)

spcInt <- ddply(res, .(outpos), summarise, pint = mean(serpos>6))
plot(spcInt$outpos, spcInt$pint, ylim=c(0,1))

transmat <- getTransMat(res, 6)
matplot(transmat, type="l")
