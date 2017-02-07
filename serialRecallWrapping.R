source("model.R")

G$retention = "immed"

G$ll <- 9

G$nPrevLists = 5;
G$prevLL = rep(9,G$nPrevLists)

G$durs <- function(ll){
  td <- c(rep(.45,ll), .9)
  td[1] <- .9
  td[4] <- .9
  td[7] <- .9
  return(td)
}

P$Gfun <- function(x){ rep(3,x)}
P$Gendfun <- function(x){3}

G$recTime <- 1000

G$openSet = FALSE

P$vSize <- 30;
P$intruderAct <- .002

P$lastProb <- 1

P$totOmmT <- 5

P$dursScale <- 0.15

#G$x = c(.005, .3, .5, .003,
#        .45, 6)

P$sigma_v =  .005#  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .5 # iPhi
P$theta = .003 # T0
P$etaO = .45 # eta^O outIntG: output interference
P$T_G = 1000 # giveUpG: omission criterion to give up on group

allSPC <- {}
alltrans <- {}

G$task <- 'wrap'

for (tcount in 1:3){
  
  print(tcount)
  
  G$startG <- tcount
  
  res <- model(G,P)
  
  spcAcc <- ddply(res[res$outpos<=9,], .(outpos), summarise, pcor = mean((serpos==cue) & (recalled==1)))
  allSPC <- cbind(allSPC,spcAcc$pcor)
  
}

matplot(allSPC, type="l", ylim=c(0,1))