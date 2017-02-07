source("model.R")

G$retention = "immed"

G$ll <- 10
ll <- G$ll

G$nPrevLists = 5;
G$prevLL = rep(10,G$nPrevLists)

G$recTime <- 100

G$openSet = TRUE

G$looseSerial <- TRUE

P$vSize <- 30;
P$intruderAct <- .002

P$iRT <- 1.5
P$totOmmT <- 10

P$dursScale <- 0.3

P$sigma_v =  .005#  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .6 # iPhi
P$theta = .005 # T0
P$etaO = .45 # eta^O outIntG: output interference
P$T_G = 6 # giveUpG: omission criterion to give up on group

G$durs <- function(ll){ c(rep(0.8,ll), 0.8)} # last element is endDur

P$Gfun <- function(x){ sample.int(5, x, TRUE, c(6,4,2,1,1))}
P$Gendfun <- function(x){ sample.int(5, 1, TRUE, c(6,4,2,1,1))}

allSPC <- matrix(0,ll,4)
allFRP <- matrix(0,ll,4)

for (tcount in 1:4){
  
  if (tcount<3){
    G$task="free"
    P$lastProb = .7; 
    P$firstCue = 0.25; 
    P$onlyLastProb = .3; 
  } else {
    P$lastProb = 0;
    P$firstCue = 0;
    P$onlyLastProb = .3;
    G$task="serialG"
  }  
  
  if ((tcount%%2)>0){
    P$phi_l <- .35
  } else {
    P$phi_l <- .5
  }
  
  res <- model(G,P)
  
  if (tcount<3){
    allSPC[,tcount] <- getAccFree(res[res$serpos<11,])$pcor
  } else {
    allSPC[,tcount] <- getGolomb(res, 10)$pcor
  }

  allFRP[,tcount] <- getFRP(res, 10)$prob
}

par(mfrow=c(1,2))
matplot(allSPC, type="b", ylim=c(0,1))
matplot(allFRP, type="b")

