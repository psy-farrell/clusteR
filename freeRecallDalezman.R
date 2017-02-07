source("model.R")

G$retention = "immed"
G$task = "free"

G$ll <- 15
ll <- G$ll

G$nPrevLists = 5;
G$prevLL = rep(15,G$nPrevLists)

G$recTime <- 90

G$openSet = TRUE


P$vSize <- 30;
P$intruderAct <- .002

P$lastProb <- 0.9
P$onlyLastProb <-  0.0
P$firstCue <- 1

P$iRT <- 1.5
P$totOmmT <- 7

P$dursScale <- 0.3

#G$x = c(.005, .3, .5, .003,
#        .45, 6)

P$sigma_v =  .005#  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .5 # iPhi
P$theta = .003 # T0
P$etaO = .4 # eta^O outIntG: output interference
P$T_G = 6 # giveUpG: omission criterion to give up on group

G$durs <- function(ll){ c(rep(1,ll), 1)} # last element is endDur

allSPC <- matrix(0,ll,4)
allFRP <- matrix(0,ll,4)

G$aboveTheLine <- FALSE

P$Gfun <- function(x){ sample.int(6, x, TRUE, c(5,4,6,5,2,1))}
P$Gendfun <- function(x){ sample.int(6, 1, TRUE, c(5,4,6,5,2,1))}


for (tcount in 1:4){ # start, middle, end, free
  
  G$dalCue <- tcount
  
  if (tcount<4){
    G$Dalezman <- TRUE
  } else {
    G$Dalezman <- FALSE
  }
  
  res <- model(G,P)
  
  allSPC[,tcount] <- getAccFree(res[res$serpos<16,])$pcor

  allFRP[,tcount] <- getFRP(res, 15)$prob
  
}

matplot(allSPC, type="b", ylim=c(0,1))
#matplot(allFRP, type="b")