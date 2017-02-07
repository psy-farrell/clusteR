source("model.R")

G$retention = "immed"
G$task = "free"

G$ll <- 18
ll <- G$ll

G$nPrevLists = 5;
G$prevLL = rep(18,G$nPrevLists)

G$recTime <- 60

G$openSet = TRUE


P$vSize <- 30;
P$intruderAct <- .002

P$lastProb <- 0.9
P$onlyLastProb <-  0.0
P$firstCue <- 0.2

P$iRT <- 1.5
P$totOmmT <- 5

P$dursScale <- 0.3

#G$x = c(.005, .3, .5, .003,
#        .45, 6)

P$sigma_v =  .005#  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .5 # iPhi
P$theta = .003 # T0
P$etaO = .45 # eta^O outIntG: output interference
P$T_G = 6 # giveUpG: omission criterion to give up on group

allSPC <- matrix(0,ll,2)
allFRP <- matrix(0,ll,2)
allCondRec <- matrix(0,8,2)


for (tcount in 1:2){
  

  if (tcount==1){
    tt <- 14.5/18
    G$durs <- function(ll){ c(rep(tt,ll), tt)} # last element is endDur
    
    P$Gfun <- function(x){ sample.int(6, x, TRUE, c(5,4,6,5,2,1))}
    P$Gendfun <- function(x){ sample.int(6, 1, TRUE, c(5,4,6,5,2,1))}
    
  } else {
    
    tt = 14.5/(18+6)
    G$durs <- function(ll){ 
      x <- c(rep(tt,ll), tt*2) # last element is endDur
      x[c(4,7,10,13,16)] <- tt*2
      return(x)
    } 
    
    P$Gfun <- function(x){ rep(3,x)}
    P$Gendfun <- function(x){ 3}
    
  }
  
  res <- model(G,P)
  
  allSPC[,tcount] <- getAccFree(res[res$serpos<19,])$pcor

  allFRP[,tcount] <- getFRP(res, 18)$prob

  allCondRec[,tcount] <- getCondRec(res, 18)$condRec[1:8]
}

matplot(allSPC, type="b")
matplot(allFRP, type="b")
matplot(allCondRec, type="b")