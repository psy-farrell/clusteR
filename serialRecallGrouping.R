source("model.R")

G$retention = "immed"

G$ll <- 9

G$nPrevLists = 5;
G$prevLL = rep(9,G$nPrevLists)

G$recTime <- 1000

G$openSet = FALSE

P$vSize <- 0;
P$intruderAct <- .002

P$lastProb <- 0

P$totOmmT <- 5

P$dursScale <- 0.5

#G$x = c(.005, .3, .5, .003,
#        .45, 6)

P$sigma_v =  .005#  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .6 # iPhi
P$theta = .003 # T0
P$etaO = .65 # eta^O outIntG: output interference
P$T_G = 1000 # giveUpG: omission criterion to give up on group

allSPC <- {}
alltrans <- {}

for (tcount in 1:3){
  
  print(tcount)
  
  if (tcount==1){
    G$task="serialUG"
  } else {
    G$task="serialG"
  }
  
  if (tcount==1 || tcount==2){
    G$durs <- function(ll){ c(rep(.6,ll), .6)}
  } else {
    G$durs <- function(ll){
      td <- c(rep(.45,ll), .9)
      td[4] <- .9
      td[7] <- .9
      return(td)
    } # last element is endDur
  }
  
  if (tcount==2){
    P$Gfun <- function(x){ sample.int(9, x, TRUE, c(0,1,1,1,1,0,0,0,1))}
    P$Gendfun <- function(x){ sample.int(9, 1, TRUE, c(0,1,1,1,1,0,0,0,1))}
  } else {
    P$Gfun <- function(x){ rep(3,x)}
    P$Gendfun <- function(x){3}
  }
  
  res <- model(G,P)
  
  spcAcc <- ddply(res[res$outpos<=9,], .(outpos), summarise, pcor = mean((serpos==outpos) & (recalled==1)))
  allSPC <- cbind(allSPC,spcAcc$pcor)
  
  tgrad <- getTransGradient(res, 9, abstrans=TRUE)
  alltrans <- cbind(alltrans, tgrad$density[2:4]/sum(tgrad$density[2:9]))
  
  
  # if there was an interposition, what is prob that interposition is next?
  tnumer <- 0
  tdenom <- 0
  res$trans <- res$outpos - res$serpos
  for (i in 1:(length(res$trans)-1)){
    if (res$outpos[i]<9){
      if (((abs(res$trans[i])==3) || (abs(res$trans[i])==6)) &&
          (res$recalled[i]==1)){
        tdenom <- tdenom + 1
        if ((res$trans[i]==res$trans[i+1]) &&
            (res$recalled[i+1]==1) ){
          tnumer <- tnumer +1
        }
      }
    }
  }
  print(tnumer/tdenom)
  
}

matplot(allSPC, type="l", ylim=c(0,1))

matplot(1:3, alltrans, type="l", ylim=c(0,0.5))