

scoreRepTrial <- function(dat, ll){
  # dat is from a single trial
  dat <- dat[order(dat$outpos),]
  
  isRep <- rep(FALSE, length(dat$outpos))
  
  if (ll>1){
    for (i in 2:ll){
      isRep[i] <- any(dat$serpos[i]==dat$serpos[1:(i-1)])
    }
  }
  
  dat$isRep <- isRep
  return(dat)
}

source("model.R")

G$retention = "immed"
G$task = "free"

G$ll <- 20

G$nPrevLists = 5;
G$prevLL = rep(8,G$nPrevLists)

G$recTime <- 60

G$openSet = TRUE

G$durs <- function(ll){ c(rep(1,ll), 1)} # last element is endDur

P$Gfun <- function(x){ sample.int(6, x, TRUE, c(5,4,6,5,2,1))}
P$Gendfun <- function(x){ 0}

P$vSize <- 30;
P$intruderAct <- .002

P$lastProb <- 0.9
P$onlyLastProb <-  0.0
P$firstCue <- 0.2

P$iRT <- 1.5
P$totOmmT <- 8

P$dursScale <- 0.5

#G$x = c(.005, .3, .5, .003,
#        .45, 6)

P$sigma_v =  .005#  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .5 # iPhi
# P$theta is specified below
P$etaO = .45 # eta^O outIntG: output interference
P$T_G = 6 # giveUpG: omission criterion to give up on group

tRange = 1:15

allSPC <- list()
allFRP <- matrix(0,15,3)
allCRP1 <- rep(NA,15)
allTC <- matrix(0,15,2)
allCondRec <- list()

for (ll in tRange){
  
  print(ll)
  G$ll <- ll
  P$theta <- 0.015/log(ll+1)
  P$lastProb <- ll/15
  
  res <- model(G,P)
  
  allSPC[[ll]] <- getAccFree(res[res$serpos<(ll+1),])$pcor
  frp <- getFRP(res, ll)$prob
  allFRP[ll,1] <- frp[1]
  if (ll>1){
    allFRP[ll,3] <- sum(frp[max(2,ll-3):ll])
  }
  allFRP[ll,2] <- 1 - (allFRP[ll,1]+allFRP[ll,3])
  
  if (ll>1){
    crp <- getlagCRP(res[res$recalled!=0,], ll)$lagrec
    allCRP1[ll] <- crp[ll+1]
  }
  
  #TC scoring
  tres <- ddply(res, .(trial), function(x) scoreRepTrial(x,ll))
  tres <- tres[!tres$isRep,]
  baba <- ddply(tres, .(trial), function(x) TCscoreTrial(x,ll = ll, crit=7))
  allTC[ll,] <- colMeans(baba)[2:3]
  
  # conditional recency
  allCondRec[[ll]] <- getCondRec(res,ll)
  
}

plot(unlist(allSPC), type="n", xlim=c(1,max(sapply(allSPC,length))))
mapply(lines,allSPC,col=seq_along(allSPC),type="b",lty=2, pch=1)
# borrowed from http://stackoverflow.com/questions/18179856/how-to-plot-a-list-of-vectors-with-different-lengths

matplot(allFRP, type="b", pch=1)

matplot(allTC, type="b", pch=1)
legend(10,4,c("prim","rec"), lty=1:2, col=1:2, pch=1)
