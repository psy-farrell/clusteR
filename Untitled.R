TCscoreTrial <- function(dat, ll, crit){
  # make sure repetitions have been excluded from the data
  prlag <- dat$outpos + (ll-dat$serpos)
  prim <- sum((prlag <= crit) & recalled==1)
  rec <- sum((prlag>crit) & recalled==1)
}

scoreRepTrial <- function(dat, ll){
  # dat is from a single trial
  dat <- dat[order(dat$outpos),]
  
  isRes <- rep(FALSE, length(dat$outpos))
  
  if (ll>1){
    for (i in 2:ll){
      isRep[i] <- any(dat$serpos[i]==dat$serpos[1:(i-1)])
    }
  }
  
  dat$isRep <- isRep
  return(dat)
}

getGolomb <- function(dat, ll){
  
  tdat <- ddply(dat, .(trial), function(x){
    tx <- x[order(x$outpos),]
    outSeq <- tx$serpos[(tx$serpos>0) & (tx$recalled==1)]
    n <- length(outSeq)
    if (n>1){
      pcor <- c(1,outSeq[2:n] > outSeq[1:(n-1)])
    } else {
      pcor <- 1
    }
    tx$pcor <- rep(0, dim(tx)[1])
    # print(sum(tx$serpos>0))
    # print(length(pcor))
    tx$pcor[tx$serpos>0 & (tx$recalled==1)] <- pcor
    return(tx)
  })
  
  summ <- ddply(tdat, .(serpos), summarise, mean=mean(pcor))
  
  return(summ[summ$serpos<=ll])
}