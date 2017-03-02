# TODO:

# 2. Model starts too often at start of group in other condition. Starting at the start of group is fine (PFR generally looks OK), but data imply that ppl then e.g. use next serial position rather than re-setting. Effect could also be explained by assuming that ppl try and recall multiple items from a group, or prob mix of recombination and standard free recall? <-- THIS, ONE OF YOUR IDEAS WAS TO SEE IF PPL SIMPLY DO FREE RECALL


# ---- functions yo ----
uniqGroups <- function(indat){
  trials <- unique(indat$trial)
  numer <- 0
  
  for (trial in trials){
    inseq <- indat$serpos[indat$trial==trial]
    gSeq <- indat$gPos[indat$trial==trial]
    uu <- length(unique(gSeq))
    numer <- numer + uu
  }
  return (numer/length(trials))
}

# ---- main-running ----
source("model.R")

G$retention = "immed"
G$task = "free"

G$ll <- 12
ll <- G$ll

G$nPrevLists = 5;
G$prevLL = rep(12,G$nPrevLists)

G$recTime <- 15

G$openSet = TRUE

P$vSize <- 30;
P$intruderAct <- .002

P$lastProb <- 0.9
P$onlyLastProb <- 0.0
P$firstCue <- 1

P$iRT <- 1.5

P$dursScale <- 0.3

P$recombine_freeprob <- 0.0 # probability that we try free recall in recombine condition
P$gMix <- 1

P$iPrimacy <- 0.03 # reduction in encoding across items within groups
# equals 1 - gamma
P$gamma <- 1-P$iPrimacy

P$etaNC = .15

#G$x = c(.005, .3, .5, .003,
#        .45, 6)

P$sigma_v = 0.005 # .005#  rNoise 
P$phi_g = .3 # gPhi
P$phi_p = .5 # iPhi
P$theta = .005 # theta (T0)
P$etaO = .0 # eta^O outIntG: output interference
P$T_G = 6 # giveUpG: omission criterion to give up on group

G$durs <- function(ll){ 
  x <- c(rep(1.25,ll), 1.25+.75) # last element is endDur
  x[c(1,5,9)] <- 1.25+1.5
  return(x)
} 

# P$Gfun <- function(x){ 
#   P$gCurrent <- runif(1)<P$gMix
#   if (P$gCurrent){
#     return(rep(4,x))
#   } else {
#     return(sample.int(6, x, TRUE, c(5,4,6,5,2,1)))
#   }
# }
# 
# P$Gendfun <- function(x){
#   if (P$gCurrent){
#     return(4)
#   } else {
#     return(sample.int(6, 1, TRUE, c(5,4,6,5,2,1)))
#   }
# }

allSPC <- {}
allFRP <- {}
allCRP <- {}
allres <- {}

#allTargs <- list(1:3, c(12,13,23), 123)
allTargs <- list(1:3, c(12,13,23), 123)

taskRange = c(1,3)
#taskRange = 3
#taskRange = 1

filePrefix = "noOI"

for (task in taskRange) { # loops across 1 vs 2 vs 3 groups to recall from
  
  if (task==1){
    G$Dalezman <- TRUE
    G$aboveTheLine <- FALSE
    G$recombine <- FALSE
    P$totOmmT <- 3
  } else if (task==2) {
    stop("Shouldn't get here")
  } else if (task==3) {
    G$Dalezman <- FALSE
    G$recombine <- TRUE
    G$aboveTheLine <- FALSE
    P$totOmmT <- 9
  }
  
  tcount <- 1

  
  for (targG in allTargs[[task]]){
    
    G$probeSeq <- targG
    
    if (task==1){
      G$dalCue <- G$probeSeq
    }
    
    res <- {}
    
    tnruns <- 2000
    
    if (P$gMix>0){
      P$Gfun <- function(x){ rep(4,x)}
      P$Gendfun <- function(x){ 4}
      
      G$nruns <- tnruns*P$gMix
      res <- rbind(res,model(G,P))
    }
    if (P$gMix<1){
      G$nruns <- tnruns*(1-P$gMix)
      P$Gfun <- function(x){ sample.int(6, x, TRUE, c(5,4,6,5,2,1))}
      P$Gendfun <- function(x){ sample.int(6, 1, TRUE, c(5,4,6,5,2,1))}
      
      res2 <- model(G,P)
      res2$trial <- res2$trial + tnruns
      res <- rbind(res,res2)
    }
    
    res$recalled[res$outpos>4] <- 0 # ppl could only recall 4 items
    res$cueType <- task
    res$cueType <- factor(res$cueType)
    res$cueCode <- factor(targG)
    
    if (!is.null(allres)){
      res$trial <- res$trial + dim(allres)[1]
    }
    allres <- rbind(allres,res)
    
    # SPC
    tdf <- getAccFree(res[res$serpos<13,])
    tdf$cueType <- task
    tdf$cueCode <- targG
    tdf$gPos <- rep(1:3, each=4)
    allSPC <- rbind(allSPC,tdf)
    
    # FRP
    tdf <- getFRP(res, 12)
    tdf$cueType <- task
    tdf$cueCode <- targG
    tdf$gPos <- rep(1:3, each=4)
    allFRP <- rbind(allFRP,tdf)
    
    # lag-CRP
    print(task)
    tdf <- getlagCRP(res[res$recalled!=0,], ll)
    tdf$cueType <- task
    tdf$cueCode <- targG
    allCRP <- rbind(allCRP,tdf)
    
    tcount <- tcount + 1
  }
  
}

allSPC$serposf <- factor(allSPC$serpos)
allSPC$cueType <- factor(allSPC$cueType)
allSPC$cueCode <- factor(allSPC$cueCode)
allSPC$gPos <- factor(allSPC$gPos)

allFRP$serposf <- factor(allFRP$serpos)
allFRP$cueType <- factor(allFRP$cueType)
allFRP$cueCode <- factor(allFRP$cueCode)
allFRP$gPos <- factor(allFRP$gPos)


allCRP$cueType <- factor(allCRP$cueType)
allCRP$cueCode <- factor(allCRP$cueCode)
allCRP <- allCRP[(allCRP$lag>-6) & (allCRP$lag<6),]
allCRP$lag <- factor(allCRP$lag)

allCRP <- ddply(allCRP, .(lag, cueType), summarise, lagrec=mean(lagrec, na.rm=TRUE))

pdf(paste("model_acc_",filePrefix,P$recombine_freeprob,".pdf", sep=""), width=8, height=4)
tp <- ggplot(allSPC, aes(x=as.numeric(serposf), y=pcor, colour=gPos, shape=cueCode)) + 
  geom_line() + geom_point(size=4) + facet_grid(. ~ cueType) + 
  ylim(c(0,1)) + xlab("Serial Position") + ylab("Proportion Recalled") +
  scale_shape_APA1() + theme_APA() + scale_colour_manual(values=rep("#000000",3)) + 
  scale_x_continuous(breaks=c(1,4,5,8,9,12)) + 
  guides(shape=guide_legend(title="Groups cued")) + guides(colour=FALSE)
print(tp)
dev.off()

pdf(paste("model_frp_",filePrefix,P$recombine_freeprob,".pdf", sep=""), width=8, height=4)
tp <- ggplot(allFRP, aes(x=as.numeric(serposf), y=prob, colour=gPos, shape=cueCode)) + 
  geom_line() + geom_point(size=4) + facet_grid(. ~ cueType) + 
  ylim(c(0,1)) + xlab("Serial Position") + ylab("First recall probability") +
  scale_shape_APA1() + theme_APA() + scale_colour_manual(values=rep("#000000",3)) + 
  scale_x_continuous(breaks=c(1,4,5,8,9,12)) + 
  guides(shape=guide_legend(title="Groups cued")) + guides(colour=FALSE)
print(tp)
dev.off()

pdf(paste("model_crp_",filePrefix,P$recombine_freeprob,".pdf", sep=""), width=8, height=4)
tp <- ggplot(allCRP, aes(x=lag, y=lagrec, group=cueType)) + 
  geom_line() + geom_point(size=4) + facet_grid(. ~ cueType) + 
  ylim(c(0,1)) + xlab("Lag") + ylab("Proportion Responses") +
  theme_APA() + scale_shape_APA1() + scale_colour_manual(values=rep("#000000",3))
print(tp)
dev.off()

# how many groups recalled?
allres$gPos <- rep(2,dim(allres)[1])
allres$gPos[allres$serpos<5] <- 1
allres$gPos[allres$serpos>8] <- 3

groupAcc <- ddply(allres[allres$serpos<=12 & allres$recalled==1,], .(cueType,trial), 
                  .fun=function(dat) uniqGroups(dat))
dnTrials <- ddply(allres[allres$serpos==1,], .(cueType), 
                  .fun=function(dat) dim(dat[1]))
groupAcc$V1 <- factor(groupAcc$V1)
groupAccTab <- ddply(groupAcc, .(cueType), 
                     function(x) data.frame(table(x$V1)))
groupAccTab <- ddply(groupAccTab, .(cueType), function(x){
  x$prop <- x$Freq/sum(x$Freq)
  return(x)
}
)

wasCued <- ((allres$cueType==1) & ((allres$cueCode==1 & allres$serpos<5) |
                                      (allres$cueCode==3 & allres$serpos>8) |
                                      (allres$cueCode==2 & allres$serpos>4 & allres$serpos<9))) |
  (allres$cueType==3)

tdat <- allres
tdat$wasCued <- wasCued

atLeastOne <- ddply(tdat, .(trial,cueCode,cueType), summarise, atleast1=any(wasCued & (recalled==1)))
atLeastOne_ss <- ddply(atLeastOne, .(cueCode,cueType), summarise, mean=mean(atleast1))
#atLeastOneSingle <- mean(atLeastOne_ss$mean[atLeastOne_ss$cueType==1])

pThree_mod <- prod(atLeastOne_ss$mean[atLeastOne_ss$cueType=="1"])

dodge <- position_dodge(width=0.9)
pdf(paste("model_uniqg_",filePrefix,P$recombine_freeprob,".pdf", sep=""), width=8, height=4)
tp <- ggplot(groupAccTab, aes(x=cueType, y=prop, fill=Var1)) + 
  geom_bar(position=dodge,stat="identity") +  
  geom_segment(aes(x=2.2, xend=2.4,y=pThree_mod,yend=pThree_mod))+
  ylim(c(0,1)) + xlab("Number of groups cued") + ylab("Proportion of Trials") + labs(fill="Unique Groups Accessed") +
  theme_APA() + scale_shape_APA1() + scale_fill_grey() + #scale_fill_CB() +
  theme(legend.position=c(0.6,0.8))
print(tp)
dev.off()
