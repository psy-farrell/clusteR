model <- function(G,P,durs=NULL){
  
  #### Initialisation ####
  if (G$init){
    set.seed(1441561)
  }
  
  # How many items do we probe for?
  if (is.null(G$nProbe)){
    nProbe <- G$ll
  } else {
    nProbe <- G$nProbe
  }
  
  #allProbe = matrix(0,G$nruns,nProbe) # records order in which items were actually probed
  
  if (G$retention=='immed'){
    listSource <- c("curr",
                    rep("prev", length(G$prevLL)))
    llStruct = c(G$ll, G$prevLL)
  } else if (G$retention=='delay'){
    listSource <- c("curr",
                    rep("prev", length(G$prevLL)),
                    "dist")
    llStruct = c(G$ll, G$prevLL, G$lldist)
  } #TODO: CD task
  totL <- length(llStruct)
  
  nItems <- ifelse(G$openSet,sum(llStruct),llStruct[1])
  
  # set this up to do parallel-wise
  
  #### Start of trial loop ####
  res <- foreach (trial=1:G$nruns, .combine='rbind') %do% {
    
  #for (trial in 1:G$nruns){
    
    groupStruct <- list()
    groupMark <- list()
    itStruct <- lapply(1:nItems, function(x) list(p=NULL, etaGV=NULL))
    cumnGroups <- 0
    nG <- 0
    nGroups <- 0
    
    eShift <- 0 # used to shift item indexes to refer to whole experiment (not just current list)
    
    # set up lists
    for (listi in 1:totL){
      
      tLL <- llStruct[listi]
      
      # not set up to do delayed SR
      if (G$task=='serialUG'){
        tGroupStruct = tLL
        tnGroups = 1
        tGroupSize = 1
        endSize <- tLL
      } else {
        
        if (listSource[listi]=="curr"){
          if (G$retention=='CD'){
            tGroupStruct = ones(1,tLL) # we assume that each item is in its own group
            endSize= 1;
          } else {
            tGroupStruct = P$Gfun(tLL)
            endSize = P$Gendfun(tLL)
          }
        } else if (listSource[listi]=="prev"){
          tGroupStruct = P$Gfun(tLL)
          endSize = P$Gendfun(tLL)
        } else if (listSource[listi]=="dist"){
          tGroupStruct <- tLL # a single group containing the distractors
          endSize <- tLL
        }
      }
      
      # at this point, tGroupStruct often contains too many groups,
      # as to be safe we make as many groups as there are list items
      # (just in case all groups happen by chance to contain a single
      # item--unlikely, but it could happen!)
      
      # the next bit finds the last group, sets it to size endSize,
      # and does some truncation if necessary 
      # to make the number of groups and group sizes
      # pack in to the list length
      
      tGroupInfo <- truncateGroups(tGroupStruct,endSize, tLL)
      cumGroupStruct <- cumsum(tGroupInfo$tGroupStruct)
      itStart <- 1
      
      if (listSource[listi]=="curr"){
        itind <- 1:tLL
        itStart <- tLL+1
      } else if (listSource[listi]=="prev"){
        if (G$openSet){
          itind <- itStart:(itStart+tLL-1)
          itStart <- itStart + tLL
        } else { #closed set
          itind <- sample.int(tLL)
        }
      } else {
        itind <- NULL
      }
      # if distractor, don't actually encode item (assume orthogonal)
      
      
      gStart <- 1
      tnGroups <- length(tGroupInfo$tGroupStruct)
      
      tDurs <- G$durs(tLL)
      
      for (groupi in 1:tnGroups){
        lOcc <- ifelse(listSource[listi]=="dist",1,listi)
        gEnd <- cumGroupStruct[groupi]
        groupProp <- tGroupInfo$tGroupProp[groupi]
        groupSize <- tGroupInfo$tGroupStruct[groupi]
        
        bump <- (tDurs[gStart]+tDurs[gEnd+1])/2
        dursBump <- (1-exp(-bump*P$dursScale))
        
        etaLC <-  dursBump + rnorm(1,sd=P$gNoise) # Equation A3
        
        # add info about the group to the list
        groupStruct[[length(groupStruct)+1]] <- list(
                              etaLC=etaLC,
                              dursBump=dursBump,
                              lOcc=lOcc,
                              gStart=gStart,
                              gEnd=gEnd,
                              groupSize=groupSize,
                              type=listSource[listi],
                              groupProp = groupProp,
                              groupPos = groupi,
                              expPos = (tnGroups - groupi)+cumnGroups
                         )

        gRange <- itind[gStart:gEnd]
        
        if (listSource[listi]!="dist"){
          
          for (i in 1:groupSize){
            
            if (G$openSet){
              iInd <- gRange[i] + eShift
            } else {
              iInd <- gRange[i]
            }
            
            if (groupSize>1){
              itStruct[[iInd]]$p <- c(itStruct[[iInd]]$p,
                                      (i-1)/(groupSize-1)*groupProp)
            } else {
              itStruct[[iInd]]$p <- c(itStruct[[iInd]]$p,0) # just put in 0 for group size 1
            }
            
            itStruct[[iInd]]$etaGV <- c(itStruct[[iInd]]$etaGV,
                                             ((P$gamma^(i-1)) + rnorm(1, sd=P$wNoise)) * P$phi_l^(listi-1))
            itStruct[[iInd]]$expPos <- c(itStruct[[iInd]]$expPos,
                                              (tnGroups - groupi)+cumnGroups)
          }
        }
        
        gStart <- cumGroupStruct[groupi] + 1
        
      }
      
      cumnGroups <- cumnGroups + tnGroups
      if (listi==1){
        nG <- tnGroups
      }
      
      eShift <- eShift + tLL
    }
    
    nGroups <- length(groupStruct)
    
    if (G$recombine){
      groupsLeft <- as.numeric(strsplit(as.character(G$probeSeq), "")[[1]])
    }
    
    #### Retrieval ####
    newGroup <- TRUE # are we recalling a new group (yes we are)?
    groupCount = 0  # how many groups (including current one) recalled so far?
    recGroups = {} # tracks groups recalled
    
    r <- rep(0, nItems + P$vSize) # initialize response suppression vector 
    
    gSuppress <- rep(FALSE,nGroups) # used to track recall of groups for suppression
    
    if ((!G$openSet) && (P$vSize>0)){
      r[(G$ll+1):length(r)] <- 1 # suppress items outside current list (expt vocab)
    }
    
    gOmm <-  0 # tracks number of successive omissions within group
    totOmm <-  0 # tracks total number of omissions
    
    currContext <- 1 # current list context
    
    rt <- P$initRT # add on time to begin recall (see section "Dynamics of recall") 
    cumrt <- P$initRT
    
    litem <- 1
    allProbep <- 1 # tracks the items that were probed for; mostly relevant for wrapping expt
    givingUp <- FALSE # if set to 1 later on, we've given up and are just recording omissions in output file
    
    # Sample the recall strategies to be used on this trial
    useLast <- runif(1) < P$lastProb # start with last group
    onlyLast <- runif(1) < P$onlyLastProb # start with last and then stop (see Simulation 8)
    useStart <- runif(1) < P$firstCue # did we nominally encode first group irrespective of preceding factors
    if (G$recombine){
      doRecombine <- runif(1) < (1-P$recombine_freeprob) # in recombination, do we actually do recombination (vs free recall)
    } else {
      doRecombine <- FALSE
    }
    
    litem <- 1
    
    tResp <- {}
    tRT <- {}
    tCumRT <- {}
    tRecalled <- {}
    tCue <- {}
    
    while ((litem <= nProbe) && (litem <= G$ll)){
  
      # ------------retrieve a new group if needed
      if (newGroup){
        cueForEnd <- FALSE # indicates whether we are recalling last group first
                          #(if so, we have intact group context and don't need to
                          # retrieve)
        groupCount <- groupCount+1;
        
        nullGroup <- FALSE
        
        c_LC <- sapply(groupStruct, function(x) sum(x$etaLC*P$phi_l^(x$lOcc-currContext)))
        c_NC <- c_LC*0 # initialise, is usually modified below
        
        if (G$task=="free"){
          
          if ((groupCount==1) && (G$retention=="immed") && G$Dalezman && (G$dalCue<4)){
            if (G$dalCue==1){
              dCue <- 1
            } else if (G$dalCue==2){
              if (nG>3){ # if there isn't one identifiable "middle group", randomly select
                if (runif(1)<0){
                  dCue <- sample(2:(nG-1),1,prob = rep(1,nG-2))
                } else {
                  dCue <- ceiling(nG/2)
                }
              } else {
                dCue <- 2
              }
            } else if (G$dalCue==3){
              dCue <- nG
            }
            
            c_NC <- sapply(groupStruct, function(x){
              ifelse(x$groupPos==dCue,P$etaNC,0)
            })
            
            if (dCue==nG){
              cueForEnd <- TRUE
            }
            
          } else if (G$recombine && doRecombine && (length(groupsLeft)>0)){
            G$dalCue <- sample(groupsLeft, 1)
            
            if (G$dalCue==1){
              dCue <- 1
            } else if (G$dalCue==3){
              dCue <- nG
            } else if (G$dalCue==2){
              if (nG>3){ # if there isn't one identifiable "middle group", randomly select
                if (runif(1)<0){
                  dCue <- sample(2:(nG-1),1,prob = rep(1,nG-2))
                } else {
                  dCue <- ceiling(nG/2)
                }
              } else {
                dCue <- 2
              }
            } 
            
            c_NC <- sapply(groupStruct, function(x){
              ifelse(x$groupPos==dCue,P$etaNC,0)
            })
            if (dCue==nG && (groupCount==1)){
              cueForEnd <- TRUE
            }
          } else if ((groupCount==1) && (G$retention=="immed") && (useLast || onlyLast)){
            cueForEnd <- TRUE
          } else if (useStart && !any(recGroups==1)) {
            c_NC <- sapply(groupStruct, function(x){
              ifelse(x$groupPos==1,P$etaNC,0)
            })
            #c_LC[groupPos==1] = c_LC[groupPos==1] + P$etaNC;
          }
        } else if (G$task=="serialUG" || G$task=="serialG"){
          groupProbe <- min(groupCount, nG) #this makes us continue cueing for last group just in case
          # we run beyond end of the list
          
          if (G$task=="serialUG" || nG==1){
            cueForEnd <- TRUE
          } else {
            
            if (onlyLast){
              cueForEnd <- TRUE # this is primarily for Simulation 8
              
            # we ignore useLast here, though I guess people could use a strategy
            # of recalling the last item and then going back to start in serial recall
            
            } else { #---use cue for the group we want
              c_NC <- sapply(groupStruct, function(x){
                ifelse(x$groupPos==groupProbe,P$etaNC,0)
              })
            }
          }
        } else if (G$task=="wrap"){
          groupProbe = (G$startG+groupCount+1)%%nG+1;
          
          if (groupCount==1 && groupProbe==nG && useLast){
            cueForEnd <- TRUE
          } else {
            c_NC <- sapply(groupStruct, function(x){
              ifelse(x$groupPos==groupProbe,P$etaNC,0)
            })
          }
        }
        
        if (cueForEnd){
          # if immediately probing for last group, no need to retrieve
          pCG <- 1; # proportion of unit activations carried over to test
          currGroup <- nG;
          nullGroup <- FALSE;
          currExpPos <- groupStruct[[nG]]$expPos
        } else {
          # winner-takes-all selection of a group
          currGroup <- which.max(c_LC + c_NC)
          
          currExpPos <- groupStruct[[currGroup]]$expPos
          
          # add on time taken to retrieve group context
          rt <- rt + P$gRT
          cumrt <- cumrt + P$gRT
          
          if (any(which(gSuppress)==currGroup)){
            pCG <- 0
            nullGroup <- TRUE
          } else {
            nullGroup <- FALSE
            pCG <- groupStruct[[currGroup]]$dursBump
          }
        } # end of cueForEnd conditional
          
        if (!nullGroup){
          if (cueForEnd){
            groupStruct[[currGroup]]$etaLC <- 
              c(groupStruct[[currGroup]]$etaLC, P$etaO * P$outIntE)
          } else {
            groupStruct[[currGroup]]$etaLC <- 
              c(groupStruct[[currGroup]]$etaLC, P$etaO)
          }
          groupStruct[[currGroup]]$lOcc <- 
            c(groupStruct[[currGroup]]$lOcc,2)
        }
        
        gSuppress[currGroup] <- TRUE
        
        recGroups <- c(recGroups, currGroup)
        if (G$recombine && doRecombine){
          groupsLeft <- setdiff(groupsLeft, currGroup)
        }
        
        newGroup <- FALSE
        
        if (groupStruct[[currGroup]]$groupSize > 1){
          # gDiff <- groupStruct[[currGroup]]$groupSize-1
          # 
          # gSeq <- seq(groupStruct[[currGroup]]$gStart,groupStruct[[currGroup]]$gEnd)-
          #   groupStruct[[currGroup]]$gStart
        
          #pVec <- gSeq/gDiff*groupStruct[[currGroup]]$groupProp
          pVec <- seq(0,1,length.out = groupStruct[[currGroup]]$groupSize)*groupStruct[[currGroup]]$groupProp
        } else {
          pVec <- 0
        }
        
        inGroupCount <- 1
        
      } else { # not a new group
        inGroupCount <- inGroupCount + 1
        
      } # end if new group conditional
      
      #------------- ITEM RETRIEVAL
      if ((G$task=='serialG') || (G$task=='serialUG')){
        targItem = litem;
      } else if  (G$task=='wrap') {
        targItem = ((litem + (G$startG-1)*3 -1)%%G$ll)+1;
      } else {
        targItem = NaN; # weren't aiming for any particular item; just fill in with NaN
      }
      
      #allProbe[datai] <- targItem
      
      v_GV <- sapply(itStruct, function(x) {
          sum(x$etaGV*P$phi_g^abs(x$expPos-currExpPos))
        })
      v_GV <- v_GV/sum(v_GV)
      
      
      v_PV <- sapply(itStruct, function(x) {
        sum(x$etaGV*P$phi_p^abs(x$p-pVec[inGroupCount]))
      })
      v_PV <- v_PV/sum(v_PV)
      
      v <- P$gWeight*pCG*v_GV + P$pWeight*v_PV #Eq A15
      v <- c(v, rep(P$intruderAct,P$vSize))
      
      a <- (v + rnorm(n=length(v), sd=P$sigma_v))*(1-r)
      
      a_ord <- order(a, decreasing=TRUE)
      
      updateLitem <- 1
      
      if ( (a[a_ord[1]]-a[a_ord[2]]) > P$theta){
        itptr <- a_ord[1]
      } else if (!givingUp) {
        itptr <- -1
        if ((G$task=="free") || G$looseSerial){
          updateLitem <- 0
        }
      }
      
      # add on item RTs
      thisRT <- rexp(1, 1/P$iRT)
      rt <- rt + thisRT
      cumrt <- cumrt + thisRT
      
      if (givingUp){
        updateLitem <- 1
      } else if (cumrt > G$recTime){ # out of time
        lRemain <- length(setdiff(1:G$ll,tResp))
        tResp <- c(tResp, setdiff(1:G$ll,tResp))
        tCue <- c(tCue, setdiff(1:G$ll,tResp))
        tRecalled <- c(tRecalled, rep(0, lRemain))
        tRT <- c(tRT, rep(NA, lRemain))
        tCumRT <- c(tCumRT, rep(NA, lRemain))
        givingUp <- TRUE # don't do any more scoring
        if (G$task=='free'){
          updateLitem <- 100000 # forces exit from trial
        }
      } else if ((itptr>0) && (!givingUp)) { # record the response
      
        tCue <- c(tCue, targItem)
        r[itptr] = 1; # suppress the item
      
        if (itptr>nItems){ # extra-experiment recall
          tResp <- c(tResp, 1000 + litem) # arbitrary numerical code
          tRecalled <- c(tRecalled, -1)
        } else { # recall from list, or possibly preceding list
          # itptr <= list length is recall from current list
          tResp <- c(tResp, itptr) # arbitrary numerical code
          if (itptr>G$ll){
            tRecalled <- c(tRecalled, -1)
          } else {
            tRecalled <- c(tRecalled, 1)
          }
        }
        tRT <- c(tRT, rt)
        tCumRT <- c(tCumRT, cumrt)
      
        rt <- 0; # we made a response, so now start IRT from scratch
      } else {
        
        # failed to recall any item
        totOmm = totOmm + 1;
        gOmm = gOmm+1;
        
        # if not free recall or "loose" serial recall, make an omit response
        if (G$task!='free' && !G$looseSerial){
          tCue <- c(tCue, targItem)
          tResp <- c(tResp, -1)
          tRT <- c(tRT, rt) # in matlab code we recorded NaN; here we record RT
          tCumRT <- c(tCumRT, cumrt)
          rt <- 0
          tRecalled <- c(tRecalled, 0)
        }
        
        if (totOmm >= P$totOmmT){
          lRemain <- length(setdiff(1:G$ll,tResp))
          tResp <- c(tResp, setdiff(1:G$ll,tResp))
          tCue <- c(tCue, setdiff(1:G$ll,tResp))
          tRecalled <- c(tRecalled, rep(0, lRemain))
          tRT <- c(tRT, rep(NA, lRemain))
          tCumRT <- c(tCumRT, rep(NA, lRemain))
          givingUp <- TRUE # don't do any more scoring
          #print("decided to give up")
          if (G$task!="free"){
            updateLitem <- 1 #keep ticking through to record probed positions for remaining omissions
          } else {
            updateLitem <- 10000 #force exit from trial
          }
        } else if (gOmm > P$T_G){
          #print("to many omissions")
          newGroup = TRUE
        }
      }
      
      # Start new group if reached length of current group (assumed to be known)
      if (inGroupCount==groupStruct[[currGroup]]$groupSize){
        #print("reached end of current")
        newGroup <- TRUE
      }
      
      # Start new group if recombining and have recalled an item from this group
      if (G$recombine && doRecombine && (length(groupsLeft)>0) && ((itptr>0) && (!givingUp))){
        #print("this happened")
        newGroup <- TRUE
      }
      
      if (((newGroup) && (groupCount==1) && (G$aboveTheLine)) ||
        (newGroup && onlyLast)){
          lRemain <- length(setdiff(1:G$ll,tResp))
          tResp <- c(tResp, setdiff(1:G$ll,tResp))
          tCue <- c(tCue, setdiff(1:G$ll,tResp))
          tRecalled <- c(tRecalled, rep(0, lRemain))
          tRT <- c(tRT, rep(NA, lRemain))
          tCumRT <- c(tCumRT, rep(NA, lRemain))
          updateLitem <- 100000
      }
      
      litem <- litem + updateLitem
    
    } # end of litem loop
    
    # resp <- c(resp, tResp)
    # rt <- c(rt, tRT)
    # recalled <- c(recalled, tRecalled)
    
    lDiff <- length(tResp)-length(tCue)
    tCue <- c(tCue, rep(0,lDiff))
    
    return(data.frame(trial=rep(trial,length(tResp)),
                      outpos=1:length(tResp),
                      cue=tCue,
                      serpos=tResp,
                      recalled=tRecalled,
                      rt=tRT))
    
  } #end of trial loop/foreach
  
  return(res)
  
} # end of function


#### OTHER FUNCTIONS

truncateGroups <- function(tGroupStruct, endSize, tLL){
  cumG <- cumsum(tGroupStruct) #cumG stores terminal serial position of each group
  
  if (endSize > 0) {#--make the last group of size endSize
    lg <- which(cumG>(tLL-endSize))[1] # find the last group that will run in
    # to the end group given endSize
    
    if (lg>1){
      tGroupStruct[lg]=tLL - endSize - cumG[lg-1];
    } else {
      tGroupStruct[lg]=tLL - endSize;
    }
    
    if (endSize>tLL){
      endSize = tLL;
    }
    
    tGroupStruct[lg+1] <- endSize
    tGroupStruct <- tGroupStruct[1:(lg+1)]
    tGroupStruct <- tGroupStruct[tGroupStruct>0]
    
    tGroupSize <- rep(1, length(tGroupStruct))
    
  } else { # -- if endsize = 0, signals that we truncate the last group
    # tGroupSize scales within-group markers to reflect the
    # fact that a group may have finished before we expected it
    # to
    lg <- which(cumG>=tLL)[1]
    if (lg>1){
      tGroupStruct[lg] <- tLL -cumG[lg-1]
    } else {
      tGroupStruct[lg] <- tLL
    }
    
    tGroupStruct = tGroupStruct[1:lg]
    
    tGroupSize = rep(1, length(tGroupStruct))
    if (lg>1){
      tGroupSize[lg] = (tLL - cumG[lg-1])/(cumG[lg]-cumG[lg-1])
    } else {
      tGroupSize[lg] = tLL/cumG[1];
    }
  }
  
  return(list(tGroupStruct=tGroupStruct,
              tGroupProp=tGroupSize))
}