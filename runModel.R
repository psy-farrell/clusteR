rm(list = ls())
graphics.off()

# ---- libraries ---- #
library(foreach)
library(plyr)
library(dplyr)
library(farrellLab)
library(farrellMem)

# ---- intialisation ---- #
G <- list()
P <- list()

G$seed <- 15437

G$nruns <- 1000
G$init <- TRUE

P$gNoise = 0.02 # 0.02 #sigma_L in paper
P$wNoise = 0.02 #.02 #%sigma_GP in paper

P$pWeight = 0.3 # weighting of positional (vs group) information
# labelled rho in the paper
P$gWeight = 1-P$pWeight

P$iPrimacy <- 0.03 # reduction in encoding across items within groups
                    # equals 1 - gamma
P$gamma <- 1-P$iPrimacy

G$vSize = 0 # default value for N_EL

P$outIntE = .9 # eta^O_end in paper (output interference for end group if recalled first)               
P$firstCue = .0 # alpha in paper (set to a default value just in case)
P$onlyLastProb = .0 # probability of cueing with last item and then giving up
                    # set to non-zero values in specific simulation scripts

G$looseSerial = FALSE # do P's record omissions (=1) or are omissions silent (=0)

P$etaNC = .15 # eta^NC

P$intruderAct = 0 # initialising value for a_EL; is set to different values elsewhere

P$totOmmT = 5 #threshold for total omissions--set to default value just in case

# Stuff about previous experience--default values
#   Other values are set in files called below
G$openSet = TRUE # are prev lists diff items (1), or permutations of items from current list (0)?
P$phi_l = .35 # phi_g--can be modified in specific simulation scripts below

# timing parameters
P$initRT = .5
P$gRT = .25
P$iRT = .4

G$Dalezman = FALSE # Are we doing Dalezman cueing? Set to 0 by default; set to 1 in Dalezman.m
G$aboveTheLine = FALSE # Dalezman's above the line scoring?

## Here are all the different simulations
# Comment out the one you want to run

# The parameters used to call model.m are:
  # [rnoise gPhi iPhi T0 
     # outintG giveUpG]
#
# In the paper they are labelled:
  # [sigma_v phi_g phi_p theta
  # eta^O T_G]
#
# In cases where T_G doesn't apply, we effect this by setting it
# to a very large value (1000).
# Similarly, where there is no threshold on item recall, we set it to
# a very negative value.

# ---- serial-ll6 ---- #
#source("serialRecallLL6.R")

# ---- free-recall-mo70 ---- #
#source("freeRecallMO70.R")

# ---- serial-recall-grouping
#source("serialRecallGrouping.R")

# ---- serial-recall-wrapping
#source("serialRecallWrapping.R")

# ---- Ward-10
#source("Ward10.R")

# ---- free-recall-grouping
#source("freeRecallGrouping.R")

# ---- free-recall-dalezman
# source("freeRecallDalezman.R")

# ---- free-versus-serial
source("freeVsSerial.R")
