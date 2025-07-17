# Main analysis of SPD use-of-force data --------------------------------------

library(dplyr)

RcppParallel::setThreadOptions(numThreads = 8)
source("code/mcmcConditionalStereotype.R")

# load SPD data
load("data/dataSPD.RData")
# make 0-based for C++
d$y     <- d$y - 1
d$idOff <- d$idOff - 1
nTotOfficers <- n_distinct(d$idOff)

for(iChain in 1:4)
{
  message("Chain: ",iChain,"\n")
  set.seed(20010618 + iChain)
  thetaInit <- c(runif(nTotOfficers,-2,2), 0.5,0.5)
  if(iChain==4)
    thetaInit <- c(rep(0, nTotOfficers), 0.5, 0.5)
  resSPD <- mcmcOrdinalStereotype(
    d, 
    lambda0 = thetaInit[1:nTotOfficers],
    sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
    nIter   = 1000000,
    thin    = 100, 
    sdProp  = 0.03)  # sdProp tuned to accept 23-25%
  resSPD$rateAccept

  thetaInit <- resSPD$draws |> tail(1) |> as.numeric()
  thetaInit[(nTotOfficers+1):(nTotOfficers+2)] <- 
    exp(thetaInit[(nTotOfficers+1):(nTotOfficers+2)])

  save(resSPD, thetaInit, 
       file=paste0("output/mcmcSampSPDchain",iChain,".RData"),
       compress = TRUE)
}
