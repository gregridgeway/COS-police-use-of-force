# Placebo test ----------------------------------------------------------------
#   Appendix F
#   Run the conditional ordinal stereotype model on permuted versions 
#   of the SPD data. Permutes y within incident

library(dplyr)

RcppParallel::setThreadOptions(numThreads = 8)

# load C++ code for conditional likelihood calcs and acceptance probability
source("code/mcmcConditionalStereotype.R")

# prep SPD data
load("data/dataSPD.RData")

# make 0-based for C++
d$y     <- d$y - 1
d$idOff <- d$idOff - 1

nTotOfficers <- n_distinct(d$idOff)

for(iPermute in 1:5)
{
  message("Permutation: ",iPermute,"\n")
  set.seed(20010618 + iPermute)
  
  # permute y within incident
  d <- d |>
    group_by(id) |>
    mutate(y=sample(y)) |>
    ungroup()
  
  thetaInit <- c(rep(0.0, nTotOfficers), 0.5, 0.5)
  
  # originally split across two different machines
  #   with different number of iterations to make sure
  #   1M was enough for convergence
  nIters <- ifelse(iPermute %in% c(1,4), 2000000, 1000000)
  
  resSPD <- mcmcOrdinalStereotype(
    d, 
    lambda0 = thetaInit[1:nTotOfficers],
    sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
    nIter   = nIters,
    thin    = 100, 
    sdProp  = 0.025) 

  thetaInit <- resSPD$draws |> tail(1) |> as.numeric()
  thetaInit[(nTotOfficers+1):(nTotOfficers+2)] <- 
    exp(thetaInit[(nTotOfficers+1):(nTotOfficers+2)])

  save(resSPD, thetaInit, 
       file=paste0("output/mcmcSampSPDPermutechain",iPermute,".RData"),
       compress = TRUE)
}
