library(foreach)

message("Compiling C++ code...")
Rcpp::sourceCpp("code/conditional_stereotype.cpp")

logPriorLambdaS <- function(lambda,logSDiff)
{
  sum(dnorm(lambda, 0, 1, log = TRUE)) +
    sum(dnorm(logSDiff, log(0.1), 1, log = TRUE))
}


mcmcOrdinalStereotype <- function(d, 
                                  lambda0, sDiff0, 
                                  nIter=10000, thin=1, sdProp=0.1,
                                  printProgress = TRUE)
{
  logSDiff0 <- log(sDiff0)

  iIncidentStart    <- which(!duplicated(d$id))-1 # -1 for C++
  nIncidentOfficers <- as.integer(table(d$id))
  nAccept <- 0
  maxVal <- 0
  
  s0 <- c(0,1,1+cumsum(exp(logSDiff0)))
  prior0 <- logPriorLambdaS(lambda0, logSDiff0)
  logCL0 <- logCLfull(lambda0, s0,
                      d$idOff,
                      d$y,
                      iIncidentStart,    # incident start index
                      nIncidentOfficers) # incident num officers
  logA <- 0.0  

  startTime <- Sys.time()
  results <- foreach(i=1:nIter, .combine = rbind) %do%
    {
      if(is.nan(logA))
      {
        return(NULL)
      }
      
      if(printProgress && 
         (i<=10 ||
          (i<=1000 && (i %% 100)==0) ||
          (i %% 2000 == 0)))
      {
        currentTime <- Sys.time()
        timeRemaining <- difftime(currentTime, startTime, units="hours") * 
          (nIter/i - 1)
        cat("Iteration: ",i,
            " Acceptance rate: ", round(nAccept/i,2), 
            "Time remaining: ", round(timeRemaining, 1), " hours ",
            "\n")
      }
      
      # symmetric proposal
      lambda1 <- rnorm(length(lambda0), lambda0, sdProp)
      logSDiff1 <- rnorm(length(logSDiff0), logSDiff0, sdProp)
      
      s1 <- c(0,1,1+cumsum(exp(logSDiff1)))
      
      prior1 <- logPriorLambdaS(lambda1, logSDiff1)
      logCL1 <- logCLfull(lambda1, s1,
                          d$idOff,
                          d$y,
                          iIncidentStart,    # incident start index
                          nIncidentOfficers) # incident num officers

      logA <- logCL1 - logCL0 + prior1 - prior0
      
      if(is.nan(logA)) 
      { # Should never end up here
        cat("logA=NaN i=",i,"\n")
        logA <- -Inf
        save(lambda00, lambda0, lambda1, 
             logDiff00, logSDiff0, logSDiff1,
             file="../crashdump.RData")
        break
      }
      A <- min(1, exp(logA))
      
      if(runif(1) < A)
      {
        # from 2 iterations ago
        lambda00   <- lambda0
        logSDiff00 <- logSDiff0
        
        # last iteration
        lambda0   <- lambda1
        logSDiff0 <- logSDiff1
        
        prior0    <- prior1
        logCL0    <- logCL1
        nAccept   <- nAccept + 1
      }
      
      if(i %% thin == 0)
      {
        value <- c(lambda=lambda0, logSDiff=logSDiff0)
      } else
      {
        value <- NULL
      }
      
      return(value)
    }
  
  return(list(draws      = results,
              rateAccept = nAccept/nIter))
}

