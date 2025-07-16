library(dplyr)
library(ggplot2)
library(tidyr)

library(PoissonMultinomial)
library(arrangements)

Rcpp::sourceCpp("code/conditional_stereotype.cpp")

# verification of log denominator calculations --------------------------------
#   the log sums in the second line of equation (14)

# generate some values
#   will compare P(y|lambda, s) with several methods to verify calculations

# individual officers' propensity to escalate
lambda <- runif(10,-2,2)
# the inverse "distances" between types of force
s <- c(0, 1, 1.5, 2)
# P = prob. that officer i uses force type j
P <- t(exp(s %*% t(lambda)))
K <- rowSums(P)
P <- P/K

# generate some force types
#   doesn't matter values, some values at which we'll compute P(y)
y <- sample(0:3, size=length(lambda), replace=TRUE)

# this is the k defined just before equation (2)
k <- factor(y,levels=0:3) |> table() |> as.numeric()


# Method 1: complete enumeration of script K
#   scriptK = all possible ways to assign force types in y to officers
scriptK <- table(factor(y, 0:(length(s)-1)))
arrangements::permutations(length(s), freq=scriptK) |>
  apply(1, function(i) sum(lambda*s[i])) |>
  exp() |>
  sum() |>
  log()

# Method 2: no-repeat Heap's algorithm
#    C++ calculation in conditional_stereotype.cpp
denomHeaps_R(y, lambda, s)

# Method 3: discrete Fourier transform from PoissonMultinomial package
log(dpmd(pmat=P, xmat=k)) + sum(log(K))

# Method 4: logDenomDFT 
#    C++ calculation in conditional_stereotype.cpp
logDenomDFT(lambda, s, y) 


# efficiency test, FFT v DFT --------------------------------------------------
#   FFT in dpmd computes entire pdf
#   DFT in logDenomDFT evaluates pdf at a single point
system.time({for(i in 1:500) log(dpmd(pmat=P, xmat=k)) + sum(log(K)) })
system.time({for(i in 1:500) logDenomDFT(lambda, s, y) })
# Conclusion: DFT much faster than FFT (100-200x faster)


# verify Heap's and DFT give same results -------------------------------------
for(i in 1:1000)
{
  n <- sample(2:8, size=1)
  lambda <- runif(n,-2,2)
  s <- c(0, 1, 1.5, 2)
  y <- sample(0:3, size=length(lambda), replace=TRUE)
  
  dHeaps <- denomHeaps_R(y,lambda,s)
  dDFT   <- logDenomDFT(lambda,s,y)

  if(abs(dHeaps-dDFT) > 0.000001)
  { # should never get to this point
    stop("Heaps and DFT give different answers")
  }
}


# Heap's/DFT log CL timing test -----------------------------------------------
#   test to determine at what size incidents to transition from Heap's to DFT
#   2.5 minute runtime
#   Heap's faster for m<=7, DFT for m>=8
results <- data.frame(m=2:10, heaps=NA, dft=NA)
for(iM in 1:nrow(results))
{
  m0 <- results$m[iM]
  message("number of officers at incident (m) = ", m0)

  lambda0 <- runif(m0, -2, 2)
  s0 <- c(0, 1, 1.5, 2)
  y <- rep(0:3, length.out=m0)
  
  tHeaps <- system.time({
    for(i in 1:30000)
    { # Heap's limit=15 forces Heap's only
      logCL_R(lambda0, length(lambda0), s0, y, iHeapsLimit = 15)
    }
  })
  tDFT <- system.time({
    for(i in 1:30000)
    {
      # Heap's limit=1 forces DFT only
      logCL_R(lambda0, length(lambda0), s0, y, iHeapsLimit = 1)
    }
  })
  results$heaps[iM] <- tHeaps["elapsed"]
  results$dft[iM]   <- tDFT["elapsed"]
}
sum(results$heaps+results$dft)/60 # total minutes to run timing test

results |>
  pivot_longer(cols = c(heaps, dft), names_to="method", values_to="seconds") |>
  mutate(method = recode(method, heaps="Heaps", dft="DFT")) |>
  ggplot(aes(x = m, y = seconds, colour = method)) +
  geom_line(linewidth = 1.2) +
  scale_y_log10(breaks = scales::log_breaks(n = 6, base = 10),
                labels = scales::comma_format()) +
  scale_colour_manual(values = c("Heaps" = "black", "DFT" = "red")) +
  labs(x = "Number of Officers",
       y = "Seconds for 30,000 log-CL calculations (log scale)",
       colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        legend.key.width = unit(1.2, "cm"))


# test acceptance probability calculations ------------------------------------
set.seed(20010618)
RcppParallel::setThreadOptions(numThreads = 16)

# simulate data from three officers
nIncidents <- 500
nOff <- sample(2:3, nIncidents, replace = TRUE)
d <- data.frame(id=rep(1:nIncidents, nOff),
                idOff=unlist(sapply(nOff, function(x) sample(1:3, x))))
lambda <- -c(0.5, 1, 2)
s <- c(0, 1, 1.5, 2)
p <- t(exp(s %*% t(lambda)))
p <- p/rowSums(p)
d$p <- p[d$idOff,]
d$y <- apply(d$p, 1, function(x) sample(1:4, 1, prob=x))

a <- aggregate(y~id, data=d, function(x) length(unique(x)))
d <- subset(d, id %in% a$id[a$y>1])
d$id <- as.numeric(as.factor(as.character(d$id)))
d <- d[order(d$id),]
d$p <- NULL

# make 0-based for C++
d$y     <- d$y - 1
d$idOff <- d$idOff - 1

# compute acceptance probability from (lambda0, s0) to (lambda1, s1)
lambda0 <- -c(0.5, 1, 2)
lambda1 <- lambda0 + c(0,0,1.5)
s0 <- c(0, 1, 1.5, 2)
s1 <- s0 + c(0, 0, 0.1, 0.1) 

# serial, compute logCL(1) - logCL(0)
a <- split(d, d$id)
sapply(a, function(d0)
{
  logCL_R(lambda1[d0$idOff+1], nrow(d0), s1, d0$y, 7) -
  logCL_R(lambda0[d0$idOff+1], nrow(d0), s0, d0$y, 7)
}) |> sum()

# in parallel, compute using logCL(1) - logCL(0)
logCLfull(lambda1, s1, d$idOff, d$y, 
          which(!duplicated(d$id))-1, 
          as.integer(table(d$id))) -
  logCLfull(lambda0, s0, d$idOff, d$y, 
            which(!duplicated(d$id))-1, 
            as.integer(table(d$id)))
