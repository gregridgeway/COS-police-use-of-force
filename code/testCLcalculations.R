library(dplyr)
library(ggplot2)
library(tidyr)

library(PoissonMultinomial)
library(arrangements)

Rcpp::sourceCpp("code/conditional_stereotype.cpp")

# verification of log denominator calculations --------------------------------
#   the log sums in the 2nd and 3rd lines of equation (14)

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

# this is k defined in equation (2)
k <- factor(y,levels=0:3) |> table() |> as.numeric()

# Method 1: complete enumeration of script K
arrangements::permutations(length(s), freq=k) |>
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

# Method 5: logDenomDP
#    2026-03-17 Updated with dynamic program
#    C++ calc in conditional_stereotype.cpp using dynamic program
logDenomDP(lambda, s, y)


# efficiency test --------------------------------------------------------------
#   FFT in dpmd computes entire pdf
#   DFT in logDenomDFT evaluates pdf at a single point
system.time({for(i in 1:500) log(dpmd(pmat=P, xmat=k)) + sum(log(K)) })
system.time({for(i in 1:500) logDenomDFT(lambda, s, y) })
# DFT much faster than FFT (100-200x faster)
system.time({for(i in 1:50000) logDenomDFT(lambda, s, y) })
system.time({for(i in 1:50000) logDenomDP(lambda, s, y) })
# DP faster than DFT (60-80x faster)


# verify Heap's, DFT, and DP give same results ---------------------------------
for(i in 1:1000)
{
  n <- sample(2:15, size=1)
  lambda <- runif(n,-2,2)
  s <- c(0, 1, 1.5, 2)
  y <- sample(0:3, size=length(lambda), replace=TRUE)
  
  dHeaps <- denomHeaps_R(y, lambda, s)
  dDFT   <- logDenomDFT(lambda, s, y)
  dDP    <- logDenomDP(lambda, s, y)
  
  if(abs(dHeaps-dDFT) > 0.000001 &
     abs(dHeaps-dDP) > 0.000001)
  { # should never get to this point
    list(n, lambda, s, y, dHeaps, dDFT, dDP) |> print()
    stop("Getting different answers")
  }
}

# check numerical accuracy of all three LogCL versions
m0 <- 5
lambda0 <- runif(m0, -2, 2)
s0 <- c(0, 1, 1.5, 2)
y <- rep(0:3, length.out=m0)
logCL_R(lambda0, length(lambda0), s0, y, iMethod = 1) # Heap
logCL_R(lambda0, length(lambda0), s0, y, iMethod = 2) # DFT
logCL_R(lambda0, length(lambda0), s0, y, iMethod = 3) # DP

# larger incident
m0 <- 20
lambda0 <- runif(m0, -2, 2)
s0 <- c(0, 1, 1.5, 2)
y <- rep(0:3, length.out=m0)
logCL_R(lambda0, length(lambda0), s0, y, iMethod = 2) # DFT
logCL_R(lambda0, length(lambda0), s0, y, iMethod = 3) # DP


# Heap's/DFT/DP log CL timing test -----------------------------------------------
#   ~7 minute runtime
results <- data.frame(m=2:20, heaps=NA, dft=NA, dp=NA)
for(iM in 1:nrow(results))
{
  m0 <- results$m[iM]
  message("number of officers at incident (m) = ", m0)

  lambda0 <- runif(m0, -2, 2)
  s0 <- c(0, 1, 1.5, 2)
  y <- rep(0:3, length.out=m0)
  
  if(m0 < 10) # Heaps too slow at this stage
  {
    tHeaps <- system.time({
      for(i in 1:30000)
      {
        logCL_R(lambda0, length(lambda0), s0, y, iMethod = 1)
      }
    })
  } else
  {
    tHeaps <- NA
  }
  
  tDFT <- system.time({
    for(i in 1:30000)
    {
      logCL_R(lambda0, length(lambda0), s0, y, iMethod = 2)
    }
  })
  tDP <- system.time({
    for(i in 1:30000)
    {
      logCL_R(lambda0, length(lambda0), s0, y, iMethod = 3)
    }
  })
  results$heaps[iM] <- tHeaps["elapsed"]
  results$dft[iM]   <- tDFT["elapsed"]
  results$dp[iM]    <- tDP["elapsed"]
  
  print(results[iM,])
}
# total minutes to run timing test
(results |> select(heaps, dft, dp) |> sum(na.rm=TRUE) ) /60

results |>
  pivot_longer(cols = c(heaps, dft, dp), names_to="method", values_to="seconds") |>
  mutate(method = recode(method, heaps="Heaps", dft="DFT", dp="DP")) |>
  ggplot(aes(x = m, y = seconds, colour = method)) +
  geom_line(linewidth = 1.2) +
  scale_y_log10(breaks = scales::log_breaks(n = 6, base = 10),
                labels = scales::comma_format()) +
  scale_colour_manual(values = c("Heaps" = "black", 
                                 "DFT" = "red",
                                 "DP" = "blue")) +
  labs(x = "Number of Officers",
       y = "Seconds for 30,000 log-CL calculations (log scale)",
       colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        legend.key.width = unit(1.2, "cm"))

# check rough time complexity
lm(log(heaps)~log(m), data=results |> filter(m>2))
lm(log(dft)~log(m), data=results |> filter(m>2))
lm(log(dp)~log(m), data=results |> filter(m>2))

# check DP polynomial scaling
resDP <- data.frame(m=rep(2:20, each=10), ProdK1=NA, dp=NA)
s0 <- c(0, 1, 1.5, 2)
for(iM in 1:nrow(resDP))
{
  m0 <- resDP$m[iM]
  lambda0 <- runif(m0, -2, 2)
  P <- t(exp(s0 %*% t(lambda0)))
  P <- P/rowSums(P)
  y <- apply(P, 1, function(x) sample(0:3, 1, prob=x))
  
  resDP$ProdK1[iM] <- prod(table(y)+1)
  
  nRep <- floor(20000000/m0)
  resDP$dp[iM] <- system.time({
    for(iRep in 1:nRep)
    {
      logCL_R(lambda0, length(lambda0), s0, y, iMethod = 3)
    }
    })["elapsed"]/nRep
}

lm(log(dp)~log(m), data=resDP)
lm1 <- lm(log(dp)~log(m) + log(ProdK1), data=resDP)
lm1
lm2 <- lm(log(dp)~log(m) + log(ProdK1), 
   data=filter(resDP, ProdK1>=100))
lm2

library(viridisLite)
m_vals <- sort(unique(resDP$m))
cols <- viridis(length(m_vals), option = "D")  # "D" is standard viridis
col_map <- cols[match(resDP$m, m_vals)]
plot(I(1e6*dp) ~ ProdK1, data = resDP, log="xy",
     ylab = "Seconds per 1,000,000 incidents",
     xlab = expression(prod (k[j] + 1)),
     col = col_map,
     pch = 16)
a <- coef(lm1)[1]
b <- coef(lm1)[3]
xseq <- seq(min(resDP$ProdK1), max(resDP$ProdK1), length.out = 200)
lines(xseq, 1000000*exp(a) * xseq^b, col = "black", lwd = 2)
a <- coef(lm2)[1]
b <- coef(lm2)[3]
xseq <- seq(100, max(resDP$ProdK1), length.out = 200)
lines(xseq, 1000000*exp(a) * xseq^b, col = "blue", lwd = 3)

legend("topleft", legend = sort(unique(resDP$m)),
       col = cols, pch = 16, title = "m", cex = 0.7)



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
  logCL_R(lambda1[d0$idOff+1], nrow(d0), s1, d0$y, 3) -
  logCL_R(lambda0[d0$idOff+1], nrow(d0), s0, d0$y, 3)
}) |> sum()

# in parallel, compute using logCL(1) - logCL(0)
logCLfull(lambda1, s1, d$idOff, d$y, 
          which(!duplicated(d$id))-1, 
          as.integer(table(d$id))) -
  logCLfull(lambda0, s0, d$idOff, d$y, 
            which(!duplicated(d$id))-1, 
            as.integer(table(d$id)))
