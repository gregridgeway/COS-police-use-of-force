
library(dplyr)

RcppParallel::setThreadOptions(numThreads = 8)
source("code/mcmcConditionalStereotype.R")

# 3+3 disconnected, 250 incidents per officers --------------------------------
set.seed(20010618)
nIncidents <- 2000
nOff <- sample(2:3, nIncidents, replace = TRUE)
d1 <- data.frame(id    = rep(1:nIncidents, nOff),
                 idOff = sapply(nOff, function(x) sample(1:3, x)) |> unlist(),
                 group = 1)
nOff <- sample(2:3, nIncidents, replace = TRUE)
d2 <- data.frame(id    = rep(1:nIncidents, nOff),
                 idOff = sapply(nOff, function(x) sample(1:3, x)) |> unlist(),
                 group = 2)
d2$id    <- d2$id    + max(d1$id) 
d2$idOff <- d2$idOff + max(d1$idOff) 
d <- rbind(d1,d2)

lambda <- -c(0.5, 1, 2, 0, 0.75, 1.5)
s <- c(0, 1, 1.5, 2)

p <- t(exp(s %*% t(lambda)))
p <- p/rowSums(p)

d$p <- p[d$idOff,]
d$y <- apply(d$p, 1, function(x) sample(1:4, 1, prob=x))

# check empirical probabilities resemble p
p
aggregate(cbind(y==1,y==2,y==3,y==4)~idOff, data=d, mean)

# keep only incidents that have more than one type of force used
d <- d |>
  right_join(d |> 
              group_by(id) |>
              summarize(nUniqueY = n_distinct(y)) |>
              filter(nUniqueY > 1) |>
              select(id),
             by = join_by(id)) |>
  select(-p) |>
  mutate(id = as.integer(as.factor(id))) |> # renumber 1, 2, 3, ...
  arrange(id, idOff)

# 250 uses-of-force per officer
idMax <- d |>
  group_by(group) |>
  count(id) |>
  mutate(nUoF = cumsum(n)) |>
  filter(nUoF < 3 * 250) |> # less than here
  slice_max(id) |>
  pull(id) + 1              # then +1 here

d <- d |> 
  filter((group==1 & id <= idMax[1]) |
         (group==2 & id <= idMax[2]))

# make 0-based for C++
d$y     <- d$y - 1
d$idOff <- d$idOff - 1


nTotOfficers <- n_distinct(d$idOff)

thetaInit <- c(rep(0, nTotOfficers), 0.5, 0.5)
system.time({
  res <- mcmcOrdinalStereotype(
    d, 
    lambda0 = thetaInit[1:nTotOfficers],
    sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
    nIter   = 200000,
    thin    = 20,
    sdProp  = 0.13)
})
res$rateAccept

thetaInit <- res$draws |> tail(1) |> as.numeric()
thetaInit[(nTotOfficers+1):(nTotOfficers+2)] <- 
  exp(thetaInit[(nTotOfficers+1):(nTotOfficers+2)])

save(res, thetaInit, d, lambda, s, 
     file="output/mcmcSampOff6.RData", compress=TRUE)



# 3+3 disconnected, 24 UoFs per officer --------------------------------

set.seed(20010618)

# 24 uses-of-force per officer, matching SPD data
idMax <- d |> 
  group_by(group) |>
  count(id) |> 
  mutate(nUoF = cumsum(n)) |> 
  filter(nUoF < 3 * 24) |>
  slice_max(id) |>
  pull(id) + 1

d0 <- d |> 
  filter((group==1 & id <= idMax[1]) |
         (group==2 & id <= idMax[2]))


system.time({
  res <- mcmcOrdinalStereotype(
    d0, 
    lambda0 = thetaInit[1:nTotOfficers],
    sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
    nIter   = 200000,
    thin    = 20, 
    sdProp  = 0.45)
})
res$rateAccept

thetaInit <- res$draws |> tail(1) |> as.numeric()
thetaInit[(nTotOfficers+1):(nTotOfficers+2)] <- 
  exp(thetaInit[(nTotOfficers+1):(nTotOfficers+2)])

save(res, thetaInit, d0, lambda, s, 
     file="output/mcmcSampOff6small.RData", compress=TRUE)

