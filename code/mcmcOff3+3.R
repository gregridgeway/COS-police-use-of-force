
library(dplyr)

# load C++ code for conditional likelihood calcs and acceptance probability
source("mcmcConditionalStereotype.R")


################################################################################
# simulate with two unconnected groups 3+3 officers
################################################################################
# Seattle example nInc = 4821, nOff = 1503, 
# 0-28807, 1-6243, 2-1885, 3-81

set.seed(20010618)
nIncidents <- 500
nOff <- sample(2:3, nIncidents, replace = TRUE)
d1 <- data.frame(id=rep(1:nIncidents, nOff),
                 idOff=sapply(nOff, function(x) sample(1:3, x)) |> unlist())
nOff <- sample(2:3, nIncidents, replace = TRUE)
d2 <- data.frame(id=rep(1:nIncidents, nOff),
                 idOff=sapply(nOff, function(x) sample(1:3, x)) |> unlist())
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
  mutate(id = as.integer(as.factor(id))) # renumber 1, 2, 3, ...

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

save(res, thetaInit, d, lambda, s, file="mcmcSampOff6.RData", compress=TRUE)


plot(res$draws[,nTotOfficers+1])
plot(res$draws[,nTotOfficers+2])
plot(res$draws[,1])

c(1 + mean(exp(res$draws[,nTotOfficers+1])),
  1 + mean(exp(res$draws[,nTotOfficers+1])) +
    mean(exp(res$draws[,nTotOfficers+2])))


################################################################################
# simulate with two unconnected groups 3+3 officers
#   with 10 incidents per officer
################################################################################

set.seed(20010618)
# use 30 with 1,2,3 and 30 with 4,5,6
i <- d |> 
  group_by(idOff <= 2) |> 
  summarize(minID = min(id),
            maxID = max(id)) 
i <- c(sample(i$minID[2]:i$maxID[2], size=30),
       sample(i$minID[1]:i$maxID[1], size=30))
d0 <- d |> 
  filter(id %in% i) |>
  mutate(id = as.integer(as.factor(id))) # renumber 1, 2, 3, ...


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

save(res, thetaInit, d0, lambda, s, file="mcmcSampOff6small.RData", compress=TRUE)


plot(res$draws[,nTotOfficers+1])
plot(res$draws[,nTotOfficers+2])
plot(res$draws[,1])

c(1 + mean(exp(res$draws[,nTotOfficers+1])),
  1 + mean(exp(res$draws[,nTotOfficers+1])) +
    mean(exp(res$draws[,nTotOfficers+2])))
