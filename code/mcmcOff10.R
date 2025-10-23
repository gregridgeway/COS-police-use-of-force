
# 10 officers in a chain, 250 incidents per officer ---------------------------

library(dplyr)

RcppParallel::setThreadOptions(numThreads = 8)
source("code/mcmcConditionalStereotype.R")


set.seed(20010618)
nIncidents <- 2000
# sample first officer
i <- sample(0:9, size=nIncidents, replace=TRUE)
a <- lapply(i, function(x) 1+(x+0:2) %% 10) # mod 10 to loop back to 1/2
nOff <- sample(2:3, nIncidents, replace = TRUE)
a <- lapply(1:nIncidents, function(i) a[[i]][1:nOff[i]])
d <- data.frame(id=rep(1:nIncidents, nOff),
                idOff=unlist(a))

# break connections between 9,10 & 1,2
a <- with(d, intersect(id[idOff %in% 1:2],
                       id[idOff %in% 9:10]))
d <- d |> filter(!(id %in% a))

lambda <- c(-0.25,1.75,-1.75,1.25,2.25,-0.75,0.25,-2.25,-1.25,0.75)
s <- c(0, 1, 1.5, 2)
# h <- 0 # try with equal baseline rate of force types

p <- t(exp(s %*% t(lambda)))
p <- p/rowSums(p)

d$p <- p[d$idOff,]
d$y <- apply(d$p, 1, function(x) sample(1:4, 1, prob=x))

# check empirical probabilities resemble p
p
aggregate(cbind(y==1,y==2,y==3,y==4)~idOff, data=d, mean)

# keep only incidents that have more than one type of force used
#   such incidents contribute 1 to likelihood
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

# 100 uses-of-force per officer
idMax <- d |>
  count(id) |>
  mutate(nUoF = cumsum(n)) |>
  filter(nUoF < 10 * 250) |> # less than here
  slice_max(id) |>
  pull(id) + 1              # then +1 here
d <- d |> filter(id <= idMax)

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
    sdProp  = 0.12)
})
res$rateAccept

thetaInit <- res$draws |> tail(1) |> as.numeric()
thetaInit[(nTotOfficers+1):(nTotOfficers+2)] <- 
  exp(thetaInit[(nTotOfficers+1):(nTotOfficers+2)])

save(res, thetaInit, d, lambda, s, 
     file="output/mcmcSampOff10.RData", compress=TRUE)



# 24 uses-of-force per officer, matching SPD data -----------------------------
set.seed(20010618)

idMax <- d |> 
  count(id) |> 
  mutate(nUoF = cumsum(n)) |> 
  filter(nUoF < 10 * 24) |>
  slice_max(id) |>
  pull(id) + 1

d0 <- d |> filter(id <= idMax)

nTotOfficers <- n_distinct(d$idOff)

#thetaInit <- c(rep(0, nTotOfficers), 0.5, 0.5)
system.time({
  res <- mcmcOrdinalStereotype(
    d0, 
    lambda0 = thetaInit[1:nTotOfficers],
    sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
    nIter   = 200000,
    thin    = 20,
    sdProp  = 0.29)
})
res$rateAccept

thetaInit <- res$draws |> tail(1) |> as.numeric()
thetaInit[(nTotOfficers+1):(nTotOfficers+2)] <- 
  exp(thetaInit[(nTotOfficers+1):(nTotOfficers+2)])

save(res, thetaInit, d0, lambda, s, 
     file="output/mcmcSampOff10small.RData", compress=TRUE)

