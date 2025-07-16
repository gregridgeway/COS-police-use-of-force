
library(dplyr)

RcppParallel::setThreadOptions(numThreads = 8)
source("code/mcmcConditionalStereotype.R")


# 3 officers, ~100 incidents per officer --------------------------------------

set.seed(20010618)
nIncidents <- 500
# number of officers per incident
nOff <- sample(2:3, nIncidents, replace = TRUE)
# randomly assign officers to incidents
d <- data.frame(id=rep(1:nIncidents, nOff),
                idOff=unlist(sapply(nOff, function(x) sample(1:3, x))))
# simulate with these parameters
lambda <- -c(0.5, 1, 2)
s <- c(0, 1, 1.5, 2)


## Figure 2 -------------------------------------------------------------------
p <- t(exp(c(-0.15, -0.75, -1.50, -3.00) + s %*% t(lambda)))
barplot(p, beside=TRUE, 
        names.arg=c("Witness", "Level 1", "Level 2", "Level 3"),
        xlab="",
        ylab="Probability of force type")
text(x=5.5, y=p[1,2] + 0.12, expression(lambda == -frac(1,2)))
text(x=6.5, y=p[1,3] + 0.12, expression(lambda == -1))
text(x=7.5, y=p[1,4] + 0.12, expression(lambda == -2))


# translate lambda & s to probabilities of using type of force
p <- t(exp(s %*% t(lambda)))
p <- p/rowSums(p)
# simulate force used
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
               select(id)) |>
  select(-p) |>
  mutate(id = as.integer(as.factor(id))) # renumber 1, 2, 3, ...


# make 0-based for C++
d$y     <- d$y - 1
d$idOff <- d$idOff - 1

lambda0 <- -c(0.5, 1, 2)
lambda1 <- lambda0 + c(0,0,1.5)
s0 <- c(0, 1, 1.5, 2)
s1 <- s0 + c(0, 0, 0.1, 0.1) 


nTotOfficers <- n_distinct(d$idOff)

thetaInit <- c(rep(0, nTotOfficers), 0.5, 0.5)
system.time({
  res <- mcmcOrdinalStereotype(
    d, 
    lambda0 = thetaInit[1:nTotOfficers],
    sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
    nIter   = 200000,
    thin    = 20, 
    sdProp  = 0.2)
})
res$rateAccept

thetaInit <- res$draws |> tail(1) |> as.numeric()
thetaInit[(nTotOfficers+1):(nTotOfficers+2)] <- 
  exp(thetaInit[(nTotOfficers+1):(nTotOfficers+2)])

save(res, thetaInit, d, lambda, s, file="mcmcSampOff3.RData", compress=TRUE)



# 3 officers, 10 incidents per officer --------------------------------------
set.seed(20010618)
d0 <- d |> filter(id <= 30)

nTotOfficers <- n_distinct(d0$idOff)

thetaInit <- c(rep(0, nTotOfficers), 0.5, 0.5)
system.time({
  res <- mcmcOrdinalStereotype(
    d0, 
    lambda0 = thetaInit[1:nTotOfficers],
    sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
    nIter   = 200000,
    thin    = 20, 
    sdProp  = 0.65)
})
res$rateAccept

thetaInit <- res$draws |> tail(1) |> as.numeric()
thetaInit[(nTotOfficers+1):(nTotOfficers+2)] <- 
  exp(thetaInit[(nTotOfficers+1):(nTotOfficers+2)])

save(res, thetaInit, d0, lambda, s, file="mcmcSampOff3small.RData", compress=TRUE)
