library(dplyr)
library(igraph)
library(xtable)
library(future)
library(doFuture)
library(progressr)
library(ggplot2)
library(ggbeeswarm)

# 3 officers ------------------------------------------------------------------
load("output/mcmcSampOff3.RData")
nOff <- 3
n_distinct(d$id)

tab <- data.frame(parm=c(paste0("\\lambda_1-\\lambda_{",2:3,"}"),
                         paste0("\\mathrm{rank}(\\lambda_",1:3,")"),
                         "s_2","s_3"),
                  true=c(lambda[1]-lambda[-1],
                         rank(lambda),
                         s[3:4]),
                  postmean=NA,
                  l95=NA,
                  u95=NA,
                  postmeanSmall=NA,
                  l95small=NA,
                  u95small=NA)


## Figure 3 -------------------------------------------------------------------
net0 <- lapply(split(d$idOff, d$id),
               function(x)
               {
                 b <- expand.grid(id1=x,id2=x) |>
                   filter(id1<id2)
               })
net0 <- bind_rows(net0) |>
  group_by(id1,id2) |>
  count()
net <- graph_from_data_frame(net0, directed=FALSE)
nodeColor <- colorRamp(heat.colors(10))
col <- nodeColor((rank(lambda)-1)/2)
col <- apply(col, 1, function(x) rgb(x[1],x[2],x[3], maxColorValue = 255))
par(mai=0.02+c(0.2,0.2,0.2,0.2))
plot(net,
     layout = layout.reingold.tilford,
     edge.width=net0$n/10,
     #vertex.color=col,
     vertex.color="white",
     vertex.label.cex=1.5,
     vertex.label=paste(1:nOff,lambda,sep="\n"),
     vertex.size=30,
     asp=0.3)

postDraws <- list(lambda = res$draws[,1:nOff],
                  sDelta = exp(res$draws[,nOff + 1:2]))
postDraws$s <- 1 + apply(postDraws$sDelta, 1, cumsum) |> t()
colnames(postDraws$s) <- c("s2","s3")

## Table 1, left panel --------------------------------------------------------
a <- apply(postDraws$s, 2, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^s_", tab$parm)
tab[i,c("postmean","l95","u95")] <- a

a <- apply(postDraws$lambda[,1]-postDraws$lambda[,2:3], 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_", tab$parm)
tab[i,c("postmean","l95","u95")] <- a

a <- apply(postDraws$lambda, 1, rank) |>
  apply(1, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(1) |> t()
i <- grep("^\\\\mathrm\\{rank\\}", tab$parm)
tab[i,c("postmean","l95","u95")] <- a


# "all the values of Var(λj | λi, data) were small, less than 0.06 for the 250 uses-of-force per officer simulation..."
Sigma <- cov(postDraws$lambda)
schurComp <- sapply(1:nOff, function(i) diag(Sigma) - Sigma[,i]^2/Sigma[i,i]) |> t() |> zapsmall()
range(schurComp[schurComp>0])

# not used in article, but show pairs plot
pairs(postDraws$lambda,
      labels=c(expression(lambda[1]),
               expression(lambda[2]),
               expression(lambda[3])),
      pch=".")
plot(postDraws$s, xlab=expression(s[2]), ylab=expression(s[3]), 
     pch=".", log="xy")


# 3 officers, 10 incidents/officer --------------------------------------------
load("output/mcmcSampOff3small.RData")
nOff <- 3
n_distinct(d0$id)

postDraws <- list(lambda = res$draws[,1:nOff],
                  sDelta = exp(res$draws[,nOff + 1:2]))
postDraws$s <- 1 + apply(postDraws$sDelta, 1, cumsum) |> t()
colnames(postDraws$s) <- c("s2","s3")

## Table 1, right panel --------------------------------------------------------
a <- apply(postDraws$s, 2, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^s_", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a

a <- apply(postDraws$lambda[,1]-postDraws$lambda[,2:3], 2, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a

a <- apply(postDraws$lambda, 1, rank) |>
  apply(1, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(1) |> t()
i <- grep("^\\\\mathrm\\{rank\\}", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a

# "...and less than 0.28 for the 24 uses-of-force per officer simulation."
Sigma <- cov(postDraws$lambda)
schurComp <- sapply(1:nOff, function(i) diag(Sigma) - Sigma[,i]^2/Sigma[i,i]) |> 
  t() |> zapsmall()
range(schurComp[schurComp>0])


## Appendix C, Figure C1 ------------------------------------------------------
pdf("output/pairs3off1.pdf", width=7, height=5)
pairs(postDraws$lambda, 
      labels=c(expression(lambda[1]),
               expression(lambda[2]),
               expression(lambda[3])),
      cex.axis = 0.8,
      pch=".")
dev.off()

## Appendix C, Figure C2 ------------------------------------------------------
pdf("output/s2s3.pdf", width=7, height=5)
plot(postDraws$s, xlab=expression(s[2]), ylab=expression(s[3]), 
     pch=".", log="xy")
dev.off()

## Appendix C, Figure C3 ------------------------------------------------------
pdf("output/sDiffs.pdf", width=7, height=5)
plot(postDraws$s[,1]-1, (postDraws$s[,2]-postDraws$s[,1]), 
     xlab=expression(s[2]-s[1]), 
     ylab=expression(s[3]-s[2]), 
     pch=".", 
     log="xy", 
     cex.axis=0.8)
dev.off()

## Table 1 ---------------------------------------------------------------------
tab |>
  mutate(interval = paste0("(",trimws(format(l95,nsmall=2)),", ",
                               trimws(format(u95,nsmall=2)),")"),
         intervalSmall = paste0("(",trimws(format(l95small,nsmall=2)),", ",
                                    trimws(format(u95small,nsmall=2)),")"),
         parm = paste0("$",parm,"$")) |>
  relocate(interval, .before=4) |>
  relocate(intervalSmall, .before=8) |>
  select(-l95, -u95, -l95small, -u95small) |>
  xtable() |>
  print(sanitize.text.function = identity,
        include.rownames = FALSE)




# 3+3 officers ----------------------------------------------------------------
load("output/mcmcSampOff6.RData")
nOff <- n_distinct(d$idOff)
n_distinct(d$id)

d$idOff <- d$idOff + 1 # back to R indexing

tab <- data.frame(parm=c(paste0("\\lambda_1-\\lambda_{",2:6,"}"),
                         paste0("\\lambda_4-\\lambda_{",5:6,"}"),
                         "s_2","s_3"),
                  true=c(lambda[1]-lambda[-1],
                         lambda[4]-lambda[5:6],
                         s[3:4]),
                  postmean=NA,
                  l95=NA,
                  u95=NA,
                  postmeanSmall=NA,
                  l95small=NA,
                  u95small=NA)

tabRanks <- data.frame(
  parm=c(paste0("\\mathrm{rank}(\\lambda_",1:3,"|\\mathcal{O}_",1:3,")"),
         paste0("\\mathrm{rank}(\\lambda_",4:6,"|\\mathcal{O}_",4:6,")")),                  true=c(rank(lambda[1:3]), rank(lambda[4:6])),
  postmean=NA,
  l95=NA,
  u95=NA,
  postmeanSmall=NA,
  l95small=NA,
  u95small=NA)

## Figure 1 -------------------------------------------------------------------
net0 <- lapply(split(d$idOff, d$id),
               function(x)
               {
                 b <- expand.grid(id1=x,id2=x) |>
                   filter(id1<id2)
               })
net0 <- bind_rows(net0) |>
  group_by(id1,id2) |>
  count()
net <- graph_from_data_frame(net0, directed=FALSE)
i <- match(names(V(net)), 1:nOff)
nodeColor <- colorRamp(heat.colors(10))
col <- nodeColor((rank(lambda[i])-1)/(nOff-1))
col <- apply(col, 1, function(x) rgb(x[1],x[2],x[3], maxColorValue = 255))
par(mai=0.02+c(0.2,0.2,0.2,0.2))
plot(net,
     layout = layout.reingold.tilford,
     edge.width=net0$n/10,
     #vertex.color=col,
     vertex.color="white",
     vertex.label.cex=1.5,
     vertex.label=paste(i,lambda[i],sep="\n"),
     vertex.size=30,
     asp=0.3)


postDraws <- list(lambda = res$draws[,1:nOff],
                  sDelta = exp(res$draws[,nOff + 1:2]))
postDraws$s <- 1 + apply(postDraws$sDelta, 1, cumsum) |> t()
colnames(postDraws$s) <- c("s2","s3")

## Table 2, left panel --------------------------------------------------------
a <- apply(postDraws$s, 2, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975))))  |>
  round(2) |> t()
i <- grep("^s_", tab$parm)
tab[i,c("postmean","l95","u95")] <- a

a <- apply(postDraws$lambda[,1]-postDraws$lambda[,2:6], 2, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_1", tab$parm)
tab[i,c("postmean","l95","u95")] <- a

a <- apply(postDraws$lambda[,4]-postDraws$lambda[,5:6], 2, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_4", tab$parm)
tab[i,c("postmean","l95","u95")] <- a



## Table D1, left panel --------------------------------------------------------
a <- apply(postDraws$lambda[,1:3], 1, rank) |>
  apply(1, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(1) |> t()
i <- grep("lambda_[123]", tabRanks$parm)
tabRanks[i,c("postmean","l95","u95")] <- a
a <- apply(postDraws$lambda[,4:6], 1, rank) |>
  apply(1, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(1) |> t()
i <- grep("lambda_[456]", tabRanks$parm)
tabRanks[i,c("postmean","l95","u95")] <- a

# "With the larger simulated dataset, Var(λ−1 | λ1, data) = [0.00 0.01 0.03 0.34 0.34 0.34], signaling which officers are Officer 1's informational peer."
Sigma <- cov(postDraws$lambda)
schurComp <- sapply(1:nOff, function(i) diag(Sigma) - Sigma[,i]^2/Sigma[i,i]) |> t() |> zapsmall()
schurComp[1,] |> round(2)

# show pairs plot, not used in article
postDraws$lambda |> 
  data.frame() |> 
  pairs(labels=c(expression(lambda[1]), expression(lambda[2]),
                 expression(lambda[3]), expression(lambda[4]),
                 expression(lambda[5]), expression(lambda[6])),
        pch=".")


# 3+3 officers, 24 incidents/officer ------------------------------------------
load("output/mcmcSampOff6small.RData")
nOff <- n_distinct(d$idOff)
d0$idOff <- d0$idOff + 1 # back to R indexing

postDraws <- list(lambda = res$draws[,1:nOff],
                  sDelta = exp(res$draws[,nOff + 1:2]))
postDraws$s <- 1 + apply(postDraws$sDelta, 1, cumsum) |> t()
colnames(postDraws$s) <- c("s2","s3")

## Table 2, right panel -------------------------------------------------------
a <- apply(postDraws$s, 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^s_", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a

a <- apply(postDraws$lambda[,1]-postDraws$lambda[,2:6], 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_1", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a

a <- apply(postDraws$lambda[,4]-postDraws$lambda[,5:6], 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_4", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a


## Table D1, right panel -------------------------------------------------------
a <- apply(postDraws$lambda[,1:3], 1, rank) |>
  apply(1, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(1) |> t()
i <- grep("lambda_[123]", tabRanks$parm)
tabRanks[i,c("postmeanSmall","l95small","u95small")] <- a
a <- apply(postDraws$lambda[,4:6], 1, rank) |>
  apply(1, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(1) |> t()
i <- grep("lambda_[456]", tabRanks$parm)
tabRanks[i,c("postmeanSmall","l95small","u95small")] <- a

# not shown in article
Sigma <- cov(postDraws$lambda)
schurComp <- sapply(1:nOff, function(i) diag(Sigma) - Sigma[,i]^2/Sigma[i,i]) |> 
  t() |> zapsmall()
schurComp[1,] |> round(2)


## Table 2 --------------------------------------------------------------------
tab |>
  mutate(interval = paste0("(",trimws(format(l95,nsmall=2)),", ",
                           trimws(format(u95,nsmall=2)),")"),
         intervalSmall = paste0("(",trimws(format(l95small,nsmall=2)),", ",
                                trimws(format(u95small,nsmall=2)),")"),
         parm = paste0("$",parm,"$")) |>
  relocate(interval, .before=4) |>
  relocate(intervalSmall, .before=8) |>
  select(-l95, -u95, -l95small, -u95small) |>
  xtable() |>
  print(sanitize.text.function = identity,
        include.rownames = FALSE)


## Table D1 --------------------------------------------------------------------
tabRanks |>
  mutate(interval = paste0("(",trimws(format(l95,nsmall=1)),", ",
                           trimws(format(u95,nsmall=1)),")"),
         intervalSmall = paste0("(",trimws(format(l95small,nsmall=1)),", ",
                                trimws(format(u95small,nsmall=1)),")"),
         parm = paste0("$",parm,"$")) |>
  relocate(interval, .before=4) |>
  relocate(intervalSmall, .before=8) |>
  select(-l95, -u95, -l95small, -u95small) |>
  xtable() |>
  print(sanitize.text.function = identity,
        include.rownames = FALSE)



# 10 officers in a chained network --------------------------------------------
load("output/mcmcSampOff10.RData")
nOff <- n_distinct(d$idOff)
n_distinct(d$id)
d$idOff <- d$idOff + 1 # back to R indexing

tab <- data.frame(parm=c(paste0("\\lambda_1-\\lambda_{",2:10,"}"),
                         paste0("\\lambda_5-\\lambda_{",(1:10)[-5],"}"),
                         "s_2","s_3"),
                  true=c(lambda[1]-lambda[-1],
                         lambda[5]-lambda[-5],
                         s[3:4]),
                  postmean=NA,
                  l95=NA,
                  u95=NA,
                  postSD=NA,
                  postmeanSmall=NA,
                  l95small=NA,
                  u95small=NA,
                  postSDsmall=NA)

tabRanks <- data.frame(parm=paste0("\\mathrm{rank}(\\lambda_",1:10,")"),
                       true=rank(lambda),
                       postmean=NA,
                       l95=NA,
                       u95=NA,
                       postmeanSmall=NA,
                       l95small=NA,
                       u95small=NA)


## Figure 4 -------------------------------------------------------------------
net0 <- lapply(split(d$idOff, d$id),
               function(x)
               {
                 b <- expand.grid(id1=x,id2=x) |>
                   filter(id1<id2)
               })
net0 <- bind_rows(net0) |>
  group_by(id1,id2) |>
  count()
net <- graph_from_data_frame(net0, directed=FALSE)
i <- match(names(V(net)), 1:nOff)
nodeColor <- colorRamp(heat.colors(10))
col <- nodeColor((rank(lambda[i])-1)/(nOff-1))
col <- apply(col, 1, function(x) rgb(x[1],x[2],x[3], maxColorValue = 255))
par(mai=0.02+c(0.2,0.2,0.2,0.2))
plot(net,
     layout = layout.reingold.tilford,
     edge.width=net0$n/10,
     vertex.color="white",
     vertex.label.cex=1.5,
     vertex.label=paste(i,lambda[i],sep="\n"),
     vertex.size=20,
     asp=0.45)


postDraws <- list(lambda = res$draws[,1:nOff],
                  sDelta = exp(res$draws[,nOff + 1:2]))
postDraws$s <- 1 + apply(postDraws$sDelta, 1, cumsum) |> t()
colnames(postDraws$s) <- c("s2","s3")

## Table 3 & Table D2, left panel ----------------------------------------------
a <- apply(postDraws$s, 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("s_", tab$parm)
tab[i,c("postmean","l95","u95")] <- a

a <- apply(postDraws$lambda[,1]-postDraws$lambda[,-1], 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_1", tab$parm)
tab[i,c("postmean","l95","u95")] <- a

a <- apply(postDraws$lambda[,5]-postDraws$lambda[,-5], 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_5", tab$parm)
tab[i,c("postmean","l95","u95")] <- a

# increasing uncertainty with distance
a <- apply(postDraws$lambda[,1]-postDraws$lambda[,-1], 2, sd)
i <- grep("^\\\\lambda_1", tab$parm)
tab[i,c("postSD")] <- a

a <- apply(postDraws$lambda[,5]-postDraws$lambda[,-5], 2, sd)
i <- grep("^\\\\lambda_5", tab$parm)
tab[i,c("postSD")] <- a

a <- apply(postDraws$s, 2, sd)
i <- grep("^s_", tab$parm)
tab[i,c("postSD")] <- a


## Table D3, left panel --------------------------------------------------------
a <- apply(postDraws$lambda, 1, rank) |>
  apply(1, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(1) |> t()
i <- grep("lambda_", tabRanks$parm)
tabRanks[i,c("postmean","l95","u95")] <- a

# "If we fixed the value of λ1, values of Var(λ−1 | λ1, data) range from around 0.05 for Officers 2 and 3 up to 0.19 for Officer 10."
Sigma <- cov(postDraws$lambda)
schurComp <- sapply(1:nOff, function(i) diag(Sigma) - Sigma[,i]^2/Sigma[i,i]) |> 
  t() |> zapsmall()
schurComp[1,] |> round(2)

# show increasing trend in conditional variance - not used in article
i <- 1
a <- schurComp[1,]
par(mai=0.02+c(0.8,0.8,0.4,0.4))
plot((1:10), a,
     pch=16,
     ylim=c(0,0.2),
     xlab="Peer Officer ID",
     ylab=bquote("Var(" ~ lambda[-1] ~ "|" ~ lambda[1] ~ ")")) 



# 10 officers in a chained network, 24 incidents/officer ----------------------
load("output/mcmcSampOff10small.RData")
nOff <- n_distinct(d0$idOff)
d0$idOff <- d0$idOff + 1 # back to R indexing

# number of incidents
n_distinct(d0$id)
postDraws <- list(lambda = res$draws[,1:nOff],
                  sDelta = exp(res$draws[,nOff + 1:2]))
postDraws$s <- 1 + apply(postDraws$sDelta, 1, cumsum) |> t()
colnames(postDraws$s) <- c("s2","s3")


## Table 3 & Table D2, right panel --------------------------------------------
a <- apply(postDraws$s, 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("s_", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a

a <- apply(postDraws$lambda[,1]-postDraws$lambda[,-1], 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_1", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a

a <- apply(postDraws$lambda[,5]-postDraws$lambda[,-5], 2, 
           function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()
i <- grep("^\\\\lambda_5", tab$parm)
tab[i,c("postmeanSmall","l95small","u95small")] <- a

a <- apply(postDraws$lambda[,1]-postDraws$lambda[,-1], 2, sd)
i <- grep("^\\\\lambda_1", tab$parm)
tab[i,c("postSDsmall")] <- a

a <- apply(postDraws$lambda[,5]-postDraws$lambda[,-5], 2, sd)
i <- grep("^\\\\lambda_5", tab$parm)
tab[i,c("postSDsmall")] <- a

a <- apply(postDraws$s, 2, sd)
i <- grep("^s_", tab$parm)
tab[i,c("postSDsmall")] <- a


## Table D3, right panel ------------------------------------------------------
a <- apply(postDraws$lambda, 1, rank) |>
  apply(1, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(1) |> t()
i <- grep("lambda_", tabRanks$parm)
tabRanks[i,c("postmeanSmall","l95small","u95small")] <- a

# not discussed in article, but show increasing trend in conditional variance
Sigma <- cov(postDraws$lambda)
schurComp <- sapply(1:nOff, function(i) diag(Sigma) - Sigma[,i]^2/Sigma[i,i]) |> 
  t() |> zapsmall()
schurComp[1,] |> round(2)
schurComp[5,] |> round(2)

# show increasing trend in conditional variance - not used in article
i <- 1
a <- schurComp[1,]
par(mai=0.02+c(0.8,0.8,0.4,0.4))
plot((1:10), a,
     pch=16,
     #ylim=c(0,0.4),
     xlab="Peer Officer ID",
     ylab=bquote("Var(" ~ lambda[-1] ~ "|" ~ lambda[1] ~ ")"))


## Table 3 ---------------------------------------------------------------------
tab |>
  mutate(interval = paste0("(",trimws(format(l95,nsmall=2)),", ",
                               trimws(format(u95,nsmall=2)),")"),
         intervalSmall = paste0("(",trimws(format(l95small,nsmall=2)),", ",
                                    trimws(format(u95small,nsmall=2)),")"),
         parm = paste0("$",parm,"$")) |>
  relocate(interval, .before=4) |>
  relocate(intervalSmall, .before=9) |>
  select(-l95, -u95, -l95small, -u95small) |>
  xtable() |>
  print(sanitize.text.function = identity,
        include.rownames = FALSE)

## Table D3 -------------------------------------------------------------------
tabRanks |>
  mutate(interval = paste0("(",trimws(format(l95,nsmall=1)),", ",
                           trimws(format(u95,nsmall=1)),")"),
         intervalSmall = paste0("(",trimws(format(l95small,nsmall=1)),", ",
                                trimws(format(u95small,nsmall=1)),")"),
         parm = paste0("$",parm,"$")) |>
  relocate(interval, .before=4) |>
  relocate(intervalSmall, .before=8) |>
  select(-l95, -u95, -l95small, -u95small) |>
  xtable(digits = 1) |>
  print(sanitize.text.function = identity,
        include.rownames = FALSE)


# Section 4.3
# "The model estimates that with 0.86 posterior probability Officer 8 is among the two officers with the smallest $\lambda$ (in truth it is the smallest). The posterior probability that Officer 5 has the largest $\lambda$ is 0.27 and a 0.71 probability that it is among the top two (in truth it is the largest)."
R <- apply(postDraws$lambda, 1, rank)
#   is 8 the smallest?
apply(R, 2, function(x) x[8] <= 2) |> mean()
#   is 5 the largest?
apply(R, 2, function(x) x[5] == 10) |> mean()
apply(R, 2, function(x) x[5] >= 9)  |> mean()



# Seattle ---------------------------------------------------------------------
load("data/dataSPD.RData")

# Section 5 stats
# number of officers
nOff <- n_distinct(d$idOff)
nOff

# number of incidents
n_distinct(d$id)
# number of uses-of-force
sum(d$y!=1)
# number of witness officers
sum(d$y==1)

# distribution of force levels
d |> count(y) |> mutate(pct=n/sum(n))


chains <- vector("list", length = 4)
for(iChain in 1:length(chains))
{
  message("Loading chain ", iChain)
  load(paste0("output/mcmcSampSPDchain",iChain,".RData"))

  postDraws <- list(lambda = resSPD$draws[,1:nOff],
                    sDelta = exp(resSPD$draws[,nOff+1:2]))
  colnames(postDraws$sDelta) <- c("sDiff1.2","sDiff2.3")

    # compute s from sDelta
  postDraws$s <- 1 + apply(postDraws$sDelta, 1, cumsum) |> t()
  colnames(postDraws$s) <- c("s2","s3")

  chains[[iChain]] <- cbind(lambda = postDraws$lambda,
                            s      = postDraws$s,
                            sDelta = postDraws$sDelta)
}

burnin <- 200
chains <- lapply(chains, function(x) x[-(1:burnin),])
postDraws <- do.call(rbind, chains)

## Posterior means of s2 and s3 in main text ----------------------------------
postDraws[,1504:1507] |> 
  apply(2, function(x) c(mean(x), quantile(x,prob=c(0.025,0.975)))) |>
  round(2) |> t()


# Find officers' peer groups, described in Section 3.4
#   compute rank of first value relative to the rest
rank1 <- function(x)
{
  1+sum(x[1]>x[-1])
}

# Schur complement calculation
#   schurComp = Var(lambda[-i]|lambda[i])
#   formula from Section 3.4
Sigma <- cov(postDraws[,1:nOff])
schurComp <- sapply(1:nOff, 
                    function(i) diag(Sigma) - Sigma[,i]^2/Sigma[i,i]) |> 
  t() |> zapsmall()

# about 20 minutes
#   use 4 cores... can be memory intensive
#   can run load("output/resultsSPD.RData") to skip this step
plan(multisession, workers = 4)
registerDoFuture()
handlers("cli")

with_progress({
  p <- progressor(steps = nOff)
  
  resultsSchur <- foreach(iID=1:nOff, .combine=rbind) %dopar%
  {
    p(sprintf("Processing officer %d", iID))
    
    # find peers
    iPeers <- (1:nOff)[schurComp[iID,] < 0.3] |> setdiff(iID)
    # find rank of officer iID relative to officers iPeers
    R <- apply(postDraws[, c(iID,iPeers)], 1, rank1)
    
    lam0 <- postDraws[,c(iID,iPeers)]
    lam0 <- lam0[,1] - rowMeans(lam0)

    res <- data.frame(idOff=iID,
                      nOff = length(iPeers),
                      nInc = sum(d$idOff==iID),
                      pRankTop5pct = mean(R >= 0.95*(1+length(iPeers))),
                      pRankBot5pct = mean(R <= 0.05*(1+length(iPeers))),
                      lamMean = mean(lam0),
                      lam025 = quantile(lam0, prob=0.025) |> as.numeric(),
                      lam975 = quantile(lam0, prob=0.975) |> as.numeric())
    res[,c("force0","force1","force2","force3")] <- table(factor(d$y)[d$idOff==iID])
    return(res)
  }
})
save(resultsSchur, file="output/resultsSPD.RData")


## Table 4 --------------------------------------------------------------------
resultsSchur |>
  select(idOff,nOff,nInc,force0,force1,force2,force3,
         pRankTop5pct,
         lamMean,lam025,lam975) |>
  filter(pRankTop5pct > 0.80) |>
  arrange(desc(pRankTop5pct)) |>
  xtable(digits = c(0,0,0,0,0,0,0,0,2,2,2,2)) |>
  print(include.rownames=FALSE)

# translate lambda to propensity to escalate
exp(1.99*diff(c(s0=0, s1=1, s2=1.05, s3=1.13)))

# low-end outliers, not shown in article
resultsSchur |>
  select(idOff,nOff,nInc,force0,force1,force2,force3,
         pRankBot5pct,
         lamMean,lam025,lam975) |>
  filter(nInc>5 & pRankBot5pct > 0.90) |>
  arrange(desc(pRankBot5pct)) |>
  xtable(digits = c(0,0,0,0,0,0,0,0,2,2,2,2)) |>
  print(include.rownames=FALSE)


# beeswarm plot, not shown in main article
#   shown as comparison to placebo tests in Figure F1
df <- resultsSchur |>
  mutate(x      = lambdaDiffPeerMean,
         group  = case_when(pRankTop5pct  > 0.8 ~ "Top 5%",
                            pRankBot5pct  > 0.8 ~ "Bottom 5%",
                            TRUE                ~ "Middle")) |>
  select(x,group)

ggplot(df, aes(x = "", y = x, colour = group)) +
  geom_quasirandom(width = 0.35, size = 1.9, alpha = 0.8) +
  scale_colour_manual(values = c("Top 5%"   = "#d95f02",
                                 "Bottom 5%" = "#0173b2",
                                 "Middle"          = "grey70"),
                      breaks = c("Top 5%", "Bottom 5%"),
                      name   = NULL) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0.15)) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0.4)) +
  labs(y = expression(lambda[i] - mean(lambda[peers])),
       x = NULL) +
  theme_minimal() +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 


## Figure 6 -------------------------------------------------------------------
iID <- 1006
par(mai=0.02+c(0.8,0.4,0.4,0.4))
hist(schurComp[iID,], axes=FALSE, ylab="", main="", 
     xlab=bquote("Var(" ~ lambda[-1006] ~ "|" ~ lambda[1006] ~ ")"),
     breaks = 30)
axis(1)

sum(schurComp[iID,]<=0.3)


# get all 1006 connected officers including incidents not involving 1006
a <- d |>
  right_join(d |> # all officers at UoF incidents with Officer 1006
               filter(id %in% d$id[d$idOff==iID]) |>
               select(idOff) |>
               distinct())


n_distinct(a$id)
n_distinct(a$idOff)
nrow(a)

# get edge counts
b <- split(a$idOff, a$id) |>
  lapply(function(x) expand.grid(x,x)) |>
  bind_rows() |>
  filter(Var1 < Var2) |>
  count(Var1,Var2) |> 
  data.frame()
net <- graph_from_data_frame(b, directed=FALSE)
col <- rep("white", length(unique(a$idOff)))
col[which(names(V(net))==iID)] <- "lightblue"

# highlight nodes with high conditional variance
i <- V(net) |> names() |> as.numeric()
i <- i[schurComp[iID,i] > 0.6]
border.col <- rep("black", length(unique(a$idOff)))
border.col[which(names(V(net))==i)] <- "red"

## Figure 5 -------------------------------------------------------------------
# Note: plot.igraph() layout changes each time generated
par(mai=c(0,0,0,0))
plot(net,
     edge.width=b$n,
     vertex.color=col,
     vertex.frame.color=border.col)
# get number of edges
igraph::E(net)

# interactive version (not used in article)
# tkplot(net, edge.width=b$n, vertex.col=col)
# tk_coords(1)

# conditional variances discussed in Section 5
schurComp[iID, c(536, 993, 1018, 1033, 232, 435)] |> round(2)

## Figure 7 --------------------------------------------------------------------
## create set O_i, with t<0.3
iPeers <- which(schurComp[iID,] < 0.3)
R <- apply(postDraws[,iPeers], 1, rank)
Rpostint <- apply(R, 1, quantile, prob=c(0.025,0.5,0.975))

j <- order(Rpostint[2,])
iZoom <- tail(1:length(iPeers), 60)
par(mai=0.02+c(0.8,0.8,0.4,0.4))
plot((1:ncol(Rpostint))[iZoom], Rpostint[2,j][iZoom], pch=16,
     ylab="Rank distribution",
     xlab="Officer ID",
     axes=FALSE,
     ylim=c(0,700))
x <- (1:length(j))[iZoom][iPeers[j][iZoom]==iID]
rect(x-0.5,par()$usr[3],x+0.5,par()$usr[4], col="gray", border=NA)

box()
axis(2)
axis(1,at=(1:length(j))[iZoom], 
     labels = iPeers[j][iZoom], 
     las=2, cex.axis=0.7)
for(k in (1:ncol(Rpostint))[iZoom])
{
  lines(c(k,k), Rpostint[c(1,3),j[k]], col="red",
        lwd = ifelse(iPeers[j[k]]==1006, 4, 1))
  points(k, Rpostint[2,j[k]], pch=16)
}
abline(h = length(iPeers)*0.95, col="gray")

# in how many cases did Officer 1006 have the highest rank?
# "Officer 1006 was the officer who used the most serious form of force in 8 of the 9 incidents for which they were on the scene."
a <- d$id[d$idOff==iID]
d |> filter(id %in% a) |>
  group_by(id) |>
  summarize(y[idOff==iID]==max(y))

# "In seven of those force incidents, Officer 1006 was the only officer to use force."
d |> filter(id %in% a) |>
  group_by(id) |>
  summarize(all(y[idOff!=iID]==1))

# which officers used a lot of force but are not flagged as outliers
#   demonstrates that traditional threshold based methods come to a different 
#   result. They ignore context

# "For example, Officer 811, who had 15 Level 2 and three Level 3 uses of force, does not appear in Table 5. Instead, the posterior probability that Officer 811 escalates less than Officer 1006, who has never exceeded Level 1 force, is 97%."
d |> 
  filter(y==4) |> 
  count(idOff) |>
  arrange(desc(n)) |>
  head(10) |>
  left_join(resultsSchur, by=join_by(idOff))

mean(postDraws[,811] < postDraws[,1006])


# make full SPD graph
b <- split(d$idOff, d$id) |>
  lapply(function(x) expand.grid(x,x)) |>
  bind_rows() |>
  filter(Var1 < Var2) |>
  count(Var1,Var2) |> 
  data.frame()
net <- graph_from_data_frame(b, directed=FALSE)

# "All but two of the 1,503 officers are members of a connected network"
net |> components() |> purrr::pluck("csize")
# "...with an average distance of 2.1 and a maximum distance of 5 separating officers,"
i <- net |> components() |> purrr::pluck("membership")
net <- net |> delete_vertices(which(i==2))
net |> distances() |> table() 
net |> distances() |> mean() |> round(1)

