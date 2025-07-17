library(dplyr)
library(doParallel)
library(future)
library(doFuture)
library(progressr)
library(ggplot2)
library(ggbeeswarm)


load("data/dataSPD.RData")

nOff <- n_distinct(d$idOff)

# rank of first value among all values
rank1 <- function(x)
{
  1+sum(x[1]>x[-1])
}


postDraws <- vector("list", length = 5)
for(iPerm in 1:length(postDraws))
{
  message(paste("Loading chain ", iPerm))
  load(paste0("output/mcmcSampSPDPermutechain",iPerm,".RData"))

  Sigma <- cov(resSPD$draws[,1:nOff])
  schurComp <- sapply(1:nOff, 
                      function(i) diag(Sigma) - Sigma[,i]^2/Sigma[i,i]) |> 
    t() |> zapsmall()

  # ~10 minutes per chain
  #   use 4 cores... can be memory intensive
  plan(multisession, workers = 4)
  registerDoFuture()
  handlers("cli")
  
  with_progress({
    p <- progressor(steps = nOff)
    
    resultsSchur <- foreach(iID=1:nOff, .combine=rbind) %dopar%
    {
      p(sprintf("Processing officer %d", iID))
      iPeers <- (1:nOff)[schurComp[iID,] < 0.3] |> setdiff(iID)
      R <- apply(resSPD$draws[,c(iID,iPeers)], 1, rank1)
      
      res <- data.frame(idOff=iID,
                        nOff = length(iPeers),
                        nInc = sum(d$idOff==iID),
                        pRankTop5pct = mean(R >= 0.95*(1+length(iPeers))),
                        pRankBot5pct = mean(R <= 0.05*(1+length(iPeers))),
                        lambdaDiffPeerMean = mean(resSPD$draws[,iID]) - 
                          mean(resSPD$draws[,c(iID,iPeers)]))
      res[,c("force0","force1","force2","force3")] <- 
        table(factor(d$y)[d$idOff==iID])
      return(res)
    }
  })
  postDraws[[iPerm]]$resultsSchur <- resultsSchur
}
save(postDraws, file = "output/resultsSPDpermute.RData", compress = TRUE)

load("output/resultsSPD.RData")
resOriginal <- resultsSchur

res <- data.frame(
  type=rep(c("True","Placebo 1","Placebo 2",
             "Placebo 3","Placebo 4","Placebo 5"), 
           each=nOff),
  pTop = c(resOriginal$pRankTop5pct, 
           sapply(postDraws, function(x) x$resultsSchur$pRankTop5pct)),
  pBot = c(resOriginal$pRankBot5pct, 
           sapply(postDraws, function(x) x$resultsSchur$pRankBot5pct)),
  lambdaDiff = c(resOriginal$lambdaDiffPeerMean, 
                 sapply(postDraws, function(x) x$resultsSchur$lambdaDiff)))
  
df <- res |>
  mutate(x    = lambdaDiff,
         flag = case_when(pTop  > 0.8 ~ "Top 5%",
                          pBot  > 0.8 ~ "Bottom 5%",
                          TRUE  ~ "Middle")) |>
  select(type,x,flag)

# Figure F1 -------------------------------------------------------------------
ggplot(df, aes(x = "", y = x, colour = flag)) +
  geom_quasirandom(width = 0.35,
                   size  = 1.9,
                   alpha = 0.8) +
  facet_wrap(~type,
             ncol = 1,
             scales = "fixed") +
  scale_colour_manual(values = c("Top 5%"    = "#d95f02",
                                 "Bottom 5%" = "#0173b2",
                                 "Middle"    = "grey70"),
                      breaks = c("Top 5%", "Bottom 5%"),
                      name   = NULL) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0.15)) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0.40)) +
  labs(y = expression(E(lambda[i]) - E(lambda[peers] ~ "|" ~ italic(O)[i])),
       x = NULL) +
  theme_minimal() +
  theme(strip.text   = element_text(face = "bold"),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()
