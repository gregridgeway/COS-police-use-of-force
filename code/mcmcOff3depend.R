
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(viridis)

RcppParallel::setThreadOptions(numThreads = 8)

# load C++ code for conditional likelihood calcs and acceptance probability
source("code/mcmcConditionalStereotype.R")

# Appendix E1 -----------------------------------------------------------------
set.seed(20010618)

nIncidents <- 1000
nOff <- sample(2:3, nIncidents, replace = TRUE)

lambda <- -c(0.5, 1, 2)
s <- c(0, 1, 1.5, 2)
nTotOfficers <- 3

thetaInit <- c(rep(0, nTotOfficers), 0.5, 0.5)

design <- expand.grid(delta=seq(-1,1,length=41),
                      ranks=NA,
                      diffs=NA,
                      lambdaMean=NA)

for(iIter in 1:nrow(design))
{
  message("Iteration: ", iIter, " of ", nrow(design))
  # randomly assign officers to incidents
  d <- data.frame(id=rep(1:nIncidents, nOff),
                  idOff=unlist(sapply(nOff, function(x) sample(1:3, x))))
  
  d$lambda <- lambda[d$idOff]
  d <- d |> 
    group_by(id) |>
    # if idOff==1 is present, add delta to all other officers
    mutate(lambda = lambda + 
             design$delta[iIter]*any(idOff==1)*(idOff!=1)) |>
    group_by(id,idOff) |>
    mutate(y = sample(1:4, 1, prob=exp(s*lambda))) |>
    ungroup()
  
  # keep only incidents that have more than one type of force used
  d <- d |>
    right_join(d |> 
                 group_by(id) |>
                 summarize(nUniqueY = n_distinct(y)) |>
                 filter(nUniqueY > 1) |>
                 select(id)) |>
    mutate(id = as.integer(as.factor(id))) |> # renumber 1, 2, 3, ...
    filter(id <= 300)
  
  # make 0-based for C++
  d$y     <- d$y - 1
  d$idOff <- d$idOff - 1
  
  system.time({
    res <- mcmcOrdinalStereotype(
      d, 
      lambda0 = thetaInit[1:nTotOfficers],
      sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
      nIter   = 200000,
      thin    = 20, 
      sdProp  = 0.2)
  })
  
  design$ranks[iIter] <- 
    res$draws[,1:nTotOfficers] |>
    apply(1, rank) |>
    apply(1, mean) |>
    list()
  
  design$diffs[iIter] <- 
    (res$draws[,2:nTotOfficers] - res$draws[,1]) |>
    apply(2, function(x) 
      quantile(x, prob=c(0.5, 0.025,0.975))) |>
    list()
  
  design$lambdaMean[iIter] <- d |> 
    group_by(idOff) |> 
    summarize(lambdamean = mean(lambda)) |> 
    pull(lambdamean) |>
    list()
  
  print(design[iIter,])
}

save(design, file="output/mcmcOff3depend.RData", compress = TRUE)

# diff plots
a    <- sapply(design$diffs, function(x) x[1,])  # medians
a025 <- sapply(design$diffs, function(x) x[2,])
a975 <- sapply(design$diffs, function(x) x[3,])

param_levels <- c("lambda2 - lambda1", "lambda3 - lambda1")

plot_df <- data.frame(
  delta = rep(design$delta, times = 2),
  parameter = factor(rep(param_levels, each = length(design$delta)), levels = param_levels),
  median = c(a[1,], a[2,]),
  lower = c(a025[1,], a025[2,]),
  upper = c(a975[1,], a975[2,]))

ref_df <- data.frame(
  parameter = factor(param_levels, levels = param_levels),
  true_diff = c(lambda[2] - lambda[1], lambda[3] - lambda[1]))

curve_df <- design$lambdaMean |>
  do.call(rbind, args=_) |>
  as.data.frame() |>
  mutate(V2=V2-V1, V3=V3-V1,
         delta = design$delta) |>
  select(-V1) |>
  pivot_longer(cols = c(1, 2),
               names_to = "parameter",
               values_to = "value") |>
  mutate(parameter = recode(parameter,
                            `V2` = "lambda2 - lambda1",
                            `V3` = "lambda3 - lambda1"))

label_expressions <- setNames(
  list(bquote(lambda[2] - lambda[1]), bquote(lambda[3] - lambda[1])),
  param_levels)

## Figure E1 ------------------------------------------------------------------
ggplot(plot_df, aes(x = delta, y = median, color = parameter, 
                    fill = parameter)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_line(data = curve_df, aes(x = delta, y = value, color = parameter),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("#1b9e77", "#d95f02"), 
                     labels = label_expressions) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02"), 
                    labels = label_expressions) +
  labs(x = expression(delta),
       y = expression(lambda[i] - lambda[1]),
       color = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())


# rank plots
param_labels <- c("lambda1", "lambda2", "lambda3")
ranks_df <- bind_rows(design$ranks) %>%
  mutate(delta = design$delta) %>%
  setNames(c(param_labels, "delta")) %>%
  pivot_longer(-delta, names_to = "parameter", values_to = "rank") %>%
  mutate(parameter = factor(parameter, levels = param_labels))

ref_lines <- data.frame(parameter = param_labels,
                        rank_value = rank(lambda))

curve_df <- design$lambdaMean |>
  do.call(rbind, args=_) |>
  apply(1, rank) |>
  t() |>
  as.data.frame() |>
  mutate(delta = design$delta) |>
  pivot_longer(-delta, names_to = "parameter", values_to = "value") |>
  mutate(parameter = recode(parameter,
                            V1 = "lambda1",
                            V2 = "lambda2",
                            V3 = "lambda3"))

label_expressions <- setNames(
  lapply(1:3, function(i) bquote(lambda[.(i)])),
  param_labels)

## Figure E2 ------------------------------------------------------------------
ggplot(ranks_df, aes(x = delta, y = rank, color = parameter)) +
  geom_line(linewidth = 1.2) +
  geom_line(data = curve_df, aes(x = delta, y = value, color = parameter),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = brewer.pal(3, "Set2"),
                     labels = label_expressions) +
  labs(x = expression(delta),
       y = "Ranks",
       color = "Parameter") +
  coord_cartesian(ylim = c(0, 4)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())





# Appendix E2 -----------------------------------------------------------------

set.seed(20010618)
nIncidents <- 10000
nTotOfficers <- 30
# number of officers per incident
nOff <- sample(2:5, nIncidents, replace = TRUE)

lambda <- seq(-2,2, length.out = nTotOfficers)
s <- c(0, 1, 1.05, 1.15)
thetaInit <- c(rep(0, nTotOfficers), 0.5, 0.5)

design <- expand.grid(delta=seq(-1,1,length=21),
                      ranks=NA,
                      diffs=NA,
                      lambdaMean=NA)

for(iIter in 1:nrow(design))
{
  message("Iteration: ", iIter, " of ", nrow(design))
  # randomly assign officers to incidents
  d <- data.frame(id=rep(1:nIncidents, nOff),
                  idOff=unlist(sapply(nOff, function(x) sample(1:nTotOfficers, x))))
  
  d$lambda0 <- lambda[d$idOff]
  d <- d |> 
    group_by(id) |>
    # if idOff==1 is present, subtract delta from all other officers
    #    Officer 1 de-escalates others
    mutate(lambda0 = lambda0 - 
             design$delta[iIter]*any(idOff==1)*(idOff!=1),
           # if idOff==30 is present, add delta to all other officers
           #    Officer 30 escalates others
           lambda0 = lambda0 + 
             design$delta[iIter]*any(idOff==30)*(idOff!=30)) |>
    group_by(id,idOff) |>
    mutate(y = sample(1:4, 1, prob=exp(s*lambda0))) |>
    ungroup()
  
  # keep only incidents that have more than one type of force used
  d <- d |>
    right_join(d |> 
                 group_by(id) |>
                 summarize(nUniqueY = n_distinct(y)) |>
                 filter(nUniqueY > 1) |>
                 select(id)) |>
    mutate(id = as.integer(as.factor(id)))   # renumber 1, 2, 3, ...
  
  # make 0-based for C++
  d$y     <- d$y - 1
  d$idOff <- d$idOff - 1
  
  system.time({
    res <- mcmcOrdinalStereotype(
      d, 
      lambda0 = thetaInit[1:nTotOfficers],
      sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
      nIter   = 400000,
      thin    = 20, 
      sdProp  = 0.035)
  })
  
  design$ranks[iIter] <- 
    res$draws[,1:nTotOfficers] |>
    apply(1, rank) |>
    apply(1, mean) |>
    list()
  
  design$diffs[iIter] <- 
    (res$draws[,2:nTotOfficers] - res$draws[,1]) |>
    apply(2, function(x) 
      quantile(x, prob=c(0.5, 0.025,0.975))) |>
    list()
  
  design$lambdaMean[iIter] <- d |> 
    group_by(idOff) |> 
    summarize(lambdamean = mean(lambda0)) |> 
    pull(lambdamean) |>
    list()
  
  print(design[iIter,])
}

save(design, file="output/mcmcOff30depend.RData", compress = TRUE)

# rank plot
param_labels <- paste0("lambda", 1:nTotOfficers)
ranks_df <- bind_rows(design$ranks) |>
  mutate(delta = design$delta) |>
  setNames(c(param_labels, "delta")) |>
  pivot_longer(-delta, names_to = "parameter", values_to = "rank") |>
  mutate(parameter = factor(parameter, levels = param_labels),
         line_size = ifelse(parameter %in% c("lambda1", "lambda30"), 1.8, 1.0))

ref_lines <- data.frame(parameter = param_labels,
                        rank_value = rank(lambda))

curve_df <- design$lambdaMean |>
  do.call(rbind, args=_) |>
  apply(1, rank) |>
  t() |>
  as.data.frame() |>
  mutate(delta = design$delta) |>
  pivot_longer(-delta, names_to = "parameter", values_to = "value") |>
  mutate(parameter = gsub("^V", "lambda", parameter),
         parameter = factor(parameter, levels = param_labels))

label_expressions <- setNames(
  lapply(1:nTotOfficers, function(i) bquote(lambda[.(i)])),
  param_labels)

palette_colors <- viridis(n = nTotOfficers, option = "D")

## Figure E3 ------------------------------------------------------------------
ggplot(ranks_df,
       aes(x = delta, y = rank, color = parameter)) +
  geom_line(aes(linewidth = line_size)) +
  geom_line(data = curve_df,
            aes(x = delta, y = value, color = parameter),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = palette_colors,
                     labels = label_expressions) +
  scale_linewidth_identity(guide = "none") +
  labs(x = expression(delta),
       y = "Ranks",
       color = NULL) +
  annotate("label", x = 0.45, y = 3,
           label = "Officer 1 de-escalates peers",
           size = 4, hjust = 0,
           label.size = 0.3, label.r = unit(0.15, "lines"),
           fill = "white") +
  annotate("label", x = 0.45, y = 29,
           label = "Officer 30 escalates peers",
           size = 4, hjust = 0,
           label.size = 0.3, label.r = unit(0.15, "lines"),
           fill = "white") +
  coord_cartesian(ylim = c(0, nTotOfficers + 1)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))



# Appendix E3 -----------------------------------------------------------------

#  "first" officer use of force escalates others
set.seed(20010618)
nIncidents <- 10000
nTotOfficers <- 3
# number of officers per incident
nOff <- sample(2:3, nIncidents, replace = TRUE)

# simulate with these parameters
lambda <- c(-1, -0.5, 0)
s <- c(0, 1, 1.05, 1.15)
thetaInit <- c(lambda, 0.5, 0.5)

design <- data.frame(delta=seq(-1,1,length=21),
                     ranks=NA,
                     diffs=NA,
                     lambdaMean=NA)

# randomly assign officers to incidents
for(iIter in 1:nrow(design))
{
  message("Iteration: ", iIter, " of ", nrow(design))

  d <- data.frame(id=rep(1:nIncidents, nOff),
                  idOff=unlist(sapply(nOff, function(x) 
                    sample(1:nTotOfficers, x)))) |> 
    mutate(lambda0 = lambda[idOff]) |>
    mutate(y0 = sample(1:4, 1, prob = exp(s * lambda0)),
           .by = c(id, idOff)) |>
    mutate(lambda00 = lambda0 + 
             (row_number() != 1) *          # not the first row in this id
             (first(y0)    > 1) *           # first officer used force
             design$delta[iIter],
           .by = id) |>
    mutate(y00 = sample(1:4, 1, prob = exp(s * lambda00)),
           .by = c(id, idOff)) |>
    mutate(y = if_else((row_number() == 1) | (first(y0)==1), y0, y00),
           .by = id)

  # keep only incidents that have more than one type of force used
  #   such incidents contribute 1 to likelihood
  d <- d |>
    right_join(d |> 
                 group_by(id) |>
                 summarize(nUniqueY = n_distinct(y)) |>
                 filter(nUniqueY > 1) |>
                 select(id)) |>
    mutate(id = as.integer(as.factor(id)))   # renumber 1, 2, 3, ...
  
  # make 0-based for C++
  d$y     <- d$y - 1
  d$idOff <- d$idOff - 1
  
  system.time({
    res <- mcmcOrdinalStereotype(
      d |> select(id,idOff,y), 
      lambda0 = thetaInit[1:nTotOfficers],
      sDiff0  = thetaInit[(nTotOfficers+1):(nTotOfficers+2)],
      nIter   = 400000,
      thin    = 20, 
      sdProp  = 0.05)
  })

  design$ranks[iIter] <- 
    res$draws[,1:nTotOfficers] |>
    apply(1, rank) |>
    apply(1, mean) |>
    list()
  
  design$diffs[iIter] <- 
    (res$draws[,2:nTotOfficers] - res$draws[,1]) |>
    apply(2, function(x) 
      quantile(x, prob=c(0.5, 0.025,0.975))) |>
    list()
  
  design$lambdaMean[iIter] <- d |> 
    group_by(idOff) |> 
    summarize(lambdamean = mean(lambda00)) |> 
    pull(lambdamean) |>
    list()
  
  print(design[iIter,])
}

save(design, file="output/mcmcOff3firstOff.RData", compress = TRUE)

# diff plots
a    <- sapply(design$diffs, function(x) x[1,])  # medians
a025 <- sapply(design$diffs, function(x) x[2,])
a975 <- sapply(design$diffs, function(x) x[3,])

curve_df <- design$lambdaMean |>
  do.call(rbind, args = _) |>
  as.data.frame() |>
  mutate(across(V2:V3, ~ .x - V1),
         delta = design$delta) |>
  select(-V1) |>
  pivot_longer(cols = starts_with("V"),
               names_to = "parameter",
               values_to = "value") |>
  mutate(index = as.integer(gsub("V", "", parameter)),
         parameter = paste0("lambda", index, " - lambda1"))

param_levels <- unique(curve_df$parameter)


plot_df <- data.frame(
  delta = rep(design$delta, times = nTotOfficers-1),
  parameter = factor(rep(param_levels, 
                         each = length(design$delta)), 
                     levels = param_levels),
  median = as.vector(t(a)),
  lower = as.vector(t(a025)),
  upper = as.vector(t(a975)))

ref_df <- data.frame(
  parameter = factor(param_levels, levels = param_levels),
  true_diff = lambda[2:length(lambda)] - lambda[1]
)

label_expressions <- setNames(
  lapply(2:3, function(i) bquote(lambda[.(i)] - lambda[1])),
  paste0("lambda", 2:3, " - lambda1")
)


palette_colors <- viridis(n = length(param_levels), option = "D")

## Figure E4 -------------------------------------------------------------------
ggplot(plot_df, aes(x = delta, y = median, color = parameter, 
                    fill = parameter)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_line(data = curve_df, aes(x = delta, y = value, color = parameter),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = palette_colors, 
                     labels = label_expressions[param_levels]) +
  scale_fill_manual(values = palette_colors, 
                    labels = label_expressions[param_levels]) +
  labs(x = expression(delta),
       y = expression(lambda[i] - lambda[1]),
       color = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

