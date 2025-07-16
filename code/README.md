CODE - Conditional Ordinal Stereotype Model to Estimate Police Officers' Propensity to Escalate Force
================

- [`conditional_stereotype.cpp`](conditional_stereotype.cpp) contains the C++ code that parallelizes the conditional likelihood calculation by incident

- [`mcmcConditionalStereotype.R`](mcmcConditionalStereotype.R) contains the Metropolis algorithm that samples from the posterior distribution of the conditional stereotype model parameters. Any script calling this script should first call `RcppParallel::setThreadOptions(numThreads = 8)`, editing `numThreads` to be appropriate for the computing environment.

- [`mcmcOff3.R`](mcmcOff3.R) Simulation #1

- [`mcmcOff3+3.R`](mcmcOff3+3.R) Simulation #2

- [`mcmcOff10.R`](mcmcOff10.R) Simulation #3

- [`mcmcSPD.R`](mcmcSPD.R) Analysis of use-of-force data from Seattle PD

- [`post mcmc analysis.R`](post mcmc analysis.R) Analysis of the posterior samples from the three simulations and the Seattle PD data. Labels blocks of code that connect to specific sections, tables, and figures in the article