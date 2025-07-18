Conditional Ordinal Stereotype Model to Estimate Police Officers' Propensity to Escalate Force
================

# Overview
This repository contains the code and data associated with "A Conditional Ordinal Stereotype Model to Estimate Police Officers' Propensity to Escalate Force." The article introduced a conditional ordinal stereotype model for estimating a police officer's latent propensity to escalate to higher-severity force options during an encounter. Propensity to escalate is the likelihood of selecting a more serious force category than peer officers confronting the same circumstances. The associated conditional likelihood depends solely on data from the times and places where multiple officers use various forms of force. All shared situational and environmental features cancel out of the conditional likelihood, eliminating the need to observe or model them explicitly.

# Core analytical components

Code for the conditional stereotype model is a blend of R and C++. The C++ code ([`code/conditional_stereotype.cpp`](code/conditional_stereotype.cpp)) contains the C++ code that parallelizes the conditional likelihood calculation by incident. For incidents with two officers, there is a simple, direct calculation. For incidents with three to seven officers, there is a no-repeats version of Heap's algorithm O(m!). For incidents with eight or more officers, there is a multi-dimensional discrete Fourier transform O(m^J). This arrangement results in the most efficient conditional likelihood calculation. The R script ([`code/mcmcConditionalStereotype.R`](code/mcmcConditionalStereotype.R)) contains the Metropolis algorithm that samples from the posterior distribution of the conditional stereotype model parameters.

# Simulations and Seattle PD data

The [`code`](code/) folder contains three scripts associated with the three simulated examples in the article. Each of the R scripts are self-contained, meaning that there is no other setup needed when running these scripts. The scripts `source()` the code from ([`code/mcmcConditionalStereotype.R`](code/mcmcConditionalStereotype.R)), which uses `Rcpp::sourceCpp()` to compile the needed C++ code and make the conditional likelihood and Metropolis acceptance probability calculation available in R. 

Any script using these should first call `RcppParallel::setThreadOptions(numThreads = 8)`, editing `numThreads` to be appropriate for the computing environment. The three simulations take a few minutes to run. The Seattle analysis can take 2-3 days depending on the computing environment.

- Simulation #1: Demonstrates proof-of-concept in a well-connected three officer network. [`code/mcmcOff3.R`](code/mcmcOff3.R)
- Simulation #2: Shows the challenge of analyzing disconnected use-of-force networks and that the conditional variance helps to partition the network into peers and non-peers. [`code/mcmcOff3+3.R`](code/mcmcOff3+3.R)
- Simulation #3: Shows that information can flow through lengthy paths in a use-of-force network. [`code/mcmcOff10.R`](code/mcmcOff10.R)

All three simulations include versions with large and small datasets. The large datasets are intended to avoid simulation error and simply show that the algorithm works as expected. The small datasets are intended to demonstrate the model's performance on more data with a more realistic incidents-per-officer rate.

This repository also contains data and code for analyzing seven years of use-of-force data from Seattle PD: [`code/mcmcSPD.R`](code/mcmcSPD.R), [`data/dataSPD.RData`](data/dataSPD.RData). As with the simulations, no other setup is needed to run this script.

The output from the three simulation scripts and the Seattle PD script are the MCMC draws from the posterior distribution.

The script [`code/postMCMCanalysis.R`](code/postMCMCanalysis.R) contains code to analyze the posterior samples from the three simulations and the Seattle PD data. The script labels blocks of code that connect to specific sections, tables, and figures in the article.

# Reproducing results in the appendix

## Appendix E: Simulations of dependence

- [`code/mcmcOff3depend.R`](code/mcmcOff3depend.R) contains simulation code for three experiments
    - Appendix E.1 Three officers where Officer 1 escalates/de-escalates others
    - Appendix E.2 Thirty officers where Officer 1 de-escalates peers and Officer 30 escalates them
    - Appendix E.3 A randomly selected "first" officer escalates or de-escalates peers

## Appendix F: Placebo test of outlier detection

- [`code/mcmcSPDpermute.R`](code/mcmcSPDpermute.R) R script for Appendix F, placebo test that runs the conditional ordinal stereotype model on permuted versions of the Seattle data. Running this script requires several days to complete, about 8 days on a 40-core machine I used and about 15 days on a 16-core machine.

- [`code/postMCMCanalysisSPDpermute.R`](code/postMCMCanalysisSPDpermute.R) R script for processing the results of the MCMC draws produced from [`code/mcmcSPDpermute.R`](code/mcmcSPDpermute.R)


# Code for verifying calculations

- [`code/testCLcalculations.R`](code/testCLcalculations.R) A range of tests to ensure that the conditional likelihood calculation is correct. This script is not used to produce any results for the article. It shows a range of tests that were used to ensure the correctness of the calculations

    - verify calculations are identical if using 
        1. Complete enumeration (using `arrangements::permutations`)
        2. C++ implementation of no-repeat Heap's algorithm (`denomHeaps_R`)
        3. Lin, Wang, Hong (2023) computation of the Poisson-Multinomial (`dpmd`)
        4. C++ implementation of the multi-dimensional discrete Fourier transform (`logDenomDFT`)

    - timing test comparing efficiency of multi-dimensional fast Fourier transform or the multi-dimension discrete Fourier transform
    
    - timing test comparing efficiency of the C++ implementation of no-repeat Heap's algorithm and the C++ implementation of the multi-dimensional discrete Fourier transform. Determines how many officers at an incident before the DFT becomes more efficient than Heap's algorithm
    
    - tests of the log conditional likelihood calculation in serial and in parallel


# R package dependencies

- For implementation of the conditional ordinal stereotype model
    - Rcpp 1.0.14
    - RcppParallel 5.1.10
    
- For testing analytic calculations
    - arrangements 1.1.9
    - PoissonMultinomial 1.1

- For tables and figures
    - ggplot2 3.5.2
    - xtable 1.8-4
    - viridis 0.6.5
    - RColorBrewer 1.1-3
    - ggbeeswarm 0.7.2
    - igraph 2.1.4

- For processing MCMC draws
    - doParallel 1.0.17
    - parallel 4.5.1
    - foreach 1.5.2
    - future 1.58.0
    - progressr 0.15.1
    - doFuture 1.1.1

- For general data manipulation
    - dplyr 1.1.4
    - tidyr 1.3.1
    - purrr 1.0.4


