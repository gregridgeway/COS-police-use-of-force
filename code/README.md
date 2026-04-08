CODE - Conditional Ordinal Stereotype Model to Estimate Police Officers' Propensity to Escalate Force
================

# Code for main text of article

- [`conditional_stereotype.cpp`](conditional_stereotype.cpp) contains the C++ code that parallelizes the conditional likelihood calculation by incident

- [`mcmcConditionalStereotype.R`](mcmcConditionalStereotype.R) contains the Metropolis algorithm that samples from the posterior distribution of the conditional stereotype model parameters. Any script calling this script should first call `RcppParallel::setThreadOptions(numThreads = 8)`, editing `numThreads` to be appropriate for the computing environment.

- [`mcmcOff3.R`](mcmcOff3.R) Simulation #1

- [`mcmcOff3+3.R`](mcmcOff3+3.R) Simulation #2

- [`mcmcOff10.R`](mcmcOff10.R) Simulation #3

- [`mcmcSPD.R`](mcmcSPD.R) Analysis of use-of-force data from Seattle PD

- [`postMCMCanalysis.R`](postMCMCanalysis.R) Analysis of the posterior samples from the three simulations and the Seattle PD data. Labels blocks of code that connect to specific sections, tables, and figures in the article


# Supplemental results code

## Appendix E: Simulations of dependence

- [`mcmcOff3depend.R`](mcmcOff3depend.R) contains simulation code for three experiments
    - Appendix E.1 Three officers where Officer 1 escalates/de-escalates others
    - Appendix E.2 Thirty officers where Officer 1 de-escalates peers and Officer 30 escalates them
    - Appendix E.3 A randomly selected "first" officer escalates or de-escalates peers

## Appendix F: Placebo test of outlier detection

- [`mcmcSPDpermute.R`](mcmcSPDpermute.R) R script for Appendix F, placebo test that runs the conditional ordinal stereotype model on permuted versions of the Seattle data

- [`postMCMCanalysisSPDpermute.R`](postMCMCanalysisSPDpermute.R) R script for processing the results of the MCMC draws produced from [`mcmcSPDpermute.R`](mcmcSPDpermute.R)


# Code for verifying calculations

- [`testCLcalculations.R`](testCLcalculations.R) A range of tests to ensure that the conditional likelihood calculation is correct. This script is not used to produce any results for the article. It shows a range of tests that were used to ensure the correctness of the calculations

    - verify calculations are identical if using 
        1. Complete enumeration (using arrangements::permutations)
        2. C++ implementation of no-repeat Heap's algorithm
        3. Lin, Wang, Hong (2023) computation of the Poisson-Multinomial (dpmd)
        4. C++ implementation of the multi-dimensional discrete Fourier transform (logDenomDFT)
        5. C++ implementation of a dynamic program (logDenomDP), added March 2026 post-publication

    - timing test comparing efficiency of multi-dimensional fast Fourier transform or the multi-dimension discrete Fourier transform
    
    - timing test comparing efficiency of the C++ implementation of no-repeat Heap's algorithm and the C++ implementation of the multi-dimensional discrete Fourier transform. Determines how many officers at an incident before the DFT becomes more efficient than Heap's algorithm. Added timing tests of dynamic program in March 2026.
    
    - empirically assessed computational complexity of dynamic program. Added C++ function logDenomDP_count that collects statistics on the computation (states visited, number of multiplications, number of exp() calls, etc.). Worst case computational complexity is O(m*J*K), where K = prod(k[j]+1). K behaves like (m/J)^J (worst case), so complexity is roughly O(m^(J+1)/J^(J-1)). Empirically with J=4, I observe O(m^1.35) (April 2026).
    
    - tests of the log conditional likelihood calculation in serial and in parallel
