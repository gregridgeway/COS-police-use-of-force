Conditional Ordinal Stereotype Model to Estimate Police Officers' Propensity to Escalate Force
================

# Overview
This repository contains the code and data associated with "A Conditional Ordinal Stereotype Model to Estimate Police Officers' Propensity to Escalate Force." The article introduced a conditional ordinal stereotype model for estimating a police officer's latent propensity to escalate to higher-severity force options during an encounter. Propensity to escalate is the likelihood of selecting a more serious force category than peer officers confronting the same circumstances. The associated conditional likelihood depends solely on data from the times and places where multiple officers use various forms of force.

Code for the conditional stereotype model is a blend of R and C++. The C++ code ([`code/conditional_stereotype.cpp`](code/conditional_stereotype.cpp)) contains the C++ code that parallelizes the conditional likelihood calculation by incident. For incidents with two officers, there is a simple, direct calculation. For incidents with three to seven officers, there is a no-repeats version of Heap's algorithm O(m!). For incidents with eight or more officers, there is a multi-dimensional discrete Fourier transform O(m^J). This arrangement results in the most efficient conditional likelihood calculation. The R script ([`code/mcmcConditionalStereotype.R`](code/mcmcConditionalStereotype.R)) contains the Metropolis algorithm that samples from the posterior distribution of the conditional stereotype model parameters.

The [`code`](code/) folder contains three scripts associated with the three simulated examples in the article.
- Simulation #1: Demonstrates proof-of-concept in a well-connected three officer network. [`code/mcmcOff3.R`](code/mcmcOff3.R)
- Simulation #2: Shows the challenge of analyzing disconnected use-of-force networks and that the conditional variance helps to partition the network into peers and non-peers. [`code/mcmcOff3+3.R`](code/mcmcOff3+3.R)
- Simulation #3: Shows that information can flow through lengthy paths in a use-of-force network. [`code/mcmcOff10.R`](code/mcmcOff10.R)

All three simulations include versions with large and small datasets. The large datasets are intended to avoid simulation error and simply show that the algorithm works as expected. The small datasets are intended to demonstrate the model's performance on more data with a more realistic incidents-per-officer rate.
This repository also contains data and code for analyzing seven years of use-of-force data from Seattle PD: [`code/mcmcSPD.R`](code/mcmcSPD.R), [`data/dataSPD.RData`](data/dataSPD.RData)

# R package dependencies
`dplyr`
`Rcpp`
`RcppParallel`
`doParallel`
`future`
`doFuture`
`progressr`
`igraph`
`xtable`
`ggplot2`
`ggbeeswarm`



**Importantly, the authors should provide an overview of how to carry
out the analyses presented in their manuscript in the `README.md` of their
repository, replacing the content in this file.** This overview would
generally refer to scripts/code files that execute the analyses and are
placed either in the main directory or the `/code` subdirectory. The
*Workflow* section of the ACC form should refer to this README.md as
containing the instructions for how to reproduce the analyses.

## Reproducibility materials file structure

The suggested components are as follows. Directories in the submission may have subdirectories to
further organize the materials.

1.  A `README.md` file - This file gives a short description of the
    paper and an overview of how to carry out the analyses presented in their manuscript.
2.  A `manuscript` directory - This directory will generally hold the source files
    (often LaTeX or Rmd) for the manuscript and any files directly related to the
    generation of the manuscript, including figure files.
3.  A `data` directory - This directory will generally hold the real data files 
    (or facsimile versions of them in place of confidential data) and simulated data files.
    See `data/README.md` for more details. 
4.  A `code` directory - This directory will generally hold 
    source code files that contain the core code to implement the method and various utility/auxiliary functions.
5.  An `output` directory - This directory will generally hold objects derived
    from computations, including results of simulations or real data analyses. See `output/README.md` for more details.

