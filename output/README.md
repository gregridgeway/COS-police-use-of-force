OUTPUT - Conditional Ordinal Stereotype Model to Estimate Police Officers' Propensity to Escalate Force
================

This folder contains all the output from the simulations and the analysis of the Seattle PD data.

# Results for main text of article

- Simulation #1, Section 4.1, three well-connected officers
    - [`output/mcmcSampOff3.RData`](mcmcSampOff3.RData): MCMC samples for the large dataset with 115 incidents per officer
    - [`output/mcmcSampOff3small.RData`](mcmcSampOff3small.RData): MCMC samples for the small dataset with 10 incidents per officer
    
- Simulation #2, Section 4.2, six officers in a disconnected network
    - [`output/mcmcSampOff6.RData`](mcmcSampOff6.RData): MCMC samples for the large dataset with 124 incidents per officer
    - [`output/mcmcSampOff6small.RData`](mcmcSampOff6small.RData): MCMC samples for the small dataset with 10 incidents per officer
    
- Simulation #3, Section 4.3, ten officers in a long network
    - [`output/mcmcSampOff10.RData`](mcmcSampOff10.RData): MCMC samples for the large dataset with 144 incidents per officer
    - [`output/mcmcSampOff10small.RData`](mcmcSampOff10small.RData): MCMC samples for the small dataset with 19 incidents per officer
    
- Seattle PD analysis, Section 5. Posterior draws from four independent MCMC chains, each with 1 million draws, 20,000 burn-in iterations, and thinning to every 100 draws (posted to github LFS)
    - [`output/mcmcSampSPDchain1.RData`](mcmcSampSPDchain1.RData)
    - [`output/mcmcSampSPDchain2.RData`](mcmcSampSPDchain2.RData)
    - [`output/mcmcSampSPDchain3.RData`](mcmcSampSPDchain3.RData)
    - [`output/mcmcSampSPDchain4.RData`](mcmcSampSPDchain4.RData)
    
- Post-processing results from the Seattle PD analysis are in [`output/resultsSPD.RData`](resultsSPD.RData). Processing the MCMC samples can take a long time, so this file is provided to avoid the need to re-run the post-processing code. The script [`code/postMCMCanalysis.R`](../code/postMCMCanalysis.R) contains the code that generated this file.

- [`code/postMCMCanalysis.R`](../code/postMCMCanalysis.R) also generated the three pdf figures here
    - [`output/pairs3off1.pdf`](pairs3off1.pdf): Figure C1
    - [`output/s2s3.pdf`](s2s3.pdf): Figure C2
    - [`output/sDiffs.pdf`](sDiffs.pdf): Figure C3
    
# Supplementary appendix output

## Appendix E: Simulations studying the effect of violations of conditional independence

- Stored results for the three simulations
    - [`output/mcmcOff3depend.RData`](mcmcOff3depend.RData): Simulation #1, three officers where Officer 1 escalates/de-escalates others
    - [`output/mcmcOff30depend.RData`](mcmcOff30depend.RData): Simulation #2, thirty officers where Officer 1 de-escalates peers and Officer 30 escalates them
    - [`output/mcmcOff3firstOff.RData`](mcmcOff3firstOff.RData): Simulation #3, a randomly selected "first" officer among three escalates or de-escalates peers

## Appendix F: Placebo test of outlier detection

- MCMC samples from the placebo test The script [`code/mcmcSPDpermute.R`](../code/mcmcSPDpermute.R) runs the conditional ordinal stereotype model on permuted versions of the Seattle PD data. The MCMC samples are stored in the following files, each containing 1-2 million draws. I originally ran these on two separate machines, chain 1 and chain 4 generating 2M draws, and chains 2, 3, and 5 generating 1M draws each  (posted to github LFS)
    - [`output/mcmcSampSPDPermutechain1.RData`](mcmcSampSPDPermutechain1.RData) (2M draws)
    - [`output/mcmcSampSPDPermutechain2.RData`](mcmcSampSPDPermutechain2.RData) (1M draws)
    - [`output/mcmcSampSPDPermutechain3.RData`](mcmcSampSPDPermutechain3.RData) (1M draws)
    - [`output/mcmcSampSPDPermutechain4.RData`](mcmcSampSPDPermutechain4.RData) (2M draws)
    - [`output/mcmcSampSPDPermutechain5.RData`](mcmcSampSPDPermutechain5.RData) (1M draws)

- Analysis of the MCMC draws for the placebo tests
    - [`output/resultsSPDpermute.RData`](resultsSPDpermute.RData): This file contains the results of the post-processing of the MCMC samples from the placebo test. The script [`code/postMCMCanalysisSPDpermute.R`](../code/postMCMCanalysisSPDpermute.R) contains the code that generated this file
