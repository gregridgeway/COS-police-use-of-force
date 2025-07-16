
This folder contains all the output from the simulations and the analysis of the Seattle PD data.

- Simulation #1, Section 4.1, three well-connected officers
    - [`output/mcmcSampOff3.RData`](output/mcmcSampOff3.RData): MCMC samples for the large dataset with 115 incidents per officer
    - [`output/mcmcSampOff3small.RData`](output/mcmcSampOff3small.RData): MCMC samples for the small dataset with 10 incidents per officer
    
- Simulation #2, Section 4.2, six officers in a disconnected network
    - [`output/mcmcSampOff6.RData`](output/mcmcSampOff6.RData): MCMC samples for the large dataset with 124 incidents per officer
    - [`output/mcmcSampOff6small.RData`](output/mcmcSampOff6small.RData): MCMC samples for the small dataset with 10 incidents per officer
    
- Simulation #3, Section 4.3, ten officers in a long network
    - [`output/mcmcSampOff10.RData`](output/mcmcSampOff10.RData): MCMC samples for the large dataset with 144 incidents per officer
    - [`output/mcmcSampOff10small.RData`](output/mcmcSampOff10small.RData): MCMC samples for the small dataset with 19 incidents per officer
    
- Seattle PD analysis, Section 5. Posterior draws from four independent MCMC chains, each with 1 million draws, 20,000 burn-in iterations, and thinning to every 100 draws
    - [`output/mcmcSampSPDchain1.RData`](output/mcmcSampSPDchain1.RData)
    - [`output/mcmcSampSPDchain2.RData`](output/mcmcSampSPDchain2.RData)
    - [`output/mcmcSampSPDchain3.RData`](output/mcmcSampSPDchain3.RData)
    - [`output/mcmcSampSPDchain4.RData`](output/mcmcSampSPDchain4.RData)
    
- Post-processing results from the Seattle PD analysis are in [`output/resultsSPD.RData`](output/resultsSPD.RData). Processing the MCMC samples can take a long time, so this file is provided to avoid the need to re-run the post-processing code. The script [`code/post mcmc analysis.R`](../code/post mcmc analysis.R) contains the code that generated this file.

- [`code/post mcmc analysis.R`](../code/post mcmc analysis.R) also generated the
