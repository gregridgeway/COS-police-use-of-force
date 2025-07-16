`dataSPD.RData` contains the Seattle PD use-of-force data already formatted for analysis, a single data frame `d`. The data frame contains the following columns:

- `id`: use-of-force incident identifier, sequential integers from 1 to 4,821

- `idOff`: officer identifier, integers from 1 to 1,503. `idOff` is a scrambled identifier and does not link to any Seattle PD officer information

- `y`: ordinal use-of-force level, the most servious level of force that the officer used at this incident. `y` contains integers from 1 to 4, where 1 is the lowest severity and 4 is the highest severity

Incidents with no variation in `y` (e.g. single officer incidents, incidents where no officer uses force) have been removed from the data since they have no information for the parameters of interest through the conditional likelihood.
