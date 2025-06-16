# Estimating the impact of myrtle rust on individual species

This subfolder contains code to estimate the impact myrtle rust on individual species, given ordinal ratings/scores (e.g., "Relatively Tolerant" or "Extremely Susceptible").

Briefly, it uses Approximate Bayesian Computation to estimate the shape of the damage/impact density curve, and estimate the location of the thresholds that separate "Relatively Tolerant" from "Moderately Susceptible" and so forth (see the paper for more details).

## Steps to run analysis

Run: **bayes_ABC.py** for the samples. (Recommended: running on a computing cluster)
In total: 20,000 accepted samples were generated.

Run: **bayes_ABC_plotting.py** for plotting, including visual checking of convergence.

Run: **bayes_ABC_generate_MR_impacts.py** to produce csvs of damage values for different plant species.

Final output: see the folder **/MR_samples** for 1000 different species-impact parameter sets.