# Simulation studies for pairwise comparisons

In this repository, R and JAGS code can be found to implement simulation studies as presented in the manuscript "A Bayesian approach to  clinical trial designs in dermatology with multiple simultaneous treatments per subject and multiple raters" by Lorenz Uhlmann, Christian Stock,
Marc Vandemeulebroecke, Christine Fink, and Meinhard Kieser (currently under review for the journal "Contemporary Clinical Trials").

There are two folders:
- "normal_endpoint" which contains R and JAGS code to evaluate models with a normally distributed endpoint
- "ordinal_endpoint" which contains R and JAGS code to evaluate models with an ordinal endpoint

To run a simulation study, the file "sim_study.R" can be used. In this file, all other files are loaded via "source()".