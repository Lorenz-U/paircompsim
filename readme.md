# Simulation studies for pairwise comparisons

In this repository, R and JAGS code can be found to implement simulation studies as presented in the manuscript "A Bayesian approach to  clinical trial designs in dermatology with multiple simultaneous treatments per subject and multiple raters" by Lorenz Uhlmann, Christian Stock,
Marc Vandemeulebroecke, Christine Mueller-Christmann, and Meinhard Kieser.

To cite the article, you can use:

Uhlmann L, Stock C, Vandemeulebroecke M, Mueller-Christmann C, Kieser M. A Bayesian approach to clinical trial designs in dermatology with multiple simultaneous treatments per subject and multiple raters. Contemp Clin Trials. 2023 May 22;131:107233. doi: 10.1016/j.cct.2023.107233. Epub ahead of print. PMID: 37225121.

There are two folders:
- "normal_endpoint" which contains R and JAGS code to evaluate models with a normally distributed endpoint
- "ordinal_endpoint" which contains R and JAGS code to evaluate models with an ordinal endpoint

To run a simulation study, the file "sim_study.R" can be used. In this file, all other files are loaded via "source()".