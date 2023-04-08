# Simulation study
##################

# Load all R files:
source("do_simulation.R")
source("generate_data.R")
source("models.R")
source("evaluate.R")

# Scenarios:
# - no heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY = 0)
# - heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY = 0.05)

# no heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY = 0)
# =============================================================================

sim.study1 <- do.sim(n.sim = 10000,
                     model = c("fixed.ordinal", "random.ordinal",
                               "wilcox", "wilcox.clust"),
                     seed = 76057, alpha = 0.05,
                     sigma_AB = 2, sigma_AC = 2, sigma_CB = 2,
                     n.burnin = 1000, n.iter = 10000, n.thin = 2,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study1, k = 5)


# heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY = 0.05)
# =============================================================================

sim.study2 <- do.sim(n.sim = 10000,
                     tau_AB = sqrt(0.05), tau_AC = sqrt(0.05), tau_CB = sqrt(0.05),
                     model = c("fixed.ordinal", "random.ordinal",
                               "wilcox", "wilcox.clust"),
                     seed = 85206, alpha = 0.05,
                     sigma_AB = 2, sigma_AC = 2, sigma_CB = 2,
                     n.burnin = 1000, n.iter = 10000, n.thin = 2,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study2, k = 5)