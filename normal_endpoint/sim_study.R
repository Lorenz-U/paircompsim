# Simulation study
##################

# Load all R files:
source("do_simulation.R")
source("generate_data.R")
source("models.R")
source("evaluate.R")

# Scenarios:
# - no heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY^2 = 0)
# - heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY^2 = 0.05)
# - heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY^2 = 0.2)
# - smaller sample size, no heterogeneity, no inconsistency, no effect
# - smaller sample size, heterogeneity, no inconsistency, no effect
# - no heterogeneity, consistency, some effect
# - heterogeneity, consistency, some effect
# - no heterogeneity, inconsistency, small direct effect, higher indirect effect
# - heterogeneity, inconsistency, small direct effect, higher indirect effect


# no heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY = 0)
# =============================================================================

sim.study1 <- do.sim(n.sim = 10000,
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 12345, alpha = 0.05,
                     n.cores = 4)
# evaluate results:
eval.sim(sim.study1, k = 5)


# heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY = 0.05)
# =============================================================================

sim.study2 <- do.sim(n.sim = 10000,
                     tau_AB = sqrt(0.05), tau_AC = sqrt(0.05), tau_CB = sqrt(0.05),
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 12345, alpha = 0.05,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study2, k = 5)


# heterogeneity, consistency, no effect (d_XY = 0, simga_XY = 1, tau_XY^2 = 0.2)
# ============================================================================

sim.study3 <- do.sim(n.sim = 10000,
                     tau_AB = sqrt(0.2), tau_AC = sqrt(0.2), tau_CB = sqrt(0.2),
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 19560, alpha = 0.05,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study3, k = 5)


# smaller sample size, no heterogeneity, no inconsistency, no effect
# ==================================================================

sim.study4 <- do.sim(n.sim = 10000,
                     n = 30, n_raters = 3,
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 12345, alpha = 0.05,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study4, k = 5)


# smaller sample size, heterogeneity, no inconsistency, no effect
# ===============================================================

sim.study5 <- do.sim(n.sim = 10000,
                     n = 30, n_raters = 3,
                     tau_AB = sqrt(0.05), tau_AC = sqrt(0.05), tau_CB = sqrt(0.05),
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 12345, alpha = 0.05,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study5, k = 5)


# no heterogeneity, consistency, some effect
# ==========================================

sim.study6 <- do.sim(n.sim = 10000,
                     d_AC = 0.122, d_AB = 0.061, d_CB = -0.061,
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 12345, alpha = 0.05,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study6, d_true = 0.122, k = 5)


# heterogeneity, consistency, some effect
# =======================================

sim.study7 <- do.sim(n.sim = 10000,
                     d_AC = 0.27, d_AB = 0.135, d_CB = -0.135,
                     tau_AB = sqrt(0.05), tau_AC = sqrt(0.05), tau_CB = sqrt(0.05),
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 12345, alpha = 0.05,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study7, d_true = 0.27, k = 5)



# no heterogeneity, inconsistency, small direct effect, higher indirect effect
# ============================================================================

# overall effect of 0.1225.
sim.study8 <- do.sim(n.sim = 10000,
                     d_AC = 0.12, d_AB = 0.07, d_CB = -0.07,
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 12345, alpha = 0.05,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study8, d_true = 0.1225, k = 5)


# heterogeneity, inconsistency, small direct effect, higher indirect effect
# =========================================================================

# overall effect of 0.273
sim.study9 <- do.sim(n.sim = 10000,
                     d_AC = 0.25, d_AB = 0.17, d_CB = -0.17,
                     tau_AB = sqrt(0.05), tau_AC = sqrt(0.05), tau_CB = sqrt(0.05),
                     model = c("fixed.normal", "random.normal",
                               "t.test", "fr.mixed"),
                     seed = 12345, alpha = 0.05,
                     n.cores = 4)

# evaluate results:
eval.sim(sim.study9, d_true = 0.28, k = 5)
