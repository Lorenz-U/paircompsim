# Generate data for simulation studies
######################################

# A: unshaved
# B: mechanic
# C: digital

# n_raters and n_questions should not be changed!
# Otherwise, the simulations will not work anymore.

# Parameters:
# n:                    number of subjects per question and rater
# n_rater:              number of raters
# n_questions:          number of questions
# d_AB, d_AC, d_CB:     assumed effects (for A, B, C, see above)
# sigma_AB, sigma_AB,
# sigma_CB:             sigma per effect (see model)
# tau_AB, tau_AC,
# tau_CB:               tau for heterogeneity between raters

# Output:
# y: data matrix of dimensions (n, n_questions, n_raters)

# truncated normal distribution (with symmetric truncation)
# Idea: draw a normally distributed value, round it, and truncate it
# if it is outside the range of mu +/- tr.
# For the simulation study, this is only meaningful if the type I error
# is assessed, with mu=0 and tr=3. In the real data example,
# the values were all in between -3 and 3.

rnorm.trunc <- function(n, mu, s, tr) {
  y <- round(rnorm(n, mu, s))
  y.w <- which((y < mu - tr) | (y > mu + tr))
  while(length(y.w) != 0) {
    y[y.w] <- round(rnorm(length(y.w), mu, s))
    y.w <- which((y < mu - tr) | (y > mu + tr))
  }
  y
}

sim.dat <- function(n = 59, n_raters = 6, n_questions = 3,
                    d_AB = 0, d_AC = 0, d_CB = 0,
                    sigma_AB = 1, sigma_AC = 1, sigma_CB = 1,
                    tau_AB = 0, tau_AC = 0, tau_CB = 0) {
  d_vec <- c(d_AB, d_AC, d_CB)
  sigma_vec <- c(sigma_AB, sigma_AC, sigma_CB)
  tau_vec <- c(tau_AB, tau_AC, tau_CB)
  y <- matrix(NA, n, n_raters*n_questions)

  if (all(tau_vec == 0)) { # if no heterogeneity:
    for (k in 1:n_raters) {
      for (j in 1:n_questions) {
        y[, n_questions*(k-1) + j] <- rnorm.trunc(n, d_vec[j], sigma_vec[j], 3)
      }
    }
  } else { # if heterogeneity:
    for (k in 1:n_raters) {
      for (j in 1:n_questions) {
        d_kj <- rnorm(1, d_vec[j], tau_vec[j])
        y[, n_questions*(k-1) + j] <- rnorm.trunc(n, d_kj, sigma_vec[j], 3)
      }
    }
  }
  return(y)
}


# sim.dat <- function(n = 59, n_raters = 6, n_questions = 3,
#                     d_AB = 0, d_AC = 0, d_CB = 0,
#                     sigma_AB = 1, sigma_AC = 1, sigma_CB = 1,
#                     tau_AB = 0, tau_AC = 0, tau_CB = 0) {
#   d_vec <- c(d_AB, d_AC, d_CB)
#   sigma_vec <- c(sigma_AB, sigma_AC, sigma_CB)
#   tau_vec <- c(tau_AB, tau_AC, tau_CB)
#   y <- matrix(NA, n, n_raters*n_questions)
#   
#   if (all(tau_vec == 0)) { # if no heterogeneity:
#     for (k in 1:n_raters) {
#       for (j in 1:n_questions) {
#         y[, n_questions*(k-1) + j] <- sample(-3:3, n, replace = T)
#       }
#     }
#   } else { # if heterogeneity:
#     for (k in 1:n_raters) {
#       for (j in 1:n_questions) {
#         d_kj <- rnorm(1, d_vec[j], tau_vec[j])
#         y[, n_questions*(k-1) + j] <- rnorm.trunc(n, d_kj, sigma_vec[j], 3)
#       }
#     }
#   }
#   return(y)
# }

# Test run:
# y <- sim.dat()
# y <- sim.dat(tau_AB = 2)
