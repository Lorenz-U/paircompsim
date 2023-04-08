# Function to implement models
##############################

# This functions only works with 3 questions!


# Function used to calculate p-values
# ===================================

pv <- function(m, sd, n = 100000) {
  sim <- rnorm(n, 0, sd)
  res <- length(which(abs(sim) > abs(m))) / n
  res
}


# Fixed model for normal data
# ===========================

# Parameters:
# y:            data matrix
# n.burnin:     number of burn-in iterations
# n.iter:       number of MCMC iterations
# n.thin:       thinning rate

# Output:
# p_AB, p_AC, p_CB: p-values for comparisons of AB, AC, and BC, respectively
# summary.fit: summary of results based on MCMC iterations (i.e., posterior distributions)

fixed.normal <- function(y, n.burnin = 5000, n.iter = 25000, n.thin = 5) {

  options(warn = -1)
  suppressMessages(require(rjags))
  
  a.c <- matrix(c(1, 1, 3, 2, 3, 2), ncol = 2)
  
  npat <- dim(y)[1]
  nfrag <- dim(y)[2]
  nrat <- dim(y)[3]
  bugs.seed <- 12345
  
  inits <- list(
    list(d = c(NA, 0.1, -0.1), .RNG.name = "base::Mersenne-Twister", .RNG.seed = bugs.seed),
    list(d = c(NA, -0.05, 0.05), .RNG.name = "base::Mersenne-Twister", .RNG.seed = bugs.seed)
  )

  fitjags <- jags.model("model_normal_fixed.txt",
                        list("y" = y, "n.p" = npat,
                             "n.c" = nfrag, "n.r" = nrat,
                             "a.c" = a.c), inits = inits,
                        n.chains = 2, quiet = T, n.adapt = 1000)
  
  update(fitjags, n.burnin, progress.bar = "none")
  fit.normal.fixed <- coda.samples(fitjags, c("d", "dd", "best", "sd"), n.iter, n.thin, progress.bar = "none")

  # Results:
  # d[1] is set to 0. It describes the effect of unshaved.
  # d[2] describes the effect of mechanically shaved. If d[2] > 0, then it is better than unshaved.
  # d[3] describes the effect of digitally shaved. If d[3] > 0, then it is better than unshaved.
  # The effects describe mean differences.
  # dd describes the effect d[2] - d[3].

  sf <- summary(fit.normal.fixed)[[1]]

  p1 <- pv(sf[5, 1], sf[5, 2]) # for d[2]
  p2 <- pv(sf[6, 1], sf[6, 2]) # for d[3]
  p3 <- pv(sf[7, 1], sf[7, 2]) # for dd
  
  return(list("p_AB" = p1, "p_AC" = p2, "p_CB" = p3,
              "summary.fit" = summary(fit.normal.fixed)))
}


# Random effects model for normal data
# ====================================

# Parameters:
# y:            data matrix
# n.burnin:     number of burn-in iterations
# n.iter:       number of MCMC iterations
# n.thin:       thinning rate

# Output:
# p_AB, p_AC, p_CB: p-values for comparisons of AB, AC, and BC, respectively
# summary.fit: summary of results based on MCMC iterations (i.e., posterior distributions)

random.normal <- function(y, n.burnin = 5000, n.iter = 25000, n.thin = 5) {

  options(warn = -1)
  suppressMessages(require(rjags))
  
  a.c <- matrix(c(1, 1, 3, 2, 3, 2), ncol = 2)
  
  npat <- dim(y)[1]
  nfrag <- dim(y)[2]
  nrat <- dim(y)[3]
  bugs.seed <- 12345
  
  inits <- list(
    list(d = c(NA, 0.1, -0.1), prec = c(0.1, 0.5, 1),
         tau = c(0.2, 0.4, 0.3), .RNG.name = "base::Mersenne-Twister", .RNG.seed = bugs.seed),
    list(d = c(NA, -0.05, 0.05), prec = c(1, 0.1, 0.5),
         tau = c(0.4, 0.3, 0.2), .RNG.name = "base::Mersenne-Twister", .RNG.seed = bugs.seed)
  )
  
  fitjags <- jags.model("model_normal_random.txt",
                        list("y" = y, "n.p" = npat,
                             "n.c" = nfrag, "n.r" = nrat,
                             "a.c" = a.c), inits = inits,
                        n.chains = 2, quiet = T, n.adapt = 1000)
  
  update(fitjags, n.burnin, progress.bar = "none")
  fit.normal.random <- coda.samples(fitjags, c("d", "dd", "best", "sd",
                                               "sdtau"),
                                    n.iter, n.thin, progress.bar = "none")
  
  # d is interpreted in the same way as in the fixed model.
  # tau is the precision of delta
  # sdtau is the standard error per comparison between raters:
  # sdtau[1] describes the standard error of unshaved vs. mechanically shaved
  # sdtau[2] describes the standard error of unshaved vs. digitally shaved
  # sdtau[3] describes the standard error of digitally vs. mechanically shaved
  # sd is the standard error within raters; it is ordered
  # in the same way as sdtau.

  sf <- summary(fit.normal.random)[[1]]

  p1 <- pv(sf[5, 1], sf[5, 2], 12345)
  p2 <- pv(sf[6, 1], sf[6, 2], 12345)
  p3 <- pv(sf[7, 1], sf[7, 2], 12345)
  
  return(list("p_AB" = p1, "p_AC" = p2, "p_CB" = p3,
              "summary.fit" = summary(fit.normal.random)))
}