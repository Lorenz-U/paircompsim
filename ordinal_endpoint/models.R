# Function to implement models
##############################

# This functions only works with 3 questions!
# IMPORTANT: d_BC is actually d_CB which means that you need to take -d_BC to get d_CB.
# Now it is updated and it says d_CB.

# source("generate_data.R")

# Functions used to calculate p-values
# ====================================

# f.r <- function(x, k) format(round(x, k), nsmall = k)
# format.r <- function(x, k = 3, cl.z = F) {
#   val <- f.r(x, k)
#   if (cl.z) {
#     val[which(val == f.r(0, k))] <- paste("<", 0.1^k, sep = "")
#   }
#   val
# }
# pv <- function(m, sd, s, n = 100000, digits = 3) {
#   set.seed(s)
#   sim <- rnorm(n, 0, sd)
#   res <- length(which(abs(sim) > abs(m))) / n
#   # format.r(res, digits, T)
#   res
# }

pv <- function(m, sd, n = 100000) {
  sim <- rnorm(n, 0, sd)
  res <- length(which(abs(sim) > abs(m))) / n
  res
}

# Transform the data into useful format
# =====================================

# We model P(y<=i|x) = f(mu_i - d*comp)
# fun.trans <- function(x, comp, rater) {
#   x.u <- sort(unique(x))
#   r <- matrix(0, nrow = length(x), ncol = length(x.u) - 1)
#   for (i in 1:length(x)) {
#     if (x[i] < x.u[length(x.u)]) {
#       w.x <- (which(x[i] == x.u)):(length(x.u) - 1)
#       r[i, w.x] <- 1
#     }
#   }
#   r.a <- n.a <- array(0, dim = c(length(unique(rater)), ncol(r), length(unique(comp))))
#   for (i in 1:ncol(r)) {
#     r.a[, i, ] <- aggregate(r[, i], by = list(rater, comp), FUN = sum)$x
#     n.a[, i, ] <- aggregate(r[, i], by = list(rater, comp), FUN = length)$x
#   }
#   return(list("r.a" = r.a, "n.a" = n.a))
# }

# For dmulti:
fun.trans <- function(x, comp, rater) {
  
  r <- sort(unique(as.vector(x)))
  x.dat <- as.data.frame(x)
  for (j in 1:ncol(x)) {
    x.dat[, j] <- factor(x.dat[, j], levels = r)
  }
  r.a <- array(NA, dim = c(length(unique(rater)), length(r), length(unique(comp))))
  
  for (i in 1:dim(r.a)[1]) {
    for (j in 1:dim(r.a)[3]) {
      r.a[i, , j] <- table(x.dat[, 3*(j-1) + i])
    }
  }
  n.ov <- nrow(x.dat)
  return(list("r.a" = r.a, "n.ov" = n.ov))
}



# Fixed model for ordinal data
# ============================

# Parameters:
# x:            data matrix
# n.burnin:     number of burn-in iterations
# n.iter:       number of MCMC iterations
# n.thin:       thinning rate

# Output:
# p_AB, p_AC, p_CB: p-values for comparisons of AB, AC, and BC, respectively
# summary.fit: summary of results based on MCMC iterations (i.e., posterior distributions)

fixed.ordinal <- function(x, n.burnin = 5000, n.iter = 25000, n.thin = 5) {

  options(warn = -1)
  suppressMessages(require(rjags))
  
  # Create one vector of all results:
  y <- c()
  for (i in 1:ncol(x)) {
    y <- c(y, x[, i])
  }
  
  # Add information of comparison and rater as separate columns:
  comp <- rep(1:(ncol(x)/3), each = nrow(x) * 3) # three columns are put below each other per rater
  rater <- rep(rep(1:3, each = nrow(x)), (ncol(x)/3)) # each question is asked in one row
  
  r <- fun.trans(x, comp, rater) # x or y?! first function: y; second function: x!
  
  # structure of r:
  # first dimension: questions (the object name "rater" above is misleading!)
  # second dimension: ratings - 1:
  # Sum of ratings:
  # If a rater says 3 (favors second option), then no count
  # If a rater says 2 (favors second option), then last element is count + 1.
  # If a rater says 1 (favors second option), then the last and the one before is count + 1.
  # ...
  # third dimension: raters (the object name "comp" above is misleading!)
  # Example:
  # r$r.a[1, 6, 2] describes the number of ratings of rater 2 in question 1 (unshaved vs.
  # mechanically shaved) that were 2 or lower (because we are looking at the 6th element of
  # the vector (-3, -2, -1, 0, 1, 2, 3).
  
  # 1=unshaved; 2=mechanically shaved; 3=digitally shaved:
  a.c <- matrix(c(1, 1, 3, 2, 3, 2), ncol = 2)
  
  require(coda)
  
  # Fixed effects model:
  inits <- list(
    list (d = c(NA, 0.1, -0.1),
          mucat0 = seq(-0.3, 0.2, by = 0.1),
          .RNG.name = "base::Mersenne-Twister", .RNG.seed = 12345),
    list (d = c(NA, -0.05, 0.05),
          mucat0 = seq(0.3, -0.2, by = -0.1),
          .RNG.name = "base::Mersenne-Twister", .RNG.seed = 12345)
  )
  fitjags <- jags.model ("model_ordinal_fixed.txt",
                         list ("n.r" = ncol(x)/3,
                               "n.c" = 3,
                               "n.cat" = 6,
                               "r.a" = r$r.a,
                               # "n.a" = r$n.a, # first fun.trans()
                               "n.ov" = r$n.ov, # second fun.trans()
                               "a.c" = a.c),
                         inits = inits,
                         n.chains = 2, quiet = T, n.adapt = 1000)
  
  update(fitjags, n.burnin, progress.bar = "none")
  fit.fixed <- coda.samples(fitjags, c("d", "dd", "best", "mucat"), n.iter, n.thin,
                            progress.bar = "none")

  sf <- summary(fit.fixed)[[1]]
  # cr1 <- format(round(sf[5, 1] + c(-1, 1) * qnorm(0.975) * sf[5, 2], 2), nsmall = 2)
  # cr2 <- format(round(sf[6, 1] + c(-1, 1) * qnorm(0.975) * sf[6, 2], 2), nsmall = 2)
  # cr3 <- format(round(sf[7, 1] + c(-1, 1) * qnorm(0.975) * sf[7, 2], 2), nsmall = 2)
  # cr1 <- paste(format(round(sf[5, 1], 2), nsmall = 2), " [", cr1[1], "; ", cr1[2], "]", sep = "")
  # cr2 <- paste(format(round(sf[6, 1], 2), nsmall = 2), " [", cr2[1], "; ", cr2[2], "]", sep = "")
  # cr3 <- paste(format(round(sf[7, 1], 2), nsmall = 2), " [", cr3[1], "; ", cr3[2], "]", sep = "")

  p1 <- pv(sf[5, 1], sf[5, 2])#, 12345) # for d[2]
  p2 <- pv(sf[6, 1], sf[6, 2])#, 12345) # for d[3]
  p3 <- pv(sf[7, 1], sf[7, 2])#, 12345) # for dd
  
  return(list("p_AB" = p1, "p_AC" = p2, "p_CB" = p3,
              "summary.fit" = summary(fit.fixed)))
}




# Random effects model for normal data
# ====================================

# Parameters:
# x:            data matrix
# n.burnin:     number of burn-in iterations
# n.iter:       number of MCMC iterations
# n.thin:       thinning rate

# Output:
# p_AB, p_AC, p_CB: p-values for comparisons of AB, AC, and BC, respectively
# summary.fit: summary of results based on MCMC iterations (i.e., posterior distributions)

random.ordinal <- function(x, n.burnin = 5000, n.iter = 25000, n.thin = 5) {

  options(warn = -1)
  suppressMessages(require(rjags))
  
  # Create one vector of all results:
  y <- c()
  for (i in 1:ncol(x)) {
    y <- c(y, x[, i])
  }
  
  # Add information of comparison and rater as separate columns:
  comp <- rep(1:(ncol(x)/3), each = nrow(x) * 3) # three columns are put below each other per rater
  rater <- rep(rep(1:3, each = nrow(x)), (ncol(x)/3)) # each question is asked in one row
  
  r <- fun.trans(x, comp, rater)
  
  # structure of r:
  # first dimension: questions (the object name "rater" above is misleading!)
  # second dimension: ratings - 1:
  # Sum of ratings:
  # If a rater says 3 (favors second option), then no count
  # If a rater says 2 (favors second option), then last element is count + 1.
  # If a rater says 1 (favors second option), then the last and the one before is count + 1.
  # ...
  # third dimension: raters (the object name "comp" above is misleading!)
  # Example:
  # r$r.a[1, 6, 2] describes the number of ratings of rater 2 in question 1 (unshaved vs.
  # mechanically shaved) that were 2 or lower (because we are looking at the 6th element of
  # the vector (-3, -2, -1, 0, 1, 2, 3).
  
  # 1=unshaved; 2=mechanically shaved; 3=digitally shaved:
  a.c <- matrix(c(1, 1, 3, 2, 3, 2), ncol = 2)
  
  require(coda)
  
  
  # random effects model:
  inits <- list(
    list (d = c(NA, 0.1, -0.1),
          mucat0 = seq(-0.3, 0.2, by = 0.1),
          # muc = c(-0.5, 0, 0.5),
          tau = c(0.9, 1, 1.1), # take this if you take "fit_ordinal_random.txt"
          # tau = c(NA, 0.1, 0.2), # delete this if you take "fit_ordinal_random.txt"
          # rho = c(NA, NA, 0.1), # delete this if you take "fit_ordinal_random.txt"
          .RNG.name = "base::Mersenne-Twister", .RNG.seed = 12345),
    list (d = c(NA, -0.05, 0.05),
          mucat0 = seq(0.3, -0.2, by = -0.1),
          # muc = c(0.5, 0, -0.5),
          tau = c(1.1, 0.9, 1), # take this if you take "fit_ordinal_random.txt"
          # tau = c(NA, 0.05, 0.1), # delete this if you take "fit_ordinal_random.txt"
          # rho = c(NA, NA, -0.1), # delete this if you take "fit_ordinal_random.txt"
          .RNG.name = "base::Mersenne-Twister", .RNG.seed = 12345)
  )
  fitjags <- jags.model("model_ordinal_random_clean.txt", # take model_ordinal_random.txt!
                        list ("n.r" = ncol(x)/3,
                              "n.c" = 3,
                              "n.cat" = 6,
                              "r.a" = r$r.a,
                              # "n.a" = r$n.a, # first fun.trans()
                              "n.ov" = r$n.ov, # second fun.trans()
                              "a.c" = a.c),
                        inits = inits,
                        n.chains = 2, quiet = T, n.adapt = 1000)
  update(fitjags, n.burnin, progress.bar = "none")
  fit.random <- coda.samples(fitjags, c("d", "dd", "best", "sdtau", "mucat"), n.iter, n.thin,
                             progress.bar = "none")
  
  sf <- summary(fit.random)[[1]]
  # cr1 <- format(round(sf[5, 1] + c(-1, 1) * qnorm(0.975) * sf[5, 2], 2), nsmall = 2)
  # cr2 <- format(round(sf[6, 1] + c(-1, 1) * qnorm(0.975) * sf[6, 2], 2), nsmall = 2)
  # cr3 <- format(round(sf[7, 1] + c(-1, 1) * qnorm(0.975) * sf[7, 2], 2), nsmall = 2)
  # cr1 <- paste(format(round(sf[5, 1], 2), nsmall = 2), " [", cr1[1], "; ", cr1[2], "]", sep = "")
  # cr2 <- paste(format(round(sf[6, 1], 2), nsmall = 2), " [", cr2[1], "; ", cr2[2], "]", sep = "")
  # cr3 <- paste(format(round(sf[7, 1], 2), nsmall = 2), " [", cr3[1], "; ", cr3[2], "]", sep = "")
  
  p1 <- pv(sf[5, 1], sf[5, 2])#, 12345) # for d[2]
  p2 <- pv(sf[6, 1], sf[6, 2])#, 12345) # for d[3]
  p3 <- pv(sf[7, 1], sf[7, 2])#, 12345) # for dd
  
  return(list("p_AB" = p1, "p_AC" = p2, "p_CB" = p3,
              "summary.fit" = summary(fit.random)))
  
}