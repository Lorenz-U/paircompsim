# Do simulation
###############

# Parameters:
# n.sim:    number of iterations in the simulation study
# model:    which model should be included ("fixed.normal", "random.normal",
#           "t.test", "fr.mixed")
# seed:     seed to make the simulation replicable
# alpha:    alpha level
# others:   see "generate_data.R" and "models.R".
# n.cores:  number of cores to be used (not in use, yet)

# Output:
# res:      a list of ouputs, one entry for each model:
#           [[1]] Bayesian fixed effects model
#           [[2]] Bayesian random effects model
#           [[3]] t-test
#           [[4]] frequentist mixed effects model
#           Each entry of the list is another list.
#           Each entry of that list represents the model results
#           based on a simulation iteration.
#           In res[[1]][[1]], the results of the Bayesian fixed effects model
#           based on the first simulation iteration is saved.
#           The whole output generated in the functions in "models.R"
#           is saved in that entry.
# p_res:    matrix of p-values; columns refer to the model, rows to the comparison/question

do.sim <- function(n.sim = 100, model = "fixed.normal", seed = 12345, alpha = 0.05,
                   n = 59, n_raters = 6, n_questions = 3,
                   d_AB = 0, d_AC = 0, d_CB = 0,
                   sigma_AB = 1, sigma_AC = 1, sigma_CB = 1,
                   tau_AB = 0, tau_AC = 0, tau_CB = 0,
                   n.burnin = 5000, n.iter = 25000, n.thin = 5,
                   n.cores = 4) {
  cat(paste0("Started at ", Sys.time(), "\n"))
  set.seed(seed)
  # create list to safe all fit summaries per model and iteration;
  # one list entry per model:
  res <- list(4)
  
  # create p-value matrix.
  p_AB_mat <- p_AC_mat <- p_CB_mat <- matrix(NA, nrow = n.sim, ncol = 4)
  
  # create object for p-value output:
  p_res <- data.frame("fixed.normal" = numeric(3), "random.normal" = numeric(3),
                      "t.test" = numeric(3), "fr.mixed" = numeric(3))
  row.names(p_res) <- c("p_AB", "p_AC", "p_CB")
  
  # Only needed for single core runs:
  # if ("fixed.normal" %in% model) {
  #   res.fix.normal <- list(n.sim)
  # }
  # if ("random.normal" %in% model) {
  #   res.ran.normal <- list(n.sim)
  # }
  
  # For parallel runs:
  require(doParallel)
  cl <- makeCluster(n.cores)#, type = "PSOCK")
  registerDoParallel(cl)
  cat(paste0("Working on ", getDoParWorkers(), " cores\n"))

  # generate data
  cat(paste0("Started data simulation; working on ", getDoParWorkers(), " cores\n"))
  y <- foreach(icount(n.sim)) %dopar% {
    source("generate_data.R")
    sim.dat(n = n, n_raters = n_raters, n_questions = n_questions,
            d_AB = d_AB, d_AC = d_AC, d_CB = d_CB,
            sigma_AB = sigma_AB, sigma_AC = sigma_AC, sigma_CB = sigma_CB,
            tau_AB = tau_AB, tau_AC = tau_AC, tau_CB = tau_CB)
  }
  stopCluster(cl)
  cat(paste0("Data simulation done, ", Sys.time(), "\n"))
  # check which models shall be evaluated and run those models:
  if ("fixed.normal" %in% model) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Started fixed models; working on ", getDoParWorkers(), " cores\n"))
    res.fix.normal <- foreach(sim.i = 1:n.sim) %dopar% {
      source("models.R")
      require(rjags)
      fixed.normal(y[[sim.i]], n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin)
    }
    stopCluster(cl)
    # save p-values: first column of the p-value matrix is reserved for fixed.normal model:
    for (sim.i in 1:n.sim) {
      p_AB_mat[sim.i, 1] <- res.fix.normal[[sim.i]]$p_AB <= alpha
      p_AC_mat[sim.i, 1] <- res.fix.normal[[sim.i]]$p_AC <= alpha
      p_CB_mat[sim.i, 1] <- res.fix.normal[[sim.i]]$p_CB <= alpha
    }
    cat(paste0("Fixed models done, ", Sys.time(), "\n"))
  }
  if ("random.normal" %in% model) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Started random models; working on ", getDoParWorkers(), " cores\n"))
    res.ran.normal <- foreach(sim.i = 1:n.sim) %dopar% {
      source("models.R")
      require(rjags)
      random.normal(y[[sim.i]], n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin)
    }
    stopCluster(cl)
    # save p-values: second column of the p-value matrix is reserved for random.normal model:
    for (sim.i in 1:n.sim) {
      p_AB_mat[sim.i, 2] <- res.ran.normal[[sim.i]]$p_AB <= alpha
      p_AC_mat[sim.i, 2] <- res.ran.normal[[sim.i]]$p_AC <= alpha
      p_CB_mat[sim.i, 2] <- res.ran.normal[[sim.i]]$p_CB <= alpha
    }
    cat(paste0("Random models done, ", Sys.time(), "\n"))
  }
  
  
  # Do frequentist t-test for A vs. C:
  if ("t.test" %in% model) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Started t-tests; working on ", getDoParWorkers(), " cores\n"))
    res.t.test <- foreach(sim.i = 1:n.sim) %dopar% {
      y.i <- y[[sim.i]]
      y.i.i <- c()
      for (vec.i in 1:dim(y.i)[3]) {
        # y.i.i <- c(y.i[, 2, 1], y.i[, 2, 2], y.i[, 2, 3], y.i[, 2, 4], y.i[, 2, 5], y.i[, 2, 6])
        y.i.i <- c(y.i.i, y.i[, 2, vec.i])
      }
      t.test(y.i.i)$p.value
    }
    stopCluster(cl)
    for (sim.i in 1:n.sim) {
      p_AB_mat[sim.i, 3] <- NA
      p_AC_mat[sim.i, 3] <- res.t.test[[sim.i]] <= alpha
      p_CB_mat[sim.i, 3] <- NA
    }
    cat(paste0("t-tests done, ", Sys.time(), "\n"))
  }
  
  # Fit frequentist mixed model for A vs. C:
  if ("fr.mixed" %in% model) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Started fr. mixed models; working on ", getDoParWorkers(), " cores\n"))
    res.fr.mixed <- foreach(sim.i = 1:n.sim) %dopar% {
      y.i <- y[[sim.i]]
      y.i.i <- c()
      for (vec.i in 1:dim(y.i)[3]) {
        y.i.i <- c(y.i.i, y.i[, 2, vec.i])
      }
      rater.i <- rep(1:dim(y.i)[3], each = dim(y.i)[1])
      require(nlme)
      fit <- lme(y.i.i ~ 1, random = ~1|rater.i) # compound symmetry structure
      coef(summary(fit))[1, 5]
    }
    stopCluster(cl)
    for (sim.i in 1:n.sim) {
      p_AB_mat[sim.i, 4] <- NA
      p_AC_mat[sim.i, 4] <- res.fr.mixed[[sim.i]] <= alpha
      p_CB_mat[sim.i, 4] <- NA
    }
    cat(paste0("Frequentist mixed models done, ", Sys.time(), "\n"))
  }
  # End of parallel run section.
  
  # For single core runs:
  # for (sim.i in 1:n.sim) {
  #   # generate data:
  #   y <- sim.dat(n = n, n_raters = n_raters, n_questions = n_questions,
  #                d_AB = d_AB, d_AC = d_AC, d_CB = d_CB,
  #                sigma_AB = sigma_AB, sigma_AC = sigma_AC, sigma_CB = sigma_CB,
  #                tau_AB = tau_AB, tau_AC = tau_AC, tau_CB = tau_CB)
  #   
  #   # check which models shall be evaluated and run those models:
  #   if ("fixed.normal" %in% model) {
  #     res.fix.normal[[sim.i]] <- fixed.normal(y, n.burnin = n.burnin,
  #                                             n.iter = n.iter, n.thin = n.thin)
  #     # save p-values: first column of the p-value matrix is reserved for fixed.normal model:
  #     p_AB_mat[sim.i, 1] <- res.fix.normal[[sim.i]]$p_AB <= alpha
  #     p_AC_mat[sim.i, 1] <- res.fix.normal[[sim.i]]$p_AC <= alpha
  #     p_CB_mat[sim.i, 1] <- res.fix.normal[[sim.i]]$p_CB <= alpha
  #   }
  #   if ("random.normal" %in% model) {
  #     res.ran.normal[[sim.i]] <- random.normal(y, n.burnin = n.burnin,
  #                                             n.iter = n.iter, n.thin = n.thin)
  #     # save p-values: second column of the p-value matrix is reserved for random.normal model:
  #     p_AB_mat[sim.i, 2] <- res.ran.normal[[sim.i]]$p_AB <= alpha
  #     p_AC_mat[sim.i, 2] <- res.ran.normal[[sim.i]]$p_AC <= alpha
  #     p_CB_mat[sim.i, 2] <- res.ran.normal[[sim.i]]$p_CB <= alpha
  #   }
  # }
  if ("fixed.normal" %in% model) {
    # save results in p-value data set and in res list:
    p_res$fixed.normal <- c(mean(p_AB_mat[, 1]), mean(p_AC_mat[, 1]), mean(p_CB_mat[, 1]))
    res[[1]] <- res.fix.normal
  }
  if ("random.normal" %in% model) {
    # save results in p-value data set and in res list:
    p_res$random.normal <- c(mean(p_AB_mat[, 2]), mean(p_AC_mat[, 2]), mean(p_CB_mat[, 2]))
    res[[2]] <- res.ran.normal
  }
  if ("t.test" %in% model) {
    # save results in p-value data set and in res list:
    p_res$t.test <- c(mean(p_AB_mat[, 3]), mean(p_AC_mat[, 3]), mean(p_CB_mat[, 3]))
    res[[3]] <- res.t.test
  }
  if ("fr.mixed" %in% model) {
    # save results in p-value data set and in res list:
    p_res$fr.mixed <- c(mean(p_AB_mat[, 4]), mean(p_AC_mat[, 4]), mean(p_CB_mat[, 4]))
    res[[4]] <- res.fr.mixed
  }
  return(list("res" = res, "p_res" = p_res))
}