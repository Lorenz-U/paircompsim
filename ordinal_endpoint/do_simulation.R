# Do simulation
###############

# Parameters:
# n.sim:    number of iterations in the simulation study
# model:    which model should be included ("fixed.ordinal", "random.ordinal",
#           "wilcox", "wilcox.clust")
# seed:     seed to make the simulation replicable
# alpha:    alpha level
# others:   see "generate_data.R" and "models.R".
# n.cores:  number of cores to be used (not in use, yet)

# Output:
# res:      a list of ouputs, one entry for each model:
#           [[1]] Bayesian fixed effects model
#           [[2]] Bayesian random effects model
#           [[3]] Wilcoxon test
#           [[4]] Wilcoxon test for clustered data
#           Each entry of the list is another list.
#           Each entry of that list represents the model results
#           based on a simulation iteration.
#           In res[[1]][[1]], the results of the Bayesian fixed effects model
#           based on the first simulation iteration is saved.
#           The whole output generated in the functions in "models.R"
#           is saved in that entry.
# p_res:    matrix of p-values; columns refer to the model, rows to the comparison/question

do.sim <- function(n.sim = 100, model = "fixed.ordinal", seed = 12345, alpha = 0.05,
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
  
  # create p-value matrix (ncol=4, because max. 4 different models will be evaluated).
  # Note: actually, we only assess 2 different models; for ordinal data, we need a
  # separate function => to be updated.
  p_AB_mat <- p_AC_mat <- p_CB_mat <- matrix(NA, nrow = n.sim, ncol = 4)
  
  # create object for p-value output:
  p_res <- data.frame("fixed.ordinal" = numeric(3), "random.ordinal" = numeric(3),
                      "wilcox" = numeric(3), "wilcox.clust" = numeric(3))
  row.names(p_res) <- c("p_AB", "p_AC", "p_CB")
  
  # Only needed for single core runs:
  # if ("fixed.ordinal" %in% model) {
  #   res.fix.ordinal <- list(n.sim)
  # }
  # if ("random.ordinal" %in% model) {
  #   res.ran.ordinal <- list(n.sim)
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
  if ("fixed.ordinal" %in% model) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Started fixed models; working on ", getDoParWorkers(), " cores\n"))
    res.fix.ordinal <- foreach(sim.i = 1:n.sim) %dopar% {
      source("models.R")
      require(rjags)
      fixed.ordinal(y[[sim.i]], n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin)
    }
    stopCluster(cl)
    # save p-values: first column of the p-value matrix is reserved for fixed.ordinal model:
    for (sim.i in 1:n.sim) {
      p_AB_mat[sim.i, 1] <- res.fix.ordinal[[sim.i]]$p_AB <= alpha
      p_AC_mat[sim.i, 1] <- res.fix.ordinal[[sim.i]]$p_AC <= alpha
      p_CB_mat[sim.i, 1] <- res.fix.ordinal[[sim.i]]$p_CB <= alpha
    }
    cat(paste0("Fixed models done, ", Sys.time(), "\n"))
  }
  if ("random.ordinal" %in% model) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Started random models; working on ", getDoParWorkers(), " cores\n"))
    res.ran.ordinal <- foreach(sim.i = 1:n.sim) %dopar% {
      source("models.R")
      require(rjags)
      random.ordinal(y[[sim.i]], n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin)
    }
    stopCluster(cl)
    # save p-values: second column of the p-value matrix is reserved for random.ordinal model:
    for (sim.i in 1:n.sim) {
      p_AB_mat[sim.i, 2] <- res.ran.ordinal[[sim.i]]$p_AB <= alpha
      p_AC_mat[sim.i, 2] <- res.ran.ordinal[[sim.i]]$p_AC <= alpha
      p_CB_mat[sim.i, 2] <- res.ran.ordinal[[sim.i]]$p_CB <= alpha
    }
    cat(paste0("Random models done, ", Sys.time(), "\n"))
  }
  
  
  # Do Wilcoxon test for A vs. C:
  
  if ("wilcox" %in% model) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Started Wilcoxon tests; working on ", getDoParWorkers(), " cores\n"))
    res.wilcox <- foreach(sim.i = 1:n.sim) %dopar% {
      y.i <- y[[sim.i]]
      y.i.i <- c()
      for (vec.i in 1:(ncol(y.i)/3)) {
        y.i.i <- c(y.i.i, y.i[, (vec.i-1)*3 + 2])
      }
      wilcox.test(y.i.i)$p.value
    }
    stopCluster(cl)
    for (sim.i in 1:n.sim) {
      p_AB_mat[sim.i, 3] <- NA
      p_AC_mat[sim.i, 3] <- res.wilcox[[sim.i]] <= alpha
      p_CB_mat[sim.i, 3] <- NA
    }
    cat(paste0("Wilcoxon tests done, ", Sys.time(), "\n"))
  }
  
  # Fit Wilcoxon cluster test for A vs. C:
  
  if ("wilcox.clust" %in% model) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Started Wilcoxon cluster tests; working on ", getDoParWorkers(), " cores\n"))
    res.wilcox.clust <- foreach(sim.i = 1:n.sim) %dopar% {
      y.i <- y[[sim.i]]
      y.i.i <- c()
      for (vec.i in 1:(ncol(y.i)/3)) {
        y.i.i <- c(y.i.i, y.i[, (vec.i-1)*3 + 2])
      }
      rater.i <- rep(1:(ncol(y.i)/3), each = nrow(y.i))
      require(clusrank)
      clusWilcox.test(y.i.i, cluster = rater.i, paired = T)$p.value
    }
    stopCluster(cl)
    for (sim.i in 1:n.sim) {
      p_AB_mat[sim.i, 4] <- NA
      p_AC_mat[sim.i, 4] <- res.wilcox.clust[[sim.i]] <= alpha
      p_CB_mat[sim.i, 4] <- NA
    }
    cat(paste0("Wilcoxon cluster tests done, ", Sys.time(), "\n"))
  }
  # stopCluster(cl)
  
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
  #   if ("fixed.ordinal" %in% model) {
  #     res.fix.ordinal[[sim.i]] <- fixed.ordinal(y, n.burnin = n.burnin,
  #                                             n.iter = n.iter, n.thin = n.thin)
  #     # save p-values: first column of the p-value matrix is reserved for fixed.ordinal model:
  #     p_AB_mat[sim.i, 1] <- res.fix.ordinal[[sim.i]]$p_AB <= alpha
  #     p_AC_mat[sim.i, 1] <- res.fix.ordinal[[sim.i]]$p_AC <= alpha
  #     p_CB_mat[sim.i, 1] <- res.fix.ordinal[[sim.i]]$p_CB <= alpha
  #   }
  #   if ("random.ordinal" %in% model) {
  #     res.ran.ordinal[[sim.i]] <- random.ordinal(y, n.burnin = n.burnin,
  #                                             n.iter = n.iter, n.thin = n.thin)
  #     # save p-values: second column of the p-value matrix is reserved for random.ordinal model:
  #     p_AB_mat[sim.i, 2] <- res.ran.ordinal[[sim.i]]$p_AB <= alpha
  #     p_AC_mat[sim.i, 2] <- res.ran.ordinal[[sim.i]]$p_AC <= alpha
  #     p_CB_mat[sim.i, 2] <- res.ran.ordinal[[sim.i]]$p_CB <= alpha
  #   }
  # }
  if ("fixed.ordinal" %in% model) {
    # save results in p-value data set and in res list:
    p_res$fixed.ordinal <- c(mean(p_AB_mat[, 1]), mean(p_AC_mat[, 1]), mean(p_CB_mat[, 1]))
    res[[1]] <- res.fix.ordinal
  }
  if ("random.ordinal" %in% model) {
    # save results in p-value data set and in res list:
    p_res$random.ordinal <- c(mean(p_AB_mat[, 2]), mean(p_AC_mat[, 2]), mean(p_CB_mat[, 2]))
    res[[2]] <- res.ran.ordinal
  }
  if ("wilcox" %in% model) {
    # save results in p-value data set and in res list:
    p_res$wilcox <- c(mean(p_AB_mat[, 3]), mean(p_AC_mat[, 3]), mean(p_CB_mat[, 3]))
    res[[3]] <- res.wilcox
  }
  if ("wilcox.clust" %in% model) {
    # save results in p-value data set and in res list:
    p_res$wilcox.clust <- c(mean(p_AB_mat[, 4]), mean(p_AC_mat[, 4]), mean(p_CB_mat[, 4]))
    res[[4]] <- res.wilcox.clust
  }
  return(list("res" = res, "p_res" = p_res))
}
