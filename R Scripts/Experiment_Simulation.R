
# Simulate Functions ------------------------------------------------------

SimulateDataset <- function(n, k, relevant, sigma_conf, seed1, seed2, seed3) {
  
  set.seed(seed1)
  sim_X <- apply(matrix(rnorm(n * k), nrow = n, k), 2, function(x) (x - mean(x)) / sd(x))

  set.seed(seed2)
  beta_t <- runif(relevant, -1, 1)
  
  set.seed(seed3)
  conf <- sigma_conf * rnorm(n, 0)
  
  sim_y <- cbind(sim_X[,1:relevant], conf) %*% c(beta_t, 1)
  sim_y <- (sim_y - mean(sim_y)) / sd(sim_y)
  
  sim <- list()
  sim[[1]] <- sim_X
  sim[[2]] <- sim_y

  return(sim)
  
}

SimulateResults <- function(sim_X, sim_y, nu_list, size = 10000L){
  
  sim_dfr <- list()
  i <- 1
  
  # t-student
  for (nu in nu_list){
    
    sim_dfr[[i]] <- SaS_t_nphi(sim_X, sim_y, nu, size)
    cat("nu = ", nu, " complete\n")
    i <- i + 1
    
  }
  
  sim_dfr[[i]] <- SaS_nphi(sim_X, sim_y, size)
  cat("Normal complete\n")
  return(sim_dfr)
  
}


# Simulate - 3 relevants --------------------------------------------------

sim_list <- list()

for (i in 1:6){
  
  sim_list[[i]] <- SimulateDataset(68, 16, 3, 0.0 + 0.75 * i, 12345, 12, 182)
  
}

nu_list <- c(4, 15)
sim_list_dfr <- list()

sim_list_dfr[[1]] <- SimulateResults(sim_list[[1]][[1]], sim_list[[1]][[2]], nu_list, 10000)
sim_list_dfr[[2]] <- SimulateResults(sim_list[[2]][[1]], sim_list[[2]][[2]], nu_list, 10000)
sim_list_dfr[[3]] <- SimulateResults(sim_list[[3]][[1]], sim_list[[3]][[2]], nu_list, 10000)
sim_list_dfr[[4]] <- SimulateResults(sim_list[[4]][[1]], sim_list[[4]][[2]], nu_list, 10000)
sim_list_dfr[[5]] <- SimulateResults(sim_list[[5]][[1]], sim_list[[5]][[2]], nu_list, 10000)
sim_list_dfr[[6]] <- SimulateResults(sim_list[[6]][[1]], sim_list[[6]][[2]], nu_list, 10000)

