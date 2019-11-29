
# Estimation Function - t student -----------------------------------------

SaS_t <- function(X, y, U, nu, S = 10000L, a = 1L, b = 1L, A = 1L, B = 1L){
  
  # Open Packages
  library(tidyverse)
  library(extraDistr)
  library(magrittr)
  library(fBasics)
  library(HDCI)
  
  # Correct Formats
  U <- as.matrix(U)
  
  # Set Start
  betan <- paste0("X", 1:ncol(X))
  colnames(X) <- betan
  
  #lambda <- escv.glmnet(X, y)$lambda.cv1se
  #lasso.mod <- Lasso(X, y, lambda = lambda * 1000)
  lasso.mod <- Lasso(X, y, fix.lambda = F)
  
  beta <- lasso.mod$beta
  names(beta) <- betan
  phi <- lasso.mod$beta0
  sigma2 <- var(mypredict(lasso.mod, X) - y)
  z <- as.integer(beta != 0)
  k <- length(z)
  Tn <- dim(X)[1]
  
  lambda2 <- rinvgamma(k, nu / 2, nu / 2)
  lambda2_bar <- nu / (nu - 2)
  beta_s <- double(k)
  names(beta_s) <- betan
  
  # Create Objects 
  x <- c(seq(0, 0.1, 0.001), seq(0.11, 0.9, 0.01), seq(0.901, 1, 0.001))
  x_c <- x[-length(x)] + (diff(x) / 2)
  xy_c <- expand.grid(x_c, x_c)
  
  grid_q    <- xy_c$Var1
  grid_r2   <- xy_c$Var2
  grid_area <- apply(expand.grid(diff(x), diff(x)), 1, prod)
  grid_lq   <- log(grid_q)
  grid_l1q  <- log(1 - grid_q)
  grid_lr2  <- log(grid_r2)
  grid_l1r2 <- log(1 - grid_r2)
  grid_c1   <- k * grid_q * ((1 - grid_r2) / (2 * grid_r2))
  grid_pr2q <- rep(0, length(grid_q))
  
  all_pos <- seq_along(grid_q)
  
  dfr <- tibble(q       = double(S),
                r2      = double(S),
                phi     = rep(list(double(k)), S),
                z       = rep(list(double(k)), S),
                beta    = rep(list(double(k)), S),
                sigma2  = double(S),
                lambda2 = rep(list(double(k)), S))
  
  for (s in seq_along(1:S)){
    
    # Sample p(r2, q | Y, phi, beta, sigma2, z, lambda2)
    cterm <- c(t(beta) %*% diag(z) %*% beta)
    grid_pr2q <- pqr(sigma2, k, grid_lq, grid_l1q, grid_lr2, grid_l1r2, grid_c1 * lambda2_bar, grid_area, cterm, z)
    
    pos <- sample(all_pos, 1, prob = grid_pr2q)
    q   <- grid_q[pos]
    r2  <- grid_r2[pos]
    
    # Sample p(phi | Y, z, beta, g, q, sigma, lambda2)
    phi <- rnorm(1, solve(t(U) %*% U) %*% t(U) %*% (y - X %*% beta), sigma2 * solve(t(U) %*% U))
    
    # Sample p(z | Y, phi, r2, q, lambda2)
    
    for (j in sample(seq_along(z))){
      
      if (!(sum(z) == 1 & z[j] == 1)){
        zp <- z
        zp[j] <- 1 - zp[j]
        probzp <- pz_t(z, zp, r2, q, lambda2_bar, y, X, U, phi, Tn, k)
        z[j] <- sample(c(z[j], zp[j]), 1, prob = c(1, probzp))
      }
      
    }
    
    # sample p(sigma2 | y, phi, r2, q, z, lambda2)
    g2 <- exp(log(r2) - (log(k) + log(q) + log(lambda2_bar) + log(1 - r2)))
    tauz <- sum(z)
    ytil <- y - U * phi
    Xtil <- X[, z == 1]
    Wtil <- t(Xtil) %*% Xtil + diag(tauz) / g2
    betatilhat <- solve(Wtil) %*% t(Xtil) %*% ytil
    
    sigma2 <- rinvgamma(1, Tn / 2, (t(ytil) %*% ytil - t(betatilhat) %*% (t(Xtil) %*% Xtil + diag(tauz) / g2) %*% betatilhat) / 2)
    
    # sample p(beta | y, phi, sigma2, r2, q, z, lambda2)
    beta_s[1:k] <- 0
    beta_s[z == 1] <- MASS::mvrnorm(1, betatilhat, sigma2 * solve(t(Xtil) %*% Xtil + diag(tauz) / g2))
    
    beta <- beta_s * sqrt(lambda2)
    
    # sample p(lambda | nu, beta, sigma2, r2)
    lambda2 <- map2_dbl(rep((nu + 1) / 2, k), (nu + (beta ^ 2 / (sigma2 * g2))) / 2, rinvgamma, n = 1)
    #lambda2_bar <- mean(lambda2)
    
    # sample p(phi, beta, sigma2, r2, z, q, lambda | y, X)
    dfr[s, 1] <- q
    dfr[s, 2] <- r2
    dfr$phi[[s]] <- phi
    dfr$z[[s]] <- z
    dfr$beta[[s]] <- beta
    dfr[s, 6] <- sigma2
    dfr$lambda2[[s]] <- lambda2
    
  }
  
  return(dfr)
  
}


# Estimation Function - t student without phi -----------------------------

SaS_t_nphi <- function(X, y, nu, S = 10000L, a = 1L, b = 1L, A = 1L, B = 1L){
  
  # Open Packages
  library(tidyverse)
  library(extraDistr)
  library(magrittr)
  library(fBasics)
  library(HDCI)
  
  # Set Start
  betan <- paste0("X", 1:ncol(X))
  colnames(X) <- betan
  
  #lambda <- escv.glmnet(X, y)$lambda.cv1se
  #lasso.mod <- Lasso(X, y, lambda = lambda * 1000)
  lasso.mod <- Lasso(X, y, fix.lambda = F)
  
  beta <- lasso.mod$beta
  names(beta) <- betan
  sigma2 <- var(mypredict(lasso.mod, X) - y)
  z <- as.integer(beta != 0)
  k <- length(z)
  Tn <- dim(X)[1]
  
  lambda2 <- rinvgamma(k, nu / 2, nu / 2)
  lambda2_bar <- nu / (nu - 2)
  beta_s <- double(k)
  names(beta_s) <- betan
  
  # Create Objects 
  x <- c(seq(0, 0.1, 0.001), seq(0.11, 0.9, 0.01), seq(0.901, 1, 0.001))
  x_c <- x[-length(x)] + (diff(x) / 2)
  xy_c <- expand.grid(x_c, x_c)
  
  grid_q    <- xy_c$Var1
  grid_r2   <- xy_c$Var2
  grid_area <- apply(expand.grid(diff(x), diff(x)), 1, prod)
  grid_lq   <- log(grid_q)
  grid_l1q  <- log(1 - grid_q)
  grid_lr2  <- log(grid_r2)
  grid_l1r2 <- log(1 - grid_r2)
  grid_c1   <- k * grid_q * ((1 - grid_r2) / (2 * grid_r2))
  grid_pr2q <- rep(0, length(grid_q))
  
  all_pos <- seq_along(grid_q)
  
  dfr <- tibble(q       = double(S),
                r2      = double(S),
                z       = rep(list(double(k)), S),
                beta    = rep(list(double(k)), S),
                sigma2  = double(S),
                lambda2 = rep(list(double(k)), S))
  
  for (s in seq_along(1:S)){
    
    # Sample p(r2, q | Y, beta, sigma2, z, lambda2)
    cterm <- c(t(beta) %*% diag(z) %*% beta)
    grid_pr2q <- pqr(sigma2, k, grid_lq, grid_l1q, grid_lr2, grid_l1r2, grid_c1 * lambda2_bar, grid_area, cterm, z)
    
    pos <- sample(all_pos, 1, prob = grid_pr2q)
    q   <- grid_q[pos]
    r2  <- grid_r2[pos]
    
    # Sample p(z | Y, r2, q, lambda2)
    
    for (j in sample(seq_along(z))){
      
      if (!(sum(z) == 1 & z[j] == 1)){
        zp <- z
        zp[j] <- 1 - zp[j]
        probzp <- pz_t_nphi(z, zp, r2, q, lambda2_bar, y, X, Tn, k)
        z[j] <- sample(c(z[j], zp[j]), 1, prob = c(1, probzp))
      }
      
    }
    
    # sample p(sigma2 | y, r2, q, z)
    #g2 <- (1 / (k * q * lambda2_bar)) * (r2 / (1 - r2))
    g2 <- exp(log(r2) - (log(k) + log(q) + log(lambda2_bar) + log(1 - r2)))
    tauz <- sum(z)
    ytil <- y
    Xtil <- X[, z == 1]
    Wtil <- t(Xtil) %*% Xtil + diag(tauz) / g2
    betatilhat <- solve(Wtil) %*% t(Xtil) %*% ytil
    
    sigma2 <- rinvgamma(1, Tn / 2, (t(ytil) %*% ytil - t(betatilhat) %*% (t(Xtil) %*% Xtil + diag(tauz) / g2) %*% betatilhat) / 2)
    
    # sample p(beta | y, sigma2, r2, q, z, lambda2)
    beta_s[1:k] <- 0
    beta_s[z == 1] <- MASS::mvrnorm(1, betatilhat, sigma2 * solve(t(Xtil) %*% Xtil + diag(tauz) / g2))
    
    beta <- beta_s * sqrt(lambda2)
    
    # sample p(lambda2 | nu, beta, sigma2, r2)
    lambda2 <- map2_dbl(rep((nu + 1) / 2, k), (nu + (beta ^ 2 / (sigma2 * g2))) / 2, rinvgamma, n = 1)
    #lambda2_bar <- mean(lambda2)
    
    # sample p(beta, sigma2, r2, z, q, lambda2 | y, X)
    dfr[s, 1] <- q
    dfr[s, 2] <- r2
    dfr$z[[s]] <- z
    dfr$beta[[s]] <- beta
    dfr[s, 5] <- sigma2
    dfr$lambda2[[s]] <- lambda2
    
  }
  
  return(dfr)
  
}
