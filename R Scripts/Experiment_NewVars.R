
# Open Packages -----------------------------------------------------------

library(tidyverse)


# Create new datasets -----------------------------------------------------

set.seed(1234)
old_datasets <- list(macro1_X, macro2_X, fin1_X, micro1_X, micro2_X)
new_datasets <- list()

n_newvars <- 2

for (i in 1:length(old_datasets)){

  new_datasets[[i]] <- old_datasets[[i]] %>% 
    {cbind(., 
           scale(
             matrix(
               rnorm(n_newvars * nrow(.)), # Draw from normal distribution
               ncol = n_newvars)))} 
  
}


# Run Spike and Slab ------------------------------------------------------

datasets_y <- list(macro1_y, macro2_y, fin1_y, micro1_y, micro2_y)
datasets_U <- list(NA, macro2_U, NA, micro1_U, NA)

randomvar_dfr <- list()

for (i in 1:length(new_datasets)){

  if(is.na(datasets_U[[i]][1])){
    
    randomvar_dfr[[i]] <- SaS_nphi(new_datasets[[i]], datasets_y[[i]], S = 10000)  
    
  } else {
    
    randomvar_dfr[[i]] <- SaS(new_datasets[[i]], datasets_y[[i]], datasets_U[[i]], S = 10000)
    
  }
  
  cat(i, "/", length(new_datasets), " concluded\n")
  
}

