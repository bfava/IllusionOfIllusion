
# Open Packages -----------------------------------------------------------

library(tidyverse)


# Open Datasets -----------------------------------------------------------

  # Datasets are now available from GLP at
  # https://www.econometricsociety.org/publications/econometrica/2021/09/01/economic-predictions-big-data-illusion-sparsity

# Macro 1
macro1_df <- R.matlab::readMat("forJunnanHe/FredMDlargeHor1.mat")
macro1_X <- macro1_df$X
macro1_y <- macro1_df$Y
macro1_U <- rep(1, nrow(macro1_X))

# Macro 2
macro2_df <- read_csv2("forJunnanHe/GrowthChern.csv")
macro2_df %<>%
  select(-intercept) %>% 
  mutate_all(function(x) (x - mean(x)) / sd(x))
macro2_X <- macro2_df[, -1] %>% as.matrix()
macro2_y <- macro2_df[[1]]
macro2_U <- rep(1, nrow(macro2_X))

# Finance 1
fin1_df <- R.matlab::readMat("forJunnanHe/Goyal.mat")
fin1_y <- fin1_df$y
fin1_X <- fin1_df$x
fin1_U <- rep(1, nrow(fin1_X))

# Finance 2

# Micro 1
micro1_df <- R.matlab::readMat("forJunnanHe/Abortion_data.mat")
micro1_y <- micro1_df$DyM %>% {. / sd(.)}
micro1_X <- micro1_df$Zmurd %>% as.tibble() %>% mutate_all(function(x) x / sd(x)) %>% as.matrix()
micro1_U <- micro1_df$DxM %>% {. / sd(.)}

# Micro 2
micro2_df <- R.matlab::readMat("forJunnanHe/Data1stStageEminentDomain.mat")
micro2_y <- micro2_df[[1]][, 1]
micro2_X <- micro2_df[[1]][, -1]
micro2_U <- rep(1, nrow(micro2_X))
