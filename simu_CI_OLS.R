#' Script for simulations NAVAE CI case of OLS
#' 
#' 

#==============================================================================#
# Paths ----------
#==============================================================================#

# Temp (later, a package), path to the R code for the CI

# Select the correct path to main folder depending on user/setting
# Setting: CREST, Lucas's computer
path_dropbox <- file.path("C:/Temp/Dropbox")
# Setting: personal computer (Lucas)
# path_dropbox <- file.path("/Users/lucasgirard/Dropbox")

path_lib_code <- file.path(path_dropbox, "finite_CI_perso", "NavaeCI", "R")

path_res_simu <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "res_simu")
if (!dir.exists(path_res_simu)) {
  dir.create(path_res_simu)
}  

#==============================================================================#
# Libraries and lib files ----------
#==============================================================================#

source(file.path(path_lib_code, "CI_mean.R"))
source(file.path(path_lib_code, "CI_OLS.R"))

library(purrr)
library(dplyr)

#==============================================================================#
# Temp to write and debug the function ----------
#==============================================================================#


n <- 1000

correlation_X1_X2 <- 0.5
variance_X1 <- 1
variance_X2 <- 2
covariance_X1_X2 <- correlation_X1_X2 * sqrt(variance_X1) * sqrt(variance_X2)
Sigma_matrix <- matrix(c(variance_X1, covariance_X1_X2, covariance_X1_X2, variance_X2), ncol = 2)

esperance_X1 <- 0
esperance_X2 <- 0



beta0 <- c(2, 1, -3)

eps <- rnorm(n = n)

X <- drawn_X 

Y <- c(beta0[[1]] + drawn_X %*% beta0[2:3] + eps)


# default argument
alpha = 0.05
a = NULL
power_of_n_for_b = NULL
omega = NULL
power_of_n_for_omega = NULL
bounds = list(lambda_reg = NULL,
              K_reg = NULL,
              K_eps = NULL,
              K_xi = NULL,
              C = NULL,
              B = NULL)
K_xi = NULL
param_BE_EE = list(
  choice = "best",
  setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
  regularity = list(C0 = 1, p = 2),
  eps = 0.1)
options = list(center = FALSE, bounded_case = FALSE, with_Exp_regime = FALSE)
matrix_u = diag(NCOL(X) + 1)



i <- 3
beta0[[2]] * drawn_X[i,1] + beta0[[3]] * drawn_X[i,2]
(drawn_X %*% beta0[2:3])[i,]

2 * -0.07 - 3 * 1.4

#==============================================================================#
# Function to draw data ----------
#==============================================================================#

Draw_data_regression_X_Normal_eps_Gumbel <- function(
  n, 
  esperance_X1 = 0, esperance_X2 = 0,
  variance_X1 = 1, variance_X2 = 2,
  correlation_X1_X2 = 0.5,
  beta0 = c(2, 1, -3))
{
  
  # Compute the variance-covariance matrix to draw X1 and X2
  covariance_X1_X2 <- correlation_X1_X2 * sqrt(variance_X1) * sqrt(variance_X2)
  Sigma_matrix <- matrix(c(variance_X1, covariance_X1_X2, covariance_X1_X2, variance_X2), ncol = 2)
  
  X <- MASS::mvrnorm(n = n, mu = c(esperance_X1, esperance_X2), Sigma = Sigma_matrix)
  
  # Draw from Gumbel distribution
  # Trick: https://en.wikipedia.org/wiki/Gumbel_distribution see generation
  # from a uniform(0, 1)
  
  # Euler-Mascheroni's constant used in our distribution 
  EM_constant <- -digamma(1)
  
  vec_beta_param_Gumbel <- sqrt((6 / pi^2) * abs(rowSums(X)))
  
  vec_mu_param_Gumbel <- -EM_constant * vec_beta_param_Gumbel

  auxiliary_unif <- runif(n = n)
  
  eps <- vec_mu_param_Gumbel - vec_beta_param_Gumbel * log(-log(auxiliary_unif))
  
  Y <- c(beta0[[1]] + X %*% beta0[2:3] + eps)
  
  return(list(X = X, Y = Y))
}


#==============================================================================#
# Trials of the function ----------
#==============================================================================#

n <- 500 * 10^3

mydata <- Draw_data_regression_X_Normal_eps_Gumbel(n = n)

res <- Navae_ci_ols(Y = mydata$Y, X = mydata$X, 
                    alpha = 0.59, K_xi = 9, verbose = TRUE,
                    options = list(center = FALSE, bounded_case = FALSE, with_Exp_regime = TRUE))

res


