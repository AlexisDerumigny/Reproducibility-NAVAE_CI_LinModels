#' Script for simulations NAVAE CI case of OLS
#' 
#' 

#==============================================================================#
# Paths ----------
#==============================================================================#

# Temp (later, a package), path to the R code for the CI

# Select the correct path to main folder depending on user/setting
# Setting: CREST, Lucas's computer
# path_dropbox <- file.path("C:/Temp/Dropbox")
# Setting: personal computer (Lucas)
path_dropbox <- file.path("/Users/lucasgirard/Dropbox")

path_lib_code <- file.path(path_dropbox, "finite_CI_perso", "NavaeCI", "R")

path_res_simu_OLS <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "res_simu_OLS")
if (!dir.exists(path_res_simu_OLS)) {
  dir.create(path_res_simu_OLS)
}  

#==============================================================================#
# Libraries and lib files ----------
#==============================================================================#

source(file.path(path_lib_code, "CI_mean.R"))
source(file.path(path_lib_code, "CI_OLS.R"))

library(purrr)
library(dplyr)

#==============================================================================#
# Temp to write and debug the function, and test ----------
#==============================================================================#

# n <- 1000
# 
# correlation_X1_X2 <- 0.5
# variance_X1 <- 1
# variance_X2 <- 2
# covariance_X1_X2 <- correlation_X1_X2 * sqrt(variance_X1) * sqrt(variance_X2)
# Sigma_matrix <- matrix(c(variance_X1, covariance_X1_X2, covariance_X1_X2, variance_X2), ncol = 2)
# 
# esperance_X1 <- 0
# esperance_X2 <- 0
# 
# 
# 
# beta0 <- c(2, 1, -3)
# 
# eps <- rnorm(n = n)
# 
# X <- drawn_X 
# 
# Y <- c(beta0[[1]] + drawn_X %*% beta0[2:3] + eps)
# 
# 
# # default argument
# alpha = 0.05
# a = NULL
# power_of_n_for_b = NULL
# omega = NULL
# power_of_n_for_omega = NULL
# bounds = list(lambda_reg = NULL,
#               K_reg = NULL,
#               K_eps = NULL,
#               K_xi = NULL,
#               C = NULL,
#               B = NULL)
# K_xi = NULL
# param_BE_EE = list(
#   choice = "best",
#   setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
#   regularity = list(C0 = 1, p = 2),
#   eps = 0.1)
# options = list(center = FALSE, bounded_case = FALSE, with_Exp_regime = FALSE)
# matrix_u = diag(NCOL(X) + 1)
# 
# 
# 
# i <- 3
# beta0[[2]] * drawn_X[i,1] + beta0[[3]] * drawn_X[i,2]
# (drawn_X %*% beta0[2:3])[i,]
# 
# 2 * -0.07 - 3 * 1.4
# 
# # importance of the choice for the power of b_n and omega_n to exit R regime
# 
# n <- 20 * 10^3
# # n <- 3 * 10^6
# 
# mydata <- Draw_data_regression_X_Normal_eps_Gumbel(n = n)
# 
# res <- Navae_ci_ols(Y = mydata$Y, X = mydata$X, 
#                     alpha = 0.05, K_xi = 9, verbose = 2,
#                     options = list(center = FALSE, bounded_case = FALSE, with_Exp_regime = FALSE))
# 
# res
# 
# 0.04136297
# 0.04136305

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
# Functions for simulations ----------------------------------------------------
#==============================================================================#

Coverage_of_CI <- function(ci, theta1)
{
  return((ci[[1]] <= theta1) && (theta1 <= ci[[2]]))
}

Length_of_CI <- function(ci)
{
  return(ci[[2]] - ci[[1]])
}

CI_OLS_one_draw_one_n <- function(simu_number, n)
{
  
  if (simu_number %% nb_simu_to_monitor == 0) {
    print(paste("Start of simulation number:", simu_number))
  }
  
  mydata <- Draw_data_regression_X_Normal_eps_Gumbel(
    n = n, 
    esperance_X1 = esperance_X1, esperance_X2 = esperance_X2,
    variance_X1 = variance_X1, variance_X2 = variance_X2, 
    correlation_X1_X2 = correlation_X1_X2, 
    beta0 = beta0)
  
  res_K9 <- Navae_ci_ols(Y = mydata$Y, X = mydata$X, 
                         alpha = alpha, K_xi = 9, verbose = 0, 
                         a = constant_for_b_n * n^power_of_n_for_b,
                         omega = n^power_of_n_for_omega,
                         options = options_navae_ci)
  
  res_Kh <- Navae_ci_ols(Y = mydata$Y, X = mydata$X, 
                         alpha = alpha, K_xi = NULL, verbose = 0,
                         a = constant_for_b_n * n^power_of_n_for_b,
                         omega = n^power_of_n_for_omega,
                         options = options_navae_ci)
  
  ci_navae_K9 <- res_K9$ci_navae[coeff_considered, ]
  ci_navae_Kh <- res_Kh$ci_navae[coeff_considered, ]
  ci_asymp <- res_K9$ci_asymp[coeff_considered, ]
  
  length_navae_K9 <- ci_navae_K9$length
  length_navae_Kh <- ci_navae_Kh$length
  length_asymp <- ci_asymp$length
  
  coverage_navae_K9 <- Coverage_of_CI(ci = ci_navae_K9, theta1 = beta0[[coeff_considered]])
  coverage_navae_Kh <- Coverage_of_CI(ci = ci_navae_Kh, theta1 = beta0[[coeff_considered]])
  coverage_asymp <- Coverage_of_CI(ci = ci_asymp, theta1 = beta0[[coeff_considered]])
  
  ratio_wrt_asymp_navae_K9 <- length_navae_K9 / length_asymp
  ratio_wrt_asymp_navae_Kh <- length_navae_Kh / length_asymp
  
  return(data.frame(simu_number = simu_number,
                    ci_navae_K9 = ci_navae_K9,
                    ci_navae_Kh = ci_navae_Kh,
                    ci_asymp = ci_asymp,
                    length_navae_K9 = length_navae_K9,
                    length_navae_Kh = length_navae_Kh,
                    length_asymp = length_asymp,
                    coverage_navae_K9 = coverage_navae_K9,
                    coverage_navae_Kh = coverage_navae_Kh,
                    coverage_asymp = coverage_asymp,
                    ratio_wrt_asymp_navae_K9 = ratio_wrt_asymp_navae_K9,
                    ratio_wrt_asymp_navae_Kh = ratio_wrt_asymp_navae_Kh))
}

CI_OLS_one_n <- function(n, S)
{
  
  list_res_across_simulations <- purrr::map(1:S, CI_OLS_one_draw_one_n, n = n)
  
  res_one_n <- purrr::reduce(list_res_across_simulations, rbind)
  
  rownames(res_one_n) <- NULL
  res_one_n$sample_size_n <- n
  res_one_n$S <- S
  
  return(res_one_n)
}

Get_filename_OLS <- function(n, S)
{
  return(paste0("res_simu_OLS",
                "_", "Gumbel_setting",
                "_", "coeff_", coeff_considered,
                "_", "alpha_", alpha,
                "_", "cst_bn_", constant_for_b_n,
                "_", "pw_bn_", power_of_n_for_b,
                "_", "pw_omegan_", power_of_n_for_omega,
                "_", "n_", n /10^3, "k",
                "_", "S_", S))
}

Perform_and_save_simu_across_n <- function(vec_n, S)
{
  
  for (n in vec_n) {
    
    cat(paste0("START OF SIMULATION FOR n = ", n))
    
    filename <- Get_filename_OLS(n = n, S = S)
    
    cat(paste0("\n", "Filename: ", filename, "\n"))
    
    res_one_n <- CI_OLS_one_n(n = n, S = S)
    
    saveRDS(object = res_one_n,
            file = file.path(path_res_simu_OLS, filename))
    
    rm(res_one_n)
  }
}

#==============================================================================#
# Global argument for simulations ----------------------------------------------
#==============================================================================#

alpha <- 0.05

coeff_considered <- 2

# Options for NAVAE CI
options_navae_ci <- list(center = FALSE, bounded_case = FALSE, with_Exp_regime = FALSE)

# constant_for_b_n <- 1
# power_of_n_for_b <- -1/5
constant_for_b_n <- 20
power_of_n_for_b <- -2/5

power_of_n_for_omega <- -1/5

# For the DGP
esperance_X1 = 0
esperance_X2 = 0
variance_X1 = 1
variance_X2 = 2
correlation_X1_X2 = 0.5
beta0 = c(2, 1, -3)

# Monitor simulations
nb_simu_to_monitor <- 10
S <- nb_MC_repetitions <- 100

# Choice of the grid of n
vec_n <- c(15, 25, 35, 50, 75, 100, 150) * 10^3

# vec_n <- c(99) * 10^3

#==============================================================================#
# Run and save the simulations -------------------------------------------------
#==============================================================================#

Perform_and_save_simu_across_n(vec_n = vec_n, S = S)

#==============================================================================#
# Coverage and graphs -------------------------------------------------
#==============================================================================#

# Load results over n
# vec_n_to_load <- c(15, 25, 35, 50, 75, 99, 100, 150) * 10^3
vec_n_to_load <- c(15, 25, 35, 50, 75, 100, 150) * 10^3
res_over_n <- NULL
for (n in vec_n_to_load) {
  res_one_n <- readRDS(file = file.path(path_res_simu_OLS,
                                        Get_filename_OLS(n = n, S = S)))
  if (is.null(res_over_n)) {
    res_over_n <- res_one_n
  } else {
    res_over_n <- rbind(res_over_n, res_one_n)
  }
}
rm(res_one_n)

colnames(res_over_n)

# Coverage

coverage <- res_over_n %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(coverage_navae_K9), mean(coverage_navae_Kh), mean(coverage_asymp))
coverage

# Proportion of R regime

prop_R_regime <- res_over_n %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(ci_navae_K9.regime != "Edg"),
                   mean(ci_navae_Kh.regime != "Edg"))
prop_R_regime

# Length

median_absolute_length <- res_over_n %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(median(length_navae_K9), median(length_navae_Kh), median(length_asymp))
median_absolute_length

mean_absolute_length <- res_over_n %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(length_navae_K9), mean(length_navae_Kh), mean(length_asymp))
mean_absolute_length

median_ratio_length <- res_over_n %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(median(ratio_wrt_asymp_navae_K9), median(ratio_wrt_asymp_navae_Kh))
median_ratio_length

mean_ratio_length <- res_over_n %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(ratio_wrt_asymp_navae_K9), mean(ratio_wrt_asymp_navae_Kh))
mean_ratio_length

# Just to check the behavior of the median with Inf, treated as a value itself, so can be Inf
# if at least half the values are Inf.
median(c(1,2,Inf,Inf))
median(c(1,2, 3,Inf,Inf))
