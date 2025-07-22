#' Script for simulations NAVAE CI case of OLS
#' 
#' 

#==============================================================================#
# Paths ------------------------------------------------------------------------
#==============================================================================#

# Select the correct path to main folder depending on user/setting
# Setting: CREST, Lucas's computer
# path_dropbox <- file.path("C:/Temp/Dropbox")
# Setting: personal computer (Lucas)
path_dropbox <- file.path("/Users/lucasgirard/Dropbox")

# Temp (later, a package), path to the R code for the CI
path_lib_code <- file.path(path_dropbox, "finite_CI_perso", "NavaeCI", "R")

path_res_simu_OLS <- file.path(path_dropbox, "finite_CI_perso",
                               "Reproducibility-NAVAE_CI_LinModels",
                               "res_simu_OLS")
if (!dir.exists(path_res_simu_OLS)) {dir.create(path_res_simu_OLS)}

#==============================================================================#
# Libraries and lib files ------------------------------------------------------
#==============================================================================#

source(file.path(path_lib_code, "CI_mean.R"))
source(file.path(path_lib_code, "CI_OLS.R"))

library(purrr)
library(dplyr)
library(tictoc)
library(foreach)
library(parallel)
library(doParallel)
# Warning: issue (at least in my configuration) with paralelization inside a function
# called; hence the structure of the code

#==============================================================================#
# Functions to draw data and for simulations -----------------------------------
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
  Sigma_matrix <- matrix(c(variance_X1, covariance_X1_X2, covariance_X1_X2, variance_X2), 
                         ncol = 2)
  
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

Coverage_of_CI <- function(ci, theta1)
{
  return((ci[[1]] <= theta1) && (theta1 <= ci[[2]]))
}

Get_filename_OLS <- function(
    coeff_considered, alpha, constant_for_b_n, power_of_n_for_b, power_of_n_for_omega, n, M)
{
  return(paste0("res_simu_OLS",
                "_", "Gumbel_setting",
                "_", "coeff_", coeff_considered,
                "_", "alpha_", alpha,
                "_", "cst_bn_", constant_for_b_n,
                "_", "pw_bn_", power_of_n_for_b,
                "_", "pw_omegan_", power_of_n_for_omega,
                "_", "n_", n /10^3, "k",
                "_", "M_", M))
}

#==============================================================================#
# Settings for simulations -----------------------------------------------------
#==============================================================================#

# For nominal level
alpha <- 0.10

# For the coeff considered
coeff_considered <- 2

# Options for NAVAE CI
options_navae_ci <- list(center = FALSE, bounded_case = FALSE, with_Exp_regime = FALSE)

# Options for delta_n
param_BE_EE = list(
  choice = "best",
  setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
  regularity = list(C0 = 1, p = 2),
  eps = 0.1)

# constant_for_b_n <- 1
# power_of_n_for_b <- -1/5
# constant_for_b_n <- 100 # with 100 * n^{-2/5}, strange behavior in the simu for n = 100k, quite weird!
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

# Choice of the grid of n
# vec_n <- c(5, 10, 25, 50, 75, 100, 125, 150, 175, 200) * 10^3 # grid used for the mean
vec_n <- c(5, 10, 25, 50, 75, 100, 125, 150, 175, 200) * 10^3

# Number of Monte-Carlo repetitions
M <- 20000

# environ 1 minute pour M = 5, assez long donc malgré parallelisation.
# 36 minutes pour M = 200
# donc ferait quand même presque 25h pour M = 20,000 !
# Peut-être faire M = quelques milliers?
# 21 minutes quand M = 50.

# M = 2,000
# durant la nuit, sur prise secteur, 
16315 / 60^2 # un peu moins de 5h

#==============================================================================#
# Some tests for minimal sample size informative CI ----------------------------
#==============================================================================#

# constant_for_b_n <- 20
# power_of_n_for_b <- -2/5
# power_of_n_for_omega <- -1/5
# n <- 3656 # minimal n to start with K = 9, a_n = 1 + 20 * n^(-2/5) and omega_n = n^(-1/5)

# constant_for_b_n <- 1
# power_of_n_for_b <- -0.4
# 
# n <- 2000 * 10^3
# 
# # constant_for_b_n <- 1
# # power_of_n_for_b <- -1/5
# # power_of_n_for_omega <- -1/5
# # n <- 3660 # minimal n to start with K = 9, a_n = 1 + 20 * n^(-2/5) and omega_n = n^(-1/5)
# 
# drawn_data <- Draw_data_regression_X_Normal_eps_Gumbel(
#   n = n,
#   esperance_X1 = esperance_X1, esperance_X2 = esperance_X2,
#   variance_X1 = variance_X1, variance_X2 = variance_X2,
#   correlation_X1_X2 = correlation_X1_X2,
#   beta0 = beta0)
# 
# res_K9 <- Navae_ci_ols(Y = drawn_data$Y, X = drawn_data$X,
#                        alpha = alpha, K_xi = 9, verbose = 2,
#                        a = 1 + constant_for_b_n * n^power_of_n_for_b,
#                        omega = n^power_of_n_for_omega,
#                        param_BE_EE = param_BE_EE,
#                        options = options_navae_ci)
# 
# res_K9$ci_navae
# 
# res_Kh <- Navae_ci_ols(Y = drawn_data$Y, X = drawn_data$X,
#                        alpha = alpha, K_xi = NULL, verbose = 2,
#                        a = 1 + constant_for_b_n * n^power_of_n_for_b,
#                        omega = n^power_of_n_for_omega,
#                        param_BE_EE = param_BE_EE,
#                        options = options_navae_ci)
# 
# 
# res_Kh$ci_navae

#==============================================================================#
# Perform simulations for different sample sizes n -----------------------------
#==============================================================================#

tic()

# Settings for parallelization
ncores <- max(1, min(detectCores() - 2, 32))
mc <- getOption("cl.cores", ncores)
cl <- makeCluster(mc)
registerDoParallel(cl)

set.seed(13)

for (n in vec_n) {
  
  filename_res_simu <- Get_filename_OLS(
    coeff_considered, alpha, constant_for_b_n, power_of_n_for_b, power_of_n_for_omega, n, M)
  
  cat(paste("START OF SIMU FOR:", filename_res_simu, "\n"))
  print(Sys.time())
  
  res_M_MC_simu_for_given_n_and_alpha <- foreach(m = 1:M, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
    
    drawn_data <- Draw_data_regression_X_Normal_eps_Gumbel(
      n = n, 
      esperance_X1 = esperance_X1, esperance_X2 = esperance_X2,
      variance_X1 = variance_X1, variance_X2 = variance_X2, 
      correlation_X1_X2 = correlation_X1_X2, 
      beta0 = beta0)
    
    res_K9 <- Navae_ci_ols(Y = drawn_data$Y, X = drawn_data$X, 
                           alpha = alpha, K_xi = 9, verbose = 0, 
                           a = 1 + constant_for_b_n * n^power_of_n_for_b,
                           omega = n^power_of_n_for_omega, 
                           param_BE_EE = param_BE_EE,
                           options = options_navae_ci)
    
    res_Kh <- Navae_ci_ols(Y = drawn_data$Y, X = drawn_data$X, 
                           alpha = alpha, K_xi = NULL, verbose = 0,
                           a = 1 + constant_for_b_n * n^power_of_n_for_b,
                           omega = n^power_of_n_for_omega,
                           param_BE_EE = param_BE_EE,
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
    
    ratio_wrt_asymp_navae_K9 <- res_K9$ratio_length_wrt_ci_asymp[[coeff_considered]]
    ratio_wrt_asymp_navae_Kh <- res_Kh$ratio_length_wrt_ci_asymp[[coeff_considered]]
    
    df_res_for_one_draw <- data.frame(ci_navae_K9 = ci_navae_K9,
                                      ci_navae_Kh = ci_navae_Kh,
                                      ci_asymp = ci_asymp,
                                      length_navae_K9 = length_navae_K9,
                                      length_navae_Kh = length_navae_Kh,
                                      length_asymp = length_asymp,
                                      coverage_navae_K9 = coverage_navae_K9,
                                      coverage_navae_Kh = coverage_navae_Kh,
                                      coverage_asymp = coverage_asymp,
                                      ratio_wrt_asymp_navae_K9 = ratio_wrt_asymp_navae_K9,
                                      ratio_wrt_asymp_navae_Kh = ratio_wrt_asymp_navae_Kh,
                                      simu_number_m = m)
    return(df_res_for_one_draw)
  }
  
  res_M_MC_simu_for_given_n_and_alpha$nb_simu_M <- M
  res_M_MC_simu_for_given_n_and_alpha$n <- n
  res_M_MC_simu_for_given_n_and_alpha$alpha <- alpha
  
  stopImplicitCluster()
  
  saveRDS(res_M_MC_simu_for_given_n_and_alpha,
          file = file.path(path_res_simu_OLS, filename_res_simu))
  
}

stopCluster(cl)

toc()

#==============================================================================#
# Load results to check --------------------------------------------------------
#==============================================================================#

res_across_n <- NULL

for (n in vec_n) {
  res_for_one_n <- readRDS(
    file = file.path(path_res_simu_OLS, 
                     Get_filename_OLS(
                       coeff_considered, alpha, constant_for_b_n, power_of_n_for_b, power_of_n_for_omega, n, M)))
  if (is.null(res_across_n)) {
    res_across_n <- res_for_one_n
  } else {
    res_across_n <- rbind(res_across_n, res_for_one_n)
  }
  rm(res_for_one_n)
}

res_across_n
