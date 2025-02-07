#' Script for simulations NAVAE CI case of a mean/expectation
#' 
#' 
#' 

#==============================================================================#
# Paths ----------
#==============================================================================#

# Temp (after a package), path to the R code for the CI

# Select the correct path to main folder depending on user/setting
# Setting: CREST, Lucas's computer
# path_dropbox <- file.path("C:/Temp/Dropbox")
# Setting: personal computer (Lucas)
path_dropbox <- file.path("/Users/lucasgirard/Dropbox")
path_lib_code <- file.path(path_dropbox, "finite_CI_perso", "NavaeCI", "R")

path_res_simu <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "res_simu")
if (!dir.exists(path_res_simu)) {
  dir.create(path_res_simu)
}  

#==============================================================================#
# Libraries and lib files ----------
#==============================================================================#

source(file.path(path_lib_code, "CI_mean.R"))

library(purrr)
library(dplyr)

#==============================================================================#
# Functions for the simulations ----------
#==============================================================================#

Perform_one_simu_for_given_n <- function(n)
{
  
  # With an arbitrary alpha, with verbose = TRUE, to recover various informations
  # using either K = 9 or plug-in, that does not depend on the particular sample
  drawn_data <- draw_data(n) # so we just draw a sample here
  res_ci_K9 <- Navae_ci_mean(data = drawn_data, alpha = 0.42,
                             a = NULL, power_of_n_for_b = power_of_n_for_b,
                             bound_K = 9, param_BE_EE = param_BE_EE, verbose = TRUE)
  res_ci_Kp <- Navae_ci_mean(data = drawn_data, alpha = 0.42,
                             a = NULL, power_of_n_for_b = power_of_n_for_b,
                             bound_K = NULL, param_BE_EE = param_BE_EE, verbose = TRUE)
  
  # Recover the minimal alpha to exit R regime with K = 9 
  # That value is independent of data, of except sample size n of course
  alpha_to_exit_R_K9 <- res_ci_K9$minimal_alpha_to_exit_R_regime
  # Round and add a bit just to be able to compute ci not R it with K = 9 and have
  # some change to compute ci not R regime with plug-in K.
  min_alpha_used <- min(0.99, ceiling(100*alpha_to_exit_R_K9)/100)
  add_to_min_alpha_used <- 0.03
  
  if (min_alpha_used > 0.10) {
    vec_alpha_used <- c(min_alpha_used, min(0.99, min_alpha_used + add_to_min_alpha_used))
  } else if (min_alpha_used > 0.05) {
    vec_alpha_used <- c(min_alpha_used, 0.10)
  } else if (min_alpha_used > 0.01) {
    vec_alpha_used <- c(min_alpha_used, 0.05, 0.10)
  } else {
    vec_alpha_used <- c(min_alpha_used, 0.01, 0.05, 0.10)
  }
  
  res <- purrr::reduce(
    purrr::map(1:M, Perform_one_simu_for_given_n_across_alpha,
                        vec_alpha_used = vec_alpha_used, n = n), rbind)
  
  res$nb_simu_M <- M
  res$sample_size_n <- n
  
  return(res)
}

Perform_one_simu_for_given_n_across_alpha <- function(simu_number, vec_alpha_used, n)
{
  if (simu_number %% nb_simu_to_monitor == 0) {
    print(paste("Start of simulation number:", simu_number))
  }
  
  drawn_data <- draw_data(n) # NB: the same sample is used for the different alphas
  
  list_res_across_alpha <- purrr::map(vec_alpha_used, 
                                      Compute_ci_for_given_alpha_and_data, 
                                      drawn_data = drawn_data)
  res_across_alpha <- purrr::reduce(list_res_across_alpha, rbind)
  res_across_alpha$simu_number <- simu_number
  
  return(res_across_alpha)
}

Compute_ci_for_given_alpha_and_data <- function(alpha, drawn_data)
{
  
  res_ci_K9 <- Navae_ci_mean(data = drawn_data, alpha = alpha,
                             a = NULL, power_of_n_for_b = power_of_n_for_b,
                             bound_K = 9, param_BE_EE = param_BE_EE, verbose = TRUE)
  res_K9 <- data.frame(res_ci_K9)[1, ] # get only the first row (two rows since ci has two elements)
  res_K9$ci <- res_K9$ci_CLT <- NULL # suppress ci and ci_CLT, later column with lower and upper ends
  res_K9$alpha <- alpha
  res_K9$ci_lo <- res_ci_K9$ci[[1]]; res_K9$ci_up <- res_ci_K9$ci[[2]]
  res_K9$ci_CLT_lo <- res_ci_K9$ci_CLT[[1]]; res_K9$ci_CLT_up <- res_ci_K9$ci_CLT[[2]]
  res_K9$length_ci <- Length_of_CI(res_ci_K9$ci)
  res_K9$length_ci_CLT <- Length_of_CI(res_ci_K9$ci_CLT)
  res_K9$coverage <- Coverage_of_CI(res_ci_K9$ci, theta1)
  res_K9$coverage_CLT <- Coverage_of_CI(res_ci_K9$ci_CLT, theta1)
  
  res_ci_Kp <- Navae_ci_mean(data = drawn_data, alpha = alpha,
                             a = NULL, power_of_n_for_b = power_of_n_for_b,
                             bound_K = NULL, param_BE_EE = param_BE_EE, verbose = TRUE)
  res_Kp <- data.frame(res_ci_Kp)[1, ] # get only the first row (two rows since ci has two elements)
  res_Kp$ci <- res_Kp$ci_CLT <- NULL # suppress ci and ci_CLT, later column with lower and upper ends
  res_Kp$alpha <- alpha
  res_Kp$ci_lo <- res_ci_Kp$ci[[1]]; res_Kp$ci_up <- res_ci_Kp$ci[[2]]
  res_Kp$ci_CLT_lo <- res_ci_Kp$ci_CLT[[1]]; res_Kp$ci_CLT_up <- res_ci_Kp$ci_CLT[[2]]
  res_Kp$length_ci <- Length_of_CI(res_ci_Kp$ci)
  res_Kp$length_ci_CLT <- Length_of_CI(res_ci_Kp$ci_CLT)
  res_Kp$coverage <- Coverage_of_CI(res_ci_Kp$ci, theta1)
  res_Kp$coverage_CLT <- Coverage_of_CI(res_ci_Kp$ci_CLT, theta1)
  
  res <- rbind(res_K9, res_Kp)
  return(res)
}

Coverage_of_CI <- function(ci, theta1)
{
  return((ci[[1]] <= theta1) && (theta1 <= ci[[2]]))
}

Length_of_CI <- function(ci)
{
  return(ci[[2]] - ci[[1]])
}

#==============================================================================#
# Perform simulations ----------
#==============================================================================#

# Choice DGP - i.i.d Exponential with expectation theta1 (set to 1).
theta1 <- 1
draw_data <- function(n){rexp(n, rate = 1/theta1)}

# Choice hyper-parameters
power_of_n_for_b <- -1/5
param_BE_EE = list(
  choice = "best",
  setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
  regularity = list(C0 = 1, p = 2),
  eps = 0.1)

# Number of MC repetitions
M <- 1000

# To monitor evolution
nb_simu_to_monitor <- 100

# Choice of the different sample sizes
vec_n <- c(1, 5, 10, 50, 100, 200) * 10^3
# vec_n <- c(1, 5) * 10^3 # to test 

set.seed(1)

for(i_n in seq_along(vec_n)) {
  
  n <- vec_n[[i_n]]
  
  print(paste("Start of simu for n =", n, "-", Sys.time()))
  
  filename <- paste0("res_simu",
                     "_rexp_", theta1,
                     "_powerb_", power_of_n_for_b,
                     "_M_", M,
                     "_n_in_k_", n/10^3,
                     ".rds")
  
  res_for_n <- Perform_one_simu_for_given_n(n)
  
  saveRDS(object = res_for_n,
          file = file.path(path_res_simu, filename))
  
  rm(res_for_n)
}
