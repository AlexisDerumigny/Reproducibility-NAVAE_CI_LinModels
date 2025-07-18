#' Script for simulations NAVAE CI case of a mean/expectation
#' for coverage table and length figures.
#' Compared to the first version, perform all the cases that we will finally
#' present: computed with K = 9 or K hat plug-in to have both.
#' known variance with b_n = n^{-1/5}, 
#' unknown variance with b_n = n^{-1/5},
#' known variance with optimized a_n (equivalent b_n),
#' and the CLT-based standard asymptotic CI.
#' For a given MC replication, all are computed on the same sample.

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

path_res_simu <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "res_simu_mean")
if (!dir.exists(path_res_simu)) {dir.create(path_res_simu)}  

#==============================================================================#
# Libraries and lib files ------------------------------------------------------
#==============================================================================#

source(file.path(path_lib_code, "CI_mean.R"))

library(purrr)
library(dplyr)
library(tictoc)
library(foreach)
library(parallel)
library(doParallel)
# Warning: issue (at least in my configuration) with paralelization inside a function
# called; hence the structure of the code

#==============================================================================#
# Functions and parameters for the simulations ---------------------------------
#==============================================================================#

# Parameters / Settings

theta1 <- 1
true_variance <- 1

alpha <- 0.10

power_of_n_for_b <- -1/5

param_BE_EE <- list(
  choice = "best",
  setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
  regularity = list(C0 = 1, p = 2),
  eps = 0.1)

# Functions

Draw_data <- function(n, theta1) 
{
  return(rexp(n, rate = 1 / theta1))
}

Length_of_CI <- function(ci)
{
  return(ci[[2]] - ci[[1]])
}

Coverage_of_CI <- function(ci, theta1)
{
  return((ci[[1]] <= theta1) && (theta1 <= ci[[2]]))
}

Get_filename_of_simu_results <- function(
    n, theta1, alpha, power_of_n_for_b, M)
{
  return(paste0("res_simu",
                "_rexp_", theta1,
                "_alpha_", alpha,
                "_on_same_sample_k_un_var_w_wo_optim_a",
                "_powerb_", power_of_n_for_b,
                "_M_", M,
                "_n_in_k_", n / 10^3,
                ".rds"))
}

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

M <- 20 * 10^3

vec_n <- c(5, 10, 25, 50, 75, 100, 125, 150, 175, 200) * 10^3

for (n in vec_n) {
  
  filename_res_simu <- Get_filename_of_simu_results(n, theta1, alpha, power_of_n_for_b, M)
  
  cat(paste("START OF SIMU FOR:", filename_res_simu, "\n"))
  print(Sys.time())
  
  res_M_MC_simu_for_given_n_and_alpha <- foreach(m = 1:M, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
    
    drawn_data <- Draw_data(n, theta1)
    
    # Known variance - K = 9 or plug-in (Kh = K hat) (no parameter a with known variance)
    # NB: it is the ratio with respect to CLT computed with known variance in this case.
    
    res_ci_K9_known_variance <- Navae_ci_mean(
      data = drawn_data, alpha = alpha, 
      bound_K = 9, known_variance = true_variance, verbose = TRUE, param_BE_EE = param_BE_EE)
    df_row_res_ci_K9_known_variance <- data.frame(
      method_ci = "navae",
      power_of_n_for_b = NA, optimize_in_a = NA, bound_K = 9, known_variance = TRUE,
      coverage = Coverage_of_CI(res_ci_K9_known_variance$ci, theta1),
      length = res_ci_K9_known_variance$length_ci,
      ratio_length_wrt_CLT = res_ci_K9_known_variance$ratio_length_wrt_CLT)
    
    res_ci_Kh_known_variance <- Navae_ci_mean(
      data = drawn_data, alpha = alpha, 
      bound_K = NULL, known_variance = true_variance, verbose = TRUE, param_BE_EE = param_BE_EE)
    df_row_res_ci_Kh_known_variance <- data.frame(
      method_ci = "navae",
      power_of_n_for_b = NA, optimize_in_a = NA, bound_K = "plug-in", known_variance = TRUE,
      coverage = Coverage_of_CI(res_ci_Kh_known_variance$ci, theta1),
      length = res_ci_Kh_known_variance$length_ci,
      ratio_length_wrt_CLT = res_ci_Kh_known_variance$ratio_length_wrt_CLT)
    
    df_row_res_ci_CLT_known_variance <- data.frame(
      method_ci = "clt",
      power_of_n_for_b = NA, optimize_in_a = NA, bound_K = NA, known_variance = TRUE,
      coverage = Coverage_of_CI(res_ci_K9_known_variance$ci_CLT, theta1),
      length = Length_of_CI(res_ci_K9_known_variance$ci_CLT),
      ratio_length_wrt_CLT = 1)
    
    # Unknown variance => tuning parameter a: optimized or not in a, K = 9 or plug-in (Kh = K hat)
    # NB: it is the ratio with respect to CLT computed with unknown variance in this case.
    
    res_ci_K9_unknown_variance_not_optimized_in_a <- Navae_ci_mean(
      data = drawn_data, alpha = alpha, a = NULL, power_of_n_for_b = power_of_n_for_b,
      optimize_in_a = FALSE, bound_K = 9, known_variance = NULL, verbose = TRUE, param_BE_EE = param_BE_EE)
    df_row_res_ci_K9_unknown_variance_not_optimized_in_a <- data.frame(
      method_ci = "navae",
      power_of_n_for_b = power_of_n_for_b, optimize_in_a = FALSE, bound_K = 9, known_variance = FALSE,
      coverage = Coverage_of_CI(res_ci_K9_unknown_variance_not_optimized_in_a$ci, theta1),
      length = res_ci_K9_unknown_variance_not_optimized_in_a$length_ci,
      ratio_length_wrt_CLT = res_ci_K9_unknown_variance_not_optimized_in_a$ratio_length_wrt_CLT)
    
    res_ci_Kh_unknown_variance_not_optimized_in_a <- Navae_ci_mean(
      data = drawn_data, alpha = alpha, a = NULL, power_of_n_for_b = power_of_n_for_b,
      optimize_in_a = FALSE, bound_K = NULL, known_variance = NULL, verbose = TRUE, param_BE_EE = param_BE_EE)
    df_row_res_ci_Kh_unknown_variance_not_optimized_in_a <- data.frame(
      method_ci = "navae",
      power_of_n_for_b = power_of_n_for_b, optimize_in_a = FALSE, bound_K = "plug-in", known_variance = FALSE,
      coverage = Coverage_of_CI(res_ci_Kh_unknown_variance_not_optimized_in_a$ci, theta1),
      length = res_ci_Kh_unknown_variance_not_optimized_in_a$length_ci,
      ratio_length_wrt_CLT = res_ci_Kh_unknown_variance_not_optimized_in_a$ratio_length_wrt_CLT)
    
    res_ci_K9_unknown_variance_optimized_in_a <- Navae_ci_mean(
      data = drawn_data, alpha = alpha, a = NULL, power_of_n_for_b = power_of_n_for_b,
      optimize_in_a = TRUE, bound_K = 9, known_variance = NULL, verbose = TRUE, param_BE_EE = param_BE_EE)
    df_row_res_ci_K9_unknown_variance_optimized_in_a <- data.frame(
      method_ci = "navae",
      power_of_n_for_b = power_of_n_for_b, optimize_in_a = TRUE, bound_K = 9, known_variance = FALSE,
      coverage = Coverage_of_CI(res_ci_K9_unknown_variance_optimized_in_a$ci, theta1),
      length = res_ci_K9_unknown_variance_optimized_in_a$length_ci,
      ratio_length_wrt_CLT = res_ci_K9_unknown_variance_optimized_in_a$ratio_length_wrt_CLT)
    
    res_ci_Kh_unknown_variance_optimized_in_a <- Navae_ci_mean(
      data = drawn_data, alpha = alpha, a = NULL, power_of_n_for_b = power_of_n_for_b,
      optimize_in_a = TRUE, bound_K = NULL, known_variance = NULL, verbose = TRUE, param_BE_EE = param_BE_EE)
    df_row_res_ci_Kh_unknown_variance_optimized_in_a <- data.frame(
      method_ci = "navae",
      power_of_n_for_b = power_of_n_for_b, optimize_in_a = TRUE, bound_K = "plug-in", known_variance = FALSE,
      coverage = Coverage_of_CI(res_ci_Kh_unknown_variance_optimized_in_a$ci, theta1),
      length = res_ci_Kh_unknown_variance_optimized_in_a$length_ci,
      ratio_length_wrt_CLT = res_ci_Kh_unknown_variance_optimized_in_a$ratio_length_wrt_CLT)
    
    df_row_res_ci_CLT_unknown_variance <- data.frame(
      method_ci = "clt",
      power_of_n_for_b = NA, optimize_in_a = NA, bound_K = NA, known_variance = FALSE,
      coverage = Coverage_of_CI(res_ci_Kh_unknown_variance_optimized_in_a$ci_CLT, theta1),
      length = Length_of_CI(res_ci_Kh_unknown_variance_optimized_in_a$ci_CLT),
      ratio_length_wrt_CLT = 1)
    
    df_res_for_one_draw <- rbind(df_row_res_ci_K9_known_variance,
                                 df_row_res_ci_Kh_known_variance,
                                 df_row_res_ci_K9_unknown_variance_not_optimized_in_a,
                                 df_row_res_ci_Kh_unknown_variance_not_optimized_in_a,
                                 df_row_res_ci_K9_unknown_variance_optimized_in_a,
                                 df_row_res_ci_Kh_unknown_variance_optimized_in_a,
                                 df_row_res_ci_CLT_known_variance,
                                 df_row_res_ci_CLT_unknown_variance)
    
    df_res_for_one_draw$simu_number_m <- m
    
    return(df_res_for_one_draw)
  }
  res_M_MC_simu_for_given_n_and_alpha$nb_simu_M <- M
  res_M_MC_simu_for_given_n_and_alpha$n <- n
  res_M_MC_simu_for_given_n_and_alpha$alpha <- alpha
  
  stopImplicitCluster()
  
  saveRDS(res_M_MC_simu_for_given_n_and_alpha,
          file = file.path(path_res_simu, filename_res_simu))
  
}

stopCluster(cl)

toc()

#==============================================================================#
# Load results to check --------------------------------------------------------
#==============================================================================#

res_across_n <- NULL

for (n in vec_n) {
  res_for_one_n <- readRDS(
    file = file.path(path_res_simu, 
                     Get_filename_of_simu_results(n, theta1, alpha, power_of_n_for_b, M)))
  if (is.null(res_across_n)) {
    res_across_n <- res_for_one_n
  } else {
    res_across_n <- rbind(res_across_n, res_for_one_n)
  }
  rm(res_for_one_n)
}

res_across_n
