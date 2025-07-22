#' Script for simulations NAVAE CI case of a mean/expectation
#' for the table of alpha min
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

path_res_simu_mean <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "res_simu_mean")
if (!dir.exists(path_res_simu_mean)) {dir.create(path_res_simu_mean)}  

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

power_of_n_for_b <- -1/5

param_BE_EE <- list(
  choice = "best",
  setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
  regularity = list(C0 = 1, p = 2),
  eps = 0.1)

M <- 20 * 10^3

vec_n <- c(0.5, 1, 5, 10, 50, 100, 200) * 10^3

# Functions

Draw_data <- function(n, theta1) 
{
  return(rexp(n, rate = 1 / theta1))
}

Get_filename_of_simu_results <- function(
    n, theta1, power_of_n_for_b, M)
{
  return(paste0("res_simu_alpha_min",
                "_rexp_", theta1,
                "_on_same_sample_k_un_var",
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

for (n in vec_n) {
  
  filename_res_simu <- Get_filename_of_simu_results(n, theta1, power_of_n_for_b, M)
  
  cat(paste("START OF SIMU FOR:", filename_res_simu, "\n"))
  print(Sys.time())
  
  # Arbitrary value for alpha, say, 0.42, just recover the alpha_min
  # With K = 9, the quantity is not data-dependent (only throung sample size n);
  # hence can be computed just once (outside of the loop over MC replications).
  
  drawn_data <- Draw_data(n, theta1)
  
  res_ci_K9_unknown_variance <- Navae_ci_mean(
    data = drawn_data, alpha = 0.42, a = 1 + n^(power_of_n_for_b), 
    bound_K = 9, known_variance = NULL, param_BE_EE = param_BE_EE)
  
  res_ci_K9_known_variance <- Navae_ci_mean(
    data = drawn_data, alpha = 0.42,
    bound_K = 9, known_variance = true_variance, param_BE_EE = param_BE_EE)

  res_M_MC_simu_for_given_n <- foreach(m = 1:M, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
    
    drawn_data <- Draw_data(n, theta1)
    
    res_ci_Kh_unknown_variance <- Navae_ci_mean(
      data = drawn_data, alpha = 0.42, a = 1 + n^(power_of_n_for_b), 
      bound_K = NULL, known_variance = NULL, param_BE_EE = param_BE_EE)
    
    df_res_for_one_draw <- data.frame(
      simu_number_m = m,
      alpha_min_Kh_unknown_variance = res_ci_Kh_unknown_variance$minimal_alpha_to_exit_R_regime)
    
    return(df_res_for_one_draw)
  }
  
  res_M_MC_simu_for_given_n$alpha_min_K9_unknown_variance = res_ci_K9_unknown_variance$minimal_alpha_to_exit_R_regime
  res_M_MC_simu_for_given_n$alpha_min_K9_known_variance = res_ci_K9_known_variance$minimal_alpha_to_exit_R_regime
  
  res_M_MC_simu_for_given_n$nb_simu_M <- M
  res_M_MC_simu_for_given_n$n <- n
  
  stopImplicitCluster()
  
  saveRDS(res_M_MC_simu_for_given_n,
          file = file.path(path_res_simu_mean, filename_res_simu))
  
}

stopCluster(cl)

toc()

#==============================================================================#
# Load results -----------------------------------------------------------------
#==============================================================================#

res_across_n <- NULL

for (n in vec_n) {
  res_for_one_n <- readRDS(
    file = file.path(path_res_simu_mean, 
                     Get_filename_of_simu_results(n, theta1, power_of_n_for_b, M)))
  if (is.null(res_across_n)) {
    res_across_n <- res_for_one_n
  } else {
    res_across_n <- rbind(res_across_n, res_for_one_n)
  }
  rm(res_for_one_n)
}

View(res_across_n)

#==============================================================================#
# Table alpha minimal to exit R regime -----------------------------------------
#==============================================================================#

minimal_alpha <- res_across_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(mean(alpha_min_Kh_unknown_variance), median(alpha_min_Kh_unknown_variance),
                   mean(alpha_min_K9_unknown_variance), mean(alpha_min_K9_known_variance))

minimal_alpha$`mean(alpha_min_Kh_unknown_variance)`

minimal_alpha$`median(alpha_min_Kh_unknown_variance)`
