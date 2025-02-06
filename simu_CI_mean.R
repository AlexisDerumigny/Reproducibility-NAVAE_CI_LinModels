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
path_dropbox <- file.path("C:/Temp/Dropbox")
# Setting: personal computer (Lucas)
# path_dropbox <- file.path("/Users/lucasgirard/Dropbox")
path_lib_code <- file.path(path_dropbox, "finite_CI_perso", "NavaeCI", "R")

#==============================================================================#
# Libraries and lib files ----------
#==============================================================================#

source(file.path(path_lib_code, "CI_mean.R"))

library(foreach)
library(parallel)
library(doParallel)
library(tictoc)

#==============================================================================#
# XXXX ----------
#==============================================================================#

# Choice DGP - i.i.d Exponential with expectation theta1 (set to 1).
theta1 <- 1
draw_data <- function(n){rexp(n, rate = 1/theta1)}

# Choice hyper-parameters
power_of_n_for_b <- -2/5
param_BE_EE = list(
  choice = "best",
  setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
  regularity = list(C0 = 1, p = 2),
  eps = 0.1)

# Number of MC repetitions
M <- 5

#==============================================================================#
# Temp tests ----------
#==============================================================================#

# Test for the sample size to exit R regime
n <- 50 * 10^3
alpha <- 0.27
res_ci_Kp <- Navae_ci_mean(data = draw_data(n), alpha = alpha,
                           a = NULL, power_of_n_for_b = power_of_n_for_b,
                           bound_K = NULL, param_BE_EE = param_BE_EE, verbose = TRUE)
res_ci_K9 <- Navae_ci_mean(data = draw_data(n), alpha = alpha,
                           a = NULL, power_of_n_for_b = power_of_n_for_b,
                           bound_K = 9, param_BE_EE = param_BE_EE, verbose = TRUE)
res_ci_Kp
res_ci_K9



#==============================================================================#
# XXXX ----------
#==============================================================================#

# Choice of the different sample sizes
vec_n <- c(1, 5, 10, 50, 100, 200) * 10^3


tic()

# Parallelization setup
ncores <- max(1, min(detectCores() - 2, 30))
mc <- getOption("cl.cores", ncores)
cl <- makeCluster(mc)
registerDoParallel(cl)

# Seed for reproducibility
set.seed(1234567)


  

for (i_n in seq_along(vec_n)) {
  
  n <- vec_n[[i_n]]

  # Filename for simulations
  file_name <- paste0("df_thetas_stats_under_H0_",
                      "case_", case,
                      "_n_", n,
                      "_rhozw_", rho_treat_instr,
                      "_rhoev_", rho_treat_error,
                      "_M_", M,
                      ".rds")
  
  # To monitor progression of simulations
  print(paste("START OF SIMU FOR:", file_name))
  print(Sys.time())
  
  # Simulation and save of the results
  
  res_MC_simu_one_dgp_setting <- foreach(m = 1:M, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
    
  

Navae_ci_mean <- function(
  data,
  alpha = 0.05,
  a = NULL, power_of_n_for_b = NULL
  bound_K = NULL,
  param_BE_EE = list(
    choice = "best",
    setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
    regularity = list(C0 = 1, p = 2),
    eps = 0.1))





