#' Analyze simulation results and make graphs
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

path_res_simu <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "res_simu")

#==============================================================================#
# Libraries and lib files ----------
#==============================================================================#

library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)

#==============================================================================#
# Temp tests ----------
#==============================================================================#

# Test for the sample size to exit R regime
n <- 100 * 10^3
alpha <- 0.1
res_ci_Kp <- Navae_ci_mean(data = draw_data(n), alpha = alpha,
                           a = NULL, power_of_n_for_b = power_of_n_for_b,
                           bound_K = NULL, param_BE_EE = param_BE_EE, verbose = TRUE)
res_ci_K9 <- Navae_ci_mean(data = draw_data(n), alpha = alpha,
                           a = NULL, power_of_n_for_b = power_of_n_for_b,
                           bound_K = 9, param_BE_EE = param_BE_EE, verbose = TRUE)
res_ci_Kp
res_ci_K9

#==============================================================================#
# Load results ----------
#==============================================================================#

# Load result
theta1 <- 1
power_of_n_for_b <- -1/5
M <- 1000
vec_n <- c(1, 5, 10, 50, 100, 200) * 10^3

res_over_n <- NULL
for (n in vec_n) {
  filename <- paste0("res_simu",
                     "_rexp_", theta1,
                     "_powerb_", power_of_n_for_b,
                     "_M_", M,
                     "_n_in_k_", n/10^3,
                     ".rds")
  
  res_for_n <- readRDS(file = file.path(path_res_simu, filename))
  
  if (is.null(res_over_n)) {
    res_over_n <- res_for_n
  } else {
    res_over_n <- rbind(res_over_n, res_for_n)
  }
}

colnames(res_over_n)

#==============================================================================#
# Table: minimal alpha to exit R regime ----------
#==============================================================================#

minimal_alpha_to_exit_R_regime_K9 <- 
  unique(res_over_n %>% 
           dplyr::filter(bound_K_method == "bound") %>%
           select(minimal_alpha_to_exit_R_regime))
minimal_alpha_to_exit_R_regime_K9_bis <- # just to check, but no variation across sample here.
  res_over_n %>% 
  dplyr::filter(bound_K_method == "bound") %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(minimal_alpha_to_exit_R_regime))
                   
minimal_alpha_to_exit_R_regime_Kp_mean <- 
  res_over_n %>% 
  dplyr::filter(bound_K_method == "plug-in") %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(minimal_alpha_to_exit_R_regime))
minimal_alpha_to_exit_R_regime_Kp_median <- 
  res_over_n %>% 
  dplyr::filter(bound_K_method == "plug-in") %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(median(minimal_alpha_to_exit_R_regime))
# Remark, it uses the rows for different alphas that, for a given simulation, use
# the same sample, so minimal_alpha_to_exit_R_regime are dupplicated, but not an issue
# as:
mean(c(1,2)) == mean(c(1, 1, 2, 2))
median(c(1,2)) == median(c(1, 1, 2, 2))

minimal_alpha_to_exit_R_regime_K9_bis
minimal_alpha_to_exit_R_regime_Kp_mean
minimal_alpha_to_exit_R_regime_Kp_median

minimal_alpha_to_exit_R_regime_K9_bis$`mean(minimal_alpha_to_exit_R_regime)`
minimal_alpha_to_exit_R_regime_Kp_mean$`mean(minimal_alpha_to_exit_R_regime)`
minimal_alpha_to_exit_R_regime_Kp_median$`median(minimal_alpha_to_exit_R_regime)`

#==============================================================================#
# Coverage ----------
#==============================================================================#

alpha_studied <- 0.1

bound_K_method_studied <- "plug-in"
# bound_K_method_studied <- "bound"

coverage <- 
  res_over_n %>% 
  dplyr::filter(bound_K_method == "plug-in") %>%
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(coverage), mean(coverage_CLT))
coverage

#==============================================================================#
# Graph lengths ----------
#==============================================================================#

alpha_studied <- 0.1
bound_K_method_studied <- "plug-in"
bound_K_method_studied <- "bound"

# Without changing the Inf into NA (to be removed from the average)
# the mean is thus Inf
res_over_n %>% 
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::filter(bound_K_method == bound_K_method_studied) %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(length_ci), mean(length_ci_CLT))

res_over_n_for_length <- res_over_n
res_over_n_for_length$length_ci[which(is.infinite(res_over_n_for_length$length_ci))] <- NA
res_over_n_for_length$relative_length <- 
  res_over_n_for_length$length_ci / res_over_n_for_length$length_ci_CLT

res_over_n_for_length %>% 
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::filter(bound_K_method == bound_K_method_studied) %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(length_ci, na.rm = TRUE), mean(length_ci_CLT))

# Longueur non distinguable à 3 ou 4 décimales entre K = 9 et K plug-in
# Logique car ici le vrai kurtosis est en effet 9 pour une exponentielle.
res_over_n %>%
  dplyr::filter(bound_K_method == "plug-in") %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(bound_K_value))

# Graph absolute length

df_absolute_length_K9 <-
  res_over_n_for_length %>% 
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::filter(bound_K_method == "bound") %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(length_ci, na.rm = TRUE), mean(length_ci_CLT))

colnames(df_absolute_length_K9) <- c("n", "Our", "CLT")

df_absolute_length_K9 <- 
  df_absolute_length_K9 %>% 
  pivot_longer(!n, names_to = "CI", 
               values_to = "Absolute_length")
  
g_abs_length_K9 <- ggplot(data = df_absolute_length_K9,
                          aes(x = n, y = Absolute_length, linetype = CI)) +
  geom_line() + 
  xlab("sample size") +
  ylab(paste0("CI length (mean over ", M, " repetitions)"))

g_abs_length_K9

# Graph relative length

df_relative_length_K9 <-
  res_over_n_for_length %>% 
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::filter(bound_K_method == "bound") %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(relative_length, na.rm = TRUE))
colnames(df_relative_length_K9) <- c("n", "Relative_length")

g_rel_length_K9 <- ggplot(data = df_relative_length_K9,
                          aes(x = n, y = Relative_length)) +
  geom_line() + 
  xlab("sample size") +
  ylab(paste0("Relative length (mean over ", M, " repetitions)")) +
  geom_hline(yintercept = 1, linetype = "dashed")

g_rel_length_K9
