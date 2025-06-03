#' Analyze simulation results and make graphs
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
path_res_simu <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "res_simu")
path_figures <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "figures")

#==============================================================================#
# Libraries and lib files ----------
#==============================================================================#

library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

source(file.path(path_lib_code, "CI_mean.R"))

#==============================================================================#
# Settings for exporting figures in Latex --------------------------------------
#==============================================================================#

size_theme_element_text <- 14

width_export_pdf_file <- 6
height_export_pdf_file <- 4.5

#==============================================================================#
# Temp tests ----------
#==============================================================================#

# Test for the sample size to exit R regime

alpha <- 0.1

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

n <- 0.5 * 10^3

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
M <- 5000 # 1000, 5000
# vec_n <- c(0.5, 1, 5, 10, 50, 100, 200) * 10^3
vec_n <- c(0.5, 1, 5, 10, 50) * 10^3

known_variance = NULL

res_over_n <- NULL
for (n in vec_n) {
  filename <- Get_filename(
    optimize_in_a, theta1, power_of_n_for_b, M, n)
  
  res_for_n <- readRDS(file = file.path(path_res_simu, filename))
  
  if (is.null(res_over_n)) {
    res_over_n <- res_for_n
  } else {
    res_over_n <- rbind(res_over_n, res_for_n)
  }
}
rm(res_for_n)
colnames(res_over_n)

View(res_over_n)

#==============================================================================#
# For both w and wo optimization in a #
#==============================================================================#

# Separate the two cases and choose with or without optimization
res_over_n_w_and_wo_optim_in_a <- res_over_n

res_over_n <- res_over_n_w_and_wo_optim_in_a %>%
  dplyr::filter(optimized_a_is_used == TRUE)

res_over_n <- res_over_n_w_and_wo_optim_in_a %>%
  dplyr::filter(optimized_a_is_used == FALSE)


#==============================================================================#
# Table: minimal alpha to exit R regime ----------------------------------------
#==============================================================================#

#' Select with or without optimization.
#' Rather without optimization is shown in the paper since the optimized a
#' tries to minimize the width, which contradicts the objective of exiting 
#' the R regime to be informative.

optimized_a_is_used_studied <- FALSE
res_over_n <- res_over_n_w_and_wo_optim_in_a %>%
  dplyr::filter(optimized_a_is_used == optimized_a_is_used_studied)

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

# Summary table
table_minimal_alpha_to_exit_R_regime <- 
  minimal_alpha_to_exit_R_regime_K9_bis %>%
  dplyr::rename(with_K_eq_9 = `mean(minimal_alpha_to_exit_R_regime)`) %>%
  dplyr::mutate(with_plug_in_mean = 
                  minimal_alpha_to_exit_R_regime_Kp_mean$`mean(minimal_alpha_to_exit_R_regime)`,
                with_plug_in_median =
                  minimal_alpha_to_exit_R_regime_Kp_median$`median(minimal_alpha_to_exit_R_regime)`)

# Compute the mininal alpha to exit R regime using optimized a and bound K = 9
# Table minimal alpha (one table only, not automated export to tex)
table_minimal_alpha_to_exit_R_regime


#==============================================================================#
# Table: coverages -------------------------------------------------------------
#==============================================================================#

alpha_studied <- 0.1

#' Get coverage of the CLT-CI, it is the same CLT-CI whatever
#' the bound used: plug-in or bound, and whatever the optimization in a
#' Therefore, we just select arbitrarily a combination to get the CLT-CI coverage.
coverage_CLT <- res_over_n_w_and_wo_optim_in_a %>%
  dplyr::filter(bound_K_method == "plug-in") %>%
  dplyr::filter(optimized_a_is_used == FALSE) %>%
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(coverage_CLT))

# bound_K_method_studied <- "plug-in"
bound_K_method_studied <- "bound"

coverage_wo_optim <- res_over_n_w_and_wo_optim_in_a %>%
  dplyr::filter(bound_K_method == bound_K_method_studied) %>%
  dplyr::filter(optimized_a_is_used == FALSE) %>%
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(coverage))

coverage_wi_optim <- res_over_n_w_and_wo_optim_in_a %>%
  dplyr::filter(bound_K_method == bound_K_method_studied) %>%
  dplyr::filter(optimized_a_is_used == TRUE) %>%
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::group_by(sample_size_n) %>%
  dplyr::summarise(mean(coverage))

# Full table
table_coverage <- coverage_CLT %>%
  dplyr::mutate(cov_wo_optim = coverage_wo_optim$`mean(coverage)`,
                cov_wi_optim = coverage_wi_optim$`mean(coverage)`) %>%
  dplyr::rename(cov_clt = `mean(coverage_CLT)`)

# To copy-paste in Latex (just one table to do... so not automated...)
nb_digits = 3
paste0(round(table_coverage$cov_clt, digits = nb_digits), collapse = " & ")
paste0(round(table_coverage$cov_wo_optim, digits = nb_digits), collapse = " & ")
paste0(round(table_coverage$cov_wi_optim, digits = nb_digits), collapse = " & ")

#==============================================================================#
# Graph; absolute lengths with or wo optimization in a --------------------------
#==============================================================================#

alpha_studied <- 0.1
# bound_K_method_studied <- "plug-in"
bound_K_method_studied <- "bound"

res_over_n <- res_over_n_w_and_wo_optim_in_a %>%
  dplyr::filter(optimized_a_is_used == FALSE)

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


g_abs_length_K9_for_latex <- ggplot(data = df_absolute_length_K9,
                                    aes(x = n, y = Absolute_length, linetype = CI)) +
  geom_line() + xlab(element_blank()) + ylab(element_blank()) +
  scale_linetype(guide = 'none') +
  theme(axis.text = element_text(size = size_theme_element_text)) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k"))
g_abs_length_K9_for_latex

figure_filename <- paste0("fig_absolute_length",
                          "_alpha_", alpha_studied,
                          "_K_", bound_K_method_studied,
                          "_b_", power_of_n_for_b,
                          "_M_", M,
                          ".pdf")
pdf(file = file.path(path_figures, figure_filename), 
    width = width_export_pdf_file, height = height_export_pdf_file)
print(g_abs_length_K9_for_latex)
dev.off()

#==============================================================================#
# Graph: absolute lengths with and wo optimization in a --------------------------
#==============================================================================#

alpha_studied <- 0.1

# bound_K_method_studied <- "plug-in"
bound_K_method_studied <- "bound"

# Put Inf values as NA
res_over_n_for_length <- res_over_n_w_and_wo_optim_in_a
res_over_n_for_length$length_ci[which(is.infinite(res_over_n_for_length$length_ci))] <- NA

res_over_n_for_length %>%
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::filter(bound_K_method == "bound")


# Create the dataframe to be used for the plot
df_absolute_length_K9 <-
  res_over_n_for_length %>% 
  dplyr::filter(alpha == alpha_studied) %>%
  dplyr::filter(bound_K_method == "bound") %>%
  dplyr::group_by(sample_size_n, optimized_a_is_used) %>%
  dplyr::summarise(mean(length_ci, na.rm = TRUE), mean(length_ci_CLT))
colnames(df_absolute_length_K9) <- c("n", "optimized_a", "Our", "CLT")
df_absolute_length_K9 <- 
  df_absolute_length_K9 %>% 
  pivot_wider(names_from = optimized_a, values_from = c("Our", "CLT")) %>%
  select(-CLT_TRUE) %>% # same values for CLT with or without 
  rename(Our_wi_optim = Our_TRUE,
         Our_wo_optim = Our_FALSE,
         CLT = CLT_FALSE) %>%
  pivot_longer(!n, names_to = "CI", 
               values_to = "Absolute_length") 

g_abs_length_K9 <- ggplot(data = df_absolute_length_K9,
                          aes(x = n, y = Absolute_length, linetype = CI)) +
  geom_line() + 
  xlab("sample size") +
  ylab(paste0("CI length (mean over ", M, " repetitions)"))

g_abs_length_K9

g_abs_length_K9_for_latex <- ggplot(data = df_absolute_length_K9,
                                    aes(x = n, y = Absolute_length, 
                                        color = CI, linetype = CI)) +
  geom_line(alpha = 0.9, linewidth = 1) + xlab(element_blank()) + ylab(element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = size_theme_element_text)) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k"))
g_abs_length_K9_for_latex

figure_filename <- paste0("fig_absolute_length",
                          "_wi_and_wo_optim_in_a",
                          "_alpha_", alpha_studied,
                          "_K_", bound_K_method_studied,
                          "_b_", power_of_n_for_b,
                          "_M_", M,
                          ".pdf")
pdf(file = file.path(path_figures, figure_filename), 
    width = width_export_pdf_file, height = height_export_pdf_file)
print(g_abs_length_K9_for_latex)
dev.off()

#==============================================================================#
# Graph relative length, fixed bound K -----------------------------------------
#==============================================================================#

# Note: when K and a_n are fixed independent of the data, the ratio of length
# is deterministic (absolute lenghts are, on the contrary, data-dependent
# through the empirical variance, but that term cancels out in the ratio).

# Wrapper function to recover ratio of length
Get_relative_length_when_data_independent <- function(
    n, alpha = 0.05, a = NULL, power_of_n_for_b = NULL, bound_K = NULL,
    param_BE_EE = list(
      choice = "best",
      setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
      regularity = list(C0 = 1, p = 2),
      eps = 0.1))
{
  res_navae_ci_mean <- Navae_ci_mean(
    data = rexp(n = n), # Fake data, just to use the function and with size n
    alpha = alpha, a = a, power_of_n_for_b = power_of_n_for_b,
    bound_K = bound_K, param_BE_EE = param_BE_EE, verbose = TRUE)
  
  return(res_navae_ci_mean$ratio_length_wrt_CLT)
}

alpha_studied <- 0.1
power_of_n_for_b <- -1/5
bound_K <- 9

Get_relative_length_when_data_independent_for_sapply <- function(n)
{
  Get_relative_length_when_data_independent(
    n, alpha = alpha_studied, a = NULL, power_of_n_for_b = power_of_n_for_b, 
    bound_K = bound_K,
    param_BE_EE = list(
      choice = "best",
      setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
      regularity = list(C0 = 1, p = 2),
      eps = 0.1))
}

grid_n <- seq(from = 5000, to = 200000, by = 100)
values <- sapply(grid_n, FUN = Get_relative_length_when_data_independent_for_sapply)
# Silly in terms of computing in the sense that we need only n, not the sample but remains a few seconds!
df_relative_length_K9 <- data.frame(
  n = grid_n, relative_length = values)

g_rel_length_K9_for_latex <- ggplot(data = df_relative_length_K9,
                                    aes(x = n, y = relative_length)) +
  geom_line() + xlab(element_blank()) + ylab(element_blank()) +
  geom_hline(yintercept = 1, linetype = "blank") +
  theme(axis.text=element_text(size = size_theme_element_text))  +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k"))
g_rel_length_K9_for_latex

figure_filename <- paste0("fig_relative_length",
                          "_alpha_", alpha_studied,
                          "_K_", bound_K_method_studied,
                          "_b_", power_of_n_for_b,
                          ".pdf")
pdf(file = file.path(path_figures, figure_filename), 
    width = width_export_pdf_file, height = height_export_pdf_file)
print(g_rel_length_K9_for_latex)
dev.off()



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
