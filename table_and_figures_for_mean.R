#' Table of coverage and graph of widths - case of expectation
#' This script:
#' Computes coverage
#' Outputs the two graphs of the width.
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

path_res_simu <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "res_simu_mean")
path_figures <- file.path(path_dropbox, "finite_CI_perso", "Reproducibility-NAVAE_CI_LinModels", "figures")

#==============================================================================#
# Libraries and lib files ------------------------------------------------------
#==============================================================================#

source(file.path(path_lib_code, "CI_mean.R"))

library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

#==============================================================================#
# Auxiliary functions ----------------------------------------------------------
#==============================================================================#

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
# Load simulations results -----------------------------------------------------
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

M <- 20 * 10^3

vec_n <- c(5, 10, 25, 50, 75, 100, 125, 150, 175, 200) * 10^3

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

#==============================================================================#
# Table Coverage ---------------------------------------------------------------
#==============================================================================#

View(res_across_n)
colnames(res_across_n)

coverage <- res_across_n %>%
  dplyr::group_by(n, method_ci, power_of_n_for_b, optimize_in_a, bound_K, known_variance) %>%
  dplyr::summarise(mean(coverage))
View(coverage)

vec_n_in_table <- c(5, 10, 50, 100, 200) * 10^3

(coverage %>% dplyr::filter(method_ci == "navae",
                            known_variance == TRUE,
                            bound_K == 9,
                            n %in% vec_n_in_table))$`mean(coverage)`

(coverage %>% dplyr::filter(method_ci == "navae",
                            known_variance == FALSE,
                            bound_K == 9,
                            optimize_in_a == FALSE,
                            n %in% vec_n_in_table))$`mean(coverage)`

(coverage %>% dplyr::filter(method_ci == "navae",
                            known_variance == FALSE,
                            bound_K == 9,
                            optimize_in_a == TRUE,
                            n %in% vec_n_in_table))$`mean(coverage)`

(coverage %>% dplyr::filter(method_ci == "clt",
                            known_variance == FALSE,
                            n %in% vec_n_in_table))$`mean(coverage)`

#==============================================================================#
# Function to get denomination of the different ci for the graphs --------------
#==============================================================================#

# Create a variable to indicate the method for the linetype and color
Denomination_ci <- function(
    method_ci, optimize_in_a, bound_K, known_variance) {
  
  denomination <- vector(mode = "character", length = length(method_ci))
  
  for (i in seq_along(denomination)) {
    
    if (method_ci[[i]] == "clt") {
      
      if (isTRUE(known_variance[[i]])) {
        denomination[[i]] <- "clt_with_known_variance"
      } else {
        denomination[[i]] <- "clt_with_unknown_variance"
      }
    }
    
    if (method_ci[[i]] == "navae") {
      
      if (isTRUE(known_variance[[i]])) {
        
        denomination[[i]] <- paste0("navae_with_known_variance_with_K_", bound_K[[i]])
        
      } else {
        
        if (isTRUE(optimize_in_a[[i]])) {
          
          denomination[[i]] <- paste0("navae_with_unknown_variance_optimized_in_a_with_K_", bound_K[[i]])
          
        } else {
          
          denomination[[i]] <- paste0("navae_with_unknown_variance_not_optimized_in_a_with_K_", bound_K[[i]])
        }
        
      }
    }
  }
  return(denomination)
}

#==============================================================================#
# Settings for exporting figures in Latex --------------------------------------
#==============================================================================#

size_theme_element_text <- 15
width_export_pdf_file <- 6
height_export_pdf_file <- 4.5

hue <- scales::hue_pal()(3)
my_red <- hue[[1]]
my_green <- hue[[2]]
my_blue <- hue[[3]]

linewidth_used <- 1.1

# Choice colors
my_colors <- c("clt_with_unknown_variance" = my_green,
               "navae_with_unknown_variance_optimized_in_a_with_K_9" = my_red,
               "navae_with_known_variance_with_K_9" = my_blue)

my_linetypes <- c("clt_with_unknown_variance" = "dotted",
                  "navae_with_unknown_variance_optimized_in_a_with_K_9" = "solid",
                  "navae_with_known_variance_with_K_9" = "dashed")

#==============================================================================#
# Graph absolute width ---------------------------------------------------------
#==============================================================================#

average_absolute_length <- res_across_n %>%
  dplyr::group_by(n, method_ci, power_of_n_for_b, optimize_in_a, bound_K, known_variance) %>%
  dplyr::summarise(mean(length))

# Choice to use plug-in or K = 9. In the paper, we show K = 9
df_for_graph_absolute_length <- average_absolute_length %>%
  dplyr::mutate(denomination_method = Denomination_ci(method_ci, optimize_in_a, bound_K, known_variance)) %>%
  dplyr::filter(is.na(bound_K) | (bound_K == 9)) %>%
  dplyr::filter(denomination_method %in% c("clt_with_unknown_variance",
                                           "navae_with_unknown_variance_optimized_in_a_with_K_9",
                                           "navae_with_known_variance_with_K_9"))

gr_absolute_length_for_latex <- ggplot(data = df_for_graph_absolute_length, 
                                       aes(x = n, y = `mean(length)`,
                                           color = denomination_method, linetype = denomination_method)) +
  geom_line(linewidth = linewidth_used) + 
  scale_color_manual(values = my_colors, guide = 'none') +
  scale_linetype_manual(values = my_linetypes, guide = 'none') +
  xlab(element_blank()) + ylab(element_blank()) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k"))

figure_filename <- paste0("fig_absolute_length_for_the_paper",
                          ".pdf")
pdf(file = file.path(path_figures, figure_filename), 
    width = width_export_pdf_file, height = height_export_pdf_file)
print(gr_absolute_length_for_latex)
dev.off()

#==============================================================================#
# Graph relative width wrt to CLT from simulations -----------------------------
#==============================================================================#

average_relative_length <- res_across_n %>%
  dplyr::group_by(n, method_ci, power_of_n_for_b, optimize_in_a, bound_K, known_variance) %>%
  dplyr::summarise(mean(ratio_length_wrt_CLT))

# Choice to use plug-in or K = 9. In the paper, we show K = 9
df_for_graph_relative_length <- average_relative_length %>%
  dplyr::mutate(denomination_method = Denomination_ci(method_ci, optimize_in_a, bound_K, known_variance)) %>%
  dplyr::filter(is.na(bound_K) | (bound_K == 9)) %>%
  dplyr::filter(denomination_method %in% c("navae_with_unknown_variance_optimized_in_a_with_K_9",
                                           "navae_with_known_variance_with_K_9"))

gr_relative_length_for_latex <- ggplot(data = df_for_graph_relative_length, 
                                       aes(x = n, y = `mean(ratio_length_wrt_CLT)`,
                                           color = denomination_method, linetype = denomination_method)) +
  geom_line(linewidth = linewidth_used) + 
  scale_color_manual(values = my_colors, guide = 'none') +
  scale_linetype_manual(values = my_linetypes, guide = 'none') +
  xlab(element_blank()) + ylab(element_blank()) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k"))

figure_filename <- paste0("fig_relative_length_for_the_paper",
                          ".pdf")
pdf(file = file.path(path_figures, figure_filename), 
    width = width_export_pdf_file, height = height_export_pdf_file)
print(gr_relative_length_for_latex)
dev.off()

#==============================================================================#
# Graph relative width wrt to CLT analytically using K = 9 ---------------------
#==============================================================================#

# Function to compute the relative length for known variance and 
# unknown variance (with optimized a_n in the case of unknown variance).
Compute_relative_length_of_NAVAE_ci_wrt_CLT_ci <- function(
    n, alpha = 0.10, K = 9,
    param_BE_EE = list(
      choice = "best",
      setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
      regularity = list(C0 = 1, p = 2),
      eps = 0.1))
{
  # Draw a sample of size n and apply our method to recover 
  # directly the ratio with respect to CLT's CI with unknown variance 
  # using an optimized a - it is not data-dependent, only throung n.
  # and recover delta_n to be used to compute the ratio in the known
  # variance case.
  res_ci <- Navae_ci_mean(data = rexp(n, rate = 1), alpha = alpha, 
                          optimize_in_a = TRUE, known_variance = NULL,
                          bound_K = K, param_BE_EE = param_BE_EE, verbose = TRUE)
  
  relative_length_with_unknown_variance <- res_ci$ratio_length_wrt_CLT
  
  relative_length_with_known_variance <- 
    qnorm(1 - alpha/2 + res_ci$delta_n) / qnorm(1 - alpha/2)
  
  return(c(with_unknown_variance = relative_length_with_unknown_variance,
           with_known_variance = relative_length_with_known_variance))
}

# Get the dataframe to plot the figures of relative length.
min_n <- 5 * 10^3
max_n <- 200 * 10^3
step_n <- 10

grid_n <- seq(from = min_n, to = max_n, by = step_n)
matrix_for_df <- matrix(NA, nrow = length(grid_n), ncol = 3)
colnames(matrix_for_df) <- c("n", 
                             "navae_with_unknown_variance_optimized_in_a_with_K_9",
                             "navae_with_known_variance_with_K_9")
# Takes 1, 2 minutes
for (i in seq_along(grid_n)) {
  n <- grid_n[[i]]
  if ((n %% 10000) == 0) {print(paste("Computing for n =", n))}
  relative_lengths <- Compute_relative_length_of_NAVAE_ci_wrt_CLT_ci(n)
  matrix_for_df[i, 1] <- n
  matrix_for_df[i, 2] <- relative_lengths["with_unknown_variance"]
  matrix_for_df[i, 3] <- relative_lengths["with_known_variance"]
}

df_for_graph_analytical_relative_length <- as.data.frame(matrix_for_df) %>% 
  pivot_longer(cols = -n, names_to = "ci", values_to = "relative_length")

# Suppress "linetype = ci": with the precise grid, dashed renders as solid,
# so distinction with colors only, ok.
gr_analytical_relative_length_for_latex <- 
  ggplot(data = df_for_graph_analytical_relative_length, 
         aes(x = n, y = relative_length,
             color = ci 
             # linetype = ci
         )) +
  geom_line(linewidth = linewidth_used) + 
  scale_color_manual(values = my_colors, guide = 'none') +
  # scale_linetype_manual(values = my_linetypes, guide = 'none') +
  xlab(element_blank()) + ylab(element_blank()) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k"))

figure_filename <- paste0("fig_analytical_relative_length_for_the_paper",
                          ".pdf")
pdf(file = file.path(path_figures, figure_filename), 
    width = width_export_pdf_file, height = height_export_pdf_file)
print(gr_analytical_relative_length_for_latex)
dev.off()

#==============================================================================#
# Study of inflexion point for delta_n with known or unknown variance ----------
#==============================================================================#

# Function to determine which bound BE or EE is used
# it is not data-dependent but we use again Navae_ci_mean function for simplicity.
# Remark that delta_n does not depend on a_n neither, just on n and K,
# so it is the same result with known or unknown variance.
Bound_delta_n_used <- function(
    n, alpha = 0.10, K = 9,
    param_BE_EE = list(
      choice = "best",
      setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
      regularity = list(C0 = 1, p = 2),
      eps = 0.1))
{
  res_ci <- Navae_ci_mean(data = rexp(n, rate = 1), alpha = alpha, 
                          known_variance = NULL, optimize_in_a = TRUE,
                          bound_K = K, param_BE_EE = param_BE_EE, verbose = TRUE)
  
  return(list(I_delta_n_from_BE = as.numeric(res_ci$delta_n_from == "BE"),
              delta_n = res_ci$delta_n))
}

grid_n <- seq(from = 35 * 10^3, to = 42 * 10^3, by = 10)
study_delta_n <- matrix(NA, nrow = length(grid_n), ncol = 3)
colnames(study_delta_n) <- c("n", "I_delta_n_from_BE", "delta_n")
for (i in seq_along(grid_n)) {
  n <- grid_n[[i]]
  if ((n %% 1000) == 0) {print(paste("Computing for n =", n))}
  res_bound_delta_n <- Bound_delta_n_used(n)
  study_delta_n[i, 1] <- n
  study_delta_n[i, 2] <- res_bound_delta_n$I_delta_n_from_BE
  study_delta_n[i, 3] <- res_bound_delta_n$delta_n
}

# Switch around:
study_delta_n[which.min(study_delta_n[,2]), 1]

# More precise grid to get exact n
grid_n <- seq(from = 38700, to = 38720, by = 1)
study_delta_n <- matrix(NA, nrow = length(grid_n), ncol = 3)
colnames(study_delta_n) <- c("n", "I_delta_n_from_BE", "delta_n")
for (i in seq_along(grid_n)) {
  n <- grid_n[[i]]
  res_bound_delta_n <- Bound_delta_n_used_with_known_variance(n)
  study_delta_n[i, 1] <- n
  study_delta_n[i, 2] <- res_bound_delta_n$I_delta_n_from_BE
  study_delta_n[i, 3] <- res_bound_delta_n$delta_n
}
study_delta_n
# Switch at:
switch_BE_to_EE <- study_delta_n[which.min(study_delta_n[,2]), 1]
switch_BE_to_EE

# Check on graph
gr_analytical_relative_length_for_latex + 
  geom_vline(xintercept = switch_BE_to_EE, colour = "black", linetype = "dashed")
