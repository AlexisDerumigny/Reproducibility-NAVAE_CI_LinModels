#' Script for creating tables and figures used in the paper 
#' for the simulations in case of OLS
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

path_figures <- file.path(path_dropbox, "finite_CI_perso", 
                          "Reproducibility-NAVAE_CI_LinModels", "figures")

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
# Load results of simulations --------------------------------------------------
#==============================================================================#

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

# For nominal level
alpha <- 0.10

# For the coeff considered
coeff_considered <- 2

# constant_for_b_n <- 1
# power_of_n_for_b <- -1/5
# constant_for_b_n <- 100 # with 100 * n^{-2/5}, strange behavior in the simu for n = 100k, quite weird!
constant_for_b_n <- 20
power_of_n_for_b <- -2/5

power_of_n_for_omega <- -1/5

vec_n <- c(5, 10, 25, 50, 75, 100, 125, 150, 175, 200) * 10^3

M <- 200

# Load results
res_over_n <- NULL
for (n in vec_n) {
  res_for_one_n <- readRDS(
    file = file.path(path_res_simu_OLS, 
                     Get_filename_OLS(
                       coeff_considered, alpha, constant_for_b_n, power_of_n_for_b, power_of_n_for_omega, n, M)))
  if (is.null(res_over_n)) {
    res_over_n <- res_for_one_n
  } else {
    res_over_n <- rbind(res_over_n, res_for_one_n)
  }
  rm(res_for_one_n)
}

#==============================================================================#
# Raw results for coverage and widths (absolute and relative) ------------------
#==============================================================================#

# Coverage

coverage <- res_over_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(mean(coverage_navae_K9), mean(coverage_navae_Kh), mean(coverage_asymp))
coverage
# We will not report the table, since 100%, just discuss it.

# Proportion of R regime

prop_R_regime <- res_over_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(mean(ci_navae_K9.regime != "Edg"),
                   mean(ci_navae_Kh.regime != "Edg"))
prop_R_regime

# Length

median_absolute_length <- res_over_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(median(length_navae_K9), median(length_navae_Kh), median(length_asymp))
median_absolute_length

mean_absolute_length <- res_over_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(mean(length_navae_K9), mean(length_navae_Kh), mean(length_asymp))
mean_absolute_length

median_ratio_length <- res_over_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(median(ratio_wrt_asymp_navae_K9), median(ratio_wrt_asymp_navae_Kh))
median_ratio_length

mean_ratio_length <- res_over_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(mean(ratio_wrt_asymp_navae_K9), mean(ratio_wrt_asymp_navae_Kh))
mean_ratio_length

# Just to check the behavior of the median with Inf, treated as a value itself, so can be Inf
# if at least half the values are Inf.
# median(c(1,2,Inf,Inf))
# median(c(1,2, 3,Inf,Inf))

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
my_colors <- c("clt" = my_green,
               "navae_K9" = my_red)

my_linetypes <- c("clt" = "dotted",
                  "navae_K9" = "solid")

#==============================================================================#
# Figure absolute width, mean over MC repetitions, K = 9 only ------------------
#==============================================================================#

df_for_gr_mean_absolute_length <- res_over_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(mean(length_navae_K9), mean(length_asymp)) %>%
  dplyr::rename(navae_K9 = `mean(length_navae_K9)`,
                clt = `mean(length_asymp)`) %>%
  pivot_longer(cols = -n, names_to = "ci", values_to = "absolute_length")

gr_ols_mean_absolute_length <- ggplot(data = df_for_gr_mean_absolute_length, 
                                      aes(x = n, y = absolute_length,
                                          color = ci, linetype = ci)) +
  geom_line(linewidth = linewidth_used) + 
  scale_color_manual(values = my_colors, guide = 'none') +
  scale_linetype_manual(values = my_linetypes, guide = 'none') +
  xlab(element_blank()) + ylab(element_blank()) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k"))

figure_filename <- paste0("fig_ols_absolute_length_for_the_paper",
                          ".pdf")
pdf(file = file.path(path_figures, figure_filename), 
    width = width_export_pdf_file, height = height_export_pdf_file)
print(gr_ols_mean_absolute_length)
dev.off()

#==============================================================================#
# Figure relative width, mean over MC repetitions, K = 9 only ------------------
#==============================================================================#

df_for_gr_mean_relative_length <- res_over_n %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(mean(ratio_wrt_asymp_navae_K9)) %>%
  dplyr::rename(navae_K9 = `mean(ratio_wrt_asymp_navae_K9)`) %>%
  pivot_longer(cols = -n, names_to = "ci", values_to = "relative_length")

gr_ols_mean_relative_length <- ggplot(data = df_for_gr_mean_relative_length, 
                                      aes(x = n, y = relative_length,
                                          color = ci, linetype = ci)) +
  geom_line(linewidth = linewidth_used) + 
  scale_color_manual(values = my_colors, guide = 'none') +
  scale_linetype_manual(values = my_linetypes, guide = 'none') +
  xlab(element_blank()) + ylab(element_blank()) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  scale_y_continuous(limits = c(1, max(df_for_gr_mean_relative_length$relative_length)), 
                     breaks = c(1, 6, 11, 16))
  # Personalized y axis to have 1 (the limiting ratio)

figure_filename <- paste0("fig_ols_relative_length_for_the_paper",
                          ".pdf")
pdf(file = file.path(path_figures, figure_filename), 
    width = width_export_pdf_file, height = height_export_pdf_file)
print(gr_ols_mean_relative_length)
dev.off()
