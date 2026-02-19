# CNPmeta_pOnly.R - script conducts meta-analysis for P addition experiments
# This is for the review paper on our theoretical and empirical understanding
# of phosphorus availability on ecophysiological traits and includes experiments
# that include a phosphorus addition treatment (as opposed to full-factorial
# experiments that add both phosphorus and nitrogen)

##############################################################################
# Libraries
##############################################################################
library(tidyverse)
library(metafor)
library(ggpubr)
library(forcats)
library(patchwork)
library(naniar) # to resolve NA/<NA> issue
library(orchaRd)

# Read compiled dataset
df_pOnly <- read.csv("../data/CNP_data_compiled.csv") %>%
  replace_with_na_all(~.x == "<NA>") %>%
  filter(fert == "p")

# Add helper function for aggregating lnRRs across experiments
# (determines meta-analytic mean)
source("../helper_fxns/analyse_meta.R")

##############################################################################
# Run meta-analysis
##############################################################################

# Select variables
use_response_p <- c("anpp_n", "anpp_p", "anpp", "leaf_n_mass", "leaf_p_mass", "anet",
                    "tbio_gm2", "bnpp", "rootshoot", "bgb", "lma",
                    "lai", "rmf", "tla", "gsw", "agb", "leaf_np", "total_biomass",
                    "asat", "leaf_wue", "leaf_n_area", "leaf_p_area", "vcmax", 
                    "jmax", "jmax_vcmax", "rd", "leaf_pnue", "leaf_ppue",
                    "leaf_structural_p", "leaf_metabolic_p", "leaf_nucleic_p",
                    "leaf_residual_p", "amax", "leaf_pi", "leaf_nre", "leaf_pre", 
                    "fine_root_biomass")

pfert_responses <- df_pOnly %>%
  filter(response %in% use_response_p) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass"),
                        "bgb", myvar))

use_vars_p <- unique(pfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
pfert_lnRR <- pfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_p <- purrr::map(as.list(use_vars_p),
                    ~analyse_meta(pfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_p) <- use_vars_p

df_box_p_k <- data.frame(var = use_vars_p,
                         k = c("(15)", "(12)", "(114)", "(160)", "(163)", "(115)",
                               "(33)", "(54)", "(61)", "(30)", "(123)",
                               "(5)", "(37)", "(17)", "(63)", "(41)", "(127)", 
                               "(33)", "(41)", "(87)", "(84)", "(52)", "(50)", 
                               "(37)", "(35)", "(65)", "(70)", "(24)", "(25)",
                               "(25)", "(24)", "(10)", "(6)", "(5)"),
                         sig.level = c("   ", "***", "***", "   ", "***", "***",
                                       ".  ", "   ", "***", "   ", "   ", 
                                       ".  ", ".  ", ".  ", "***", "***", "***",
                                       "***", "   ", "   ", "***", "***", "***", 
                                       "   ", "   ", ".  ", "*  ", "   ", "** ",
                                       "   ", "   ", "***", "   ", "   "))

df_box_p_k$k_sig <- str_c(df_box_p_k$k, df_box_p_k$sig.level)

df_box_p <- purrr::map_dfr(out_p, "df_box") |> 
  full_join(df_box_p_k) |>
  left_join(
    pfert_lnRR |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var") |>
  mutate(var = factor(var, 
                      levels = c("rootshoot", "rmf", "bgb", "bnpp", "agb_p", 
                                 "agb_p_mass", "agb_n", "agb_n_mass",
                                 "agb", "anpp", "total_biomass",
                                 "leaf_wue", "leaf_ppue", "leaf_pnue", "jmax_vcmax",
                                 "jmax", "vcmax", "gsw", "asat", "leaf_residual_p", 
                                 "leaf_structural_p", "leaf_nucleic_p", 
                                 "leaf_metabolic_p", "leaf_pi", 
                                 "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                 "leaf_n_area", "leaf_n_mass", "lma")),
         middle_fixed = ifelse(middle >= 0,
                               sprintf(" %.3f", middle),
                               sprintf("%.3f", middle)),
         plot_label = paste0("paste('", middle_fixed, "'^'", sig.level, "')"),
         bold_yn = ifelse(sig.level %in% c("*", "**", "***"),
                          "bold", "plain"),
         middle_perc = (exp(middle) - 1) * 100,
         upper_perc = (exp(ymax) - 1) * 100,
         lower_perc = (exp(ymin) - 1) * 100)

##############################################################################
# Make some plots
##############################################################################

padd_chemistry_plot <- ggplot(data = df_box_p %>% 
                                drop_na(var) %>% 
                                filter(var %in% c("lma", 
                                                  "leaf_n_mass", 
                                                  "leaf_n_area", 
                                                  "leaf_p_mass", 
                                                  "leaf_p_area",
                                                  "leaf_np")),
                              aes(x = var, y = middle)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]),
                              expression("P"["mass"]),
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  geom_text(aes(label = k_sig), y = 1, fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = "", 
       y = NULL) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
padd_chemistry_plot

padd_photo_plot <- ggplot(data = df_box_p %>% 
                            drop_na(var) %>% 
                            filter(var %in% c("asat", 
                                              "vcmax", 
                                              "jmax",
                                              "leaf_pnue", 
                                              "leaf_ppue")),
                          aes(x = var, y = middle)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_text(aes(label = k_sig), y = 1.5, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-1, 1.5), breaks = seq(-1, 1.5, 0.5)) +
  scale_fill_manual(values = "blue") +
  labs(x = "", y = "Log response to P addition") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 18, color = "black"))
padd_photo_plot

padd_bio_plot <- ggplot(data = df_box_p %>% 
                          drop_na(var) %>% 
                          filter(var %in% c("rootshoot", 
                                            "rmf", 
                                            "bnpp", 
                                            "anpp", 
                                            "total_biomass")),
                        aes(x = var, y = middle)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_text(aes(label = k_sig), y = 1.5, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "BNPP",
                              "ANPP",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-1, 1.5), breaks = seq(-1, 1.5, 0.5)) +
  labs(x = "", 
       y = NULL)  +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
padd_bio_plot


padd_pfract_plot <- ggplot(data = df_box_p %>% 
                             drop_na(var) %>% 
                             filter(var %in% c("leaf_residual_p", 
                                               "leaf_structural_p", 
                                               "leaf_nucleic_p", 
                                               "leaf_metabolic_p", 
                                               "leaf_pi")),
                           aes(x = var, y = middle)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_text(aes(label = k_sig), y = 3, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf residual P",
                              "Leaf structural P",
                              "Leaf nucleic P",
                              "Leaf metabolic P",
                              "Leaf inorganic P")) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  labs(x = "", 
       y = NULL)  +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
padd_pfract_plot





