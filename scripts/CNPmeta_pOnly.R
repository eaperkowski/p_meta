# CNPmeta_pOnly.R - script conducts meta-analysis for P addition experiments
# This is for the review paper Nick is leading on our theoretical and empirical
# understanding of phosphorus availability on ecophysiological traits and 
# includes experiments that include a phosphorus addition treatment
#
# Note: all paths assume that the root directory is the location of this script.

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

#####################################################################
# Set up and carry out P fertilization meta-analysis
#####################################################################

# What traits are included?
unique(df_pOnly$response)

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
                        "bgb", myvar),
         myvar = ifelse(myvar %in% "lai", "tla", myvar))

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

# Extract summary statistics for each model
out_p_modl <- lapply(out_p, `[[`, "modl")

out_p_modl_df <- data.frame(
  var = names(out_p_modl),
  nut_add = "p",
  k = sapply(out_p_modl, \(x) x$k),
  estimate = sapply(out_p_modl, \(x) as.numeric(coef(x))[1]),
  SE = sapply(out_p_modl, \(x) x$se[1]),
  zval = sapply(out_p_modl, \(x) x$zval[1]),
  pval = sapply(out_p_modl, \(x) x$pval[1]),
  ci.lb = sapply(out_p_modl, \(x) x$ci.lb[1]),
  ci.ub = sapply(out_p_modl, \(x) x$ci.ub[1]),
  row.names = NULL) |>
  mutate(across(estimate:zval, \(x) round(x, 3)),
         across(ci.lb:ci.ub, \(x) round(x, 3)),
         pval = ifelse(pval < 0.001, "<0.001",
                       round(pval, digits = 3)),
         k_plot = str_c("(", k, ")"))

##############################################################################
# Make some plots
##############################################################################

padd_chemistry_plot <- ggplot(data = out_p_modl_df %>% 
                                drop_na(var) %>% 
                                filter(var %in% c("lma", 
                                                  "leaf_n_mass", 
                                                  "leaf_n_area", 
                                                  "leaf_p_mass", 
                                                  "leaf_p_area",
                                                  "leaf_np")),
                              aes(x = var, y = estimate)) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]),
                              expression("P"["mass"]),
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  geom_text(aes(label = k_plot), y = 1.3, fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(-1, 1.5), breaks = seq(-1, 1.5, 0.5)) +
  labs(x = "", 
       y = NULL) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
padd_chemistry_plot

padd_photo_plot <- ggplot(data = out_p_modl_df %>% 
                            drop_na(var) %>% 
                            filter(var %in% c("asat",
                                              "rd",
                                              "vcmax", 
                                              "jmax",
                                              "leaf_pnue", 
                                              "leaf_ppue")),
                          aes(x = factor(var, levels = c("leaf_ppue",
                                                         "leaf_pnue",
                                                         "jmax",
                                                         "vcmax",
                                                         "rd",
                                                         "asat")),
                              
                              y = estimate)) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_text(aes(label = k_plot), y = 1.3, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("R"["d"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-1, 1.5), breaks = seq(-1, 1.5, 0.5)) +
  scale_fill_manual(values = "blue") +
  labs(x = "", y = NULL) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 18, color = "black"))
padd_photo_plot

padd_bio_plot <- ggplot(data = out_p_modl_df %>% 
                          drop_na(var) %>% 
                          filter(var %in% c("rootshoot", 
                                            "rmf", 
                                            "bnpp", 
                                            "anpp", 
                                            "tbio_gm2")),
                        aes(x = var, y = estimate)) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_text(aes(label = k_plot), y = 1.3, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "BNPP",
                              "ANPP",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-1, 1.5), breaks = seq(-1, 1.5, 0.5)) +
  labs(x = "", 
       y = "Log response to P addition")  +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
padd_bio_plot


padd_pfract_plot <- ggplot(data = out_p_modl_df %>% 
                             drop_na(var) %>% 
                             filter(var %in% c("leaf_residual_p", 
                                               "leaf_structural_p", 
                                               "leaf_nucleic_p", 
                                               "leaf_metabolic_p", 
                                               "leaf_pi")),
                           aes(x = factor(var,
                                          levels = c("leaf_residual_p", 
                                          "leaf_structural_p", 
                                          "leaf_nucleic_p", 
                                          "leaf_metabolic_p", 
                                          "leaf_pi")), 
                               y = estimate)) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_text(aes(label = k_plot), y = 2.5, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf residual P",
                              "Leaf structural P",
                              "Leaf nucleic P",
                              "Leaf metabolic P",
                              "Leaf inorganic P")) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  labs(x = "", 
       y = "Log response to P addition")  +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
padd_pfract_plot



# Figure 1 - individual effects
#png("../plots/CNP_pOnly_indEffects.png", height = 10, width = 6, 
#    units = "in", res = 600)
ggarrange(padd_chemistry_plot, padd_photo_plot, padd_bio_plot,
          nrow = 3, ncol = 1, labels = c("(a)", "(b)", "(c)"), 
          font.label = list(size = 18),
          align = "hv")
#dev.off()

# Figure 2 - P fractionation
#png("../plots/CNP_pOnly_Pfraction.png", height = 5, width = 8, 
#    units = "in", res = 600)
padd_pfract_plot
#dev.off()


