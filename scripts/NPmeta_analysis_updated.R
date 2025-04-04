# Script that explores the effect of N, P, and N+P addition on leaf- 
# and whole-plant functional traits. Meta-analysis uses data compiled
# from full-factorial N*P manipulation experiments. Whole-plant
# measurements are aggregated at the experiment level and leaf-level
# measurements are aggregated at the species level.

# Libraries
library(tidyverse)
library(metafor)
library(MAd)
library(ggpubr)
library(forcats)
library(patchwork)

# Read data sources (MESI, NutNet, EAP manual compilation)
mesi <- read.csv("../data/mesi_main_manual.csv")
nutnet <- read.csv("../data/nutnet_main.csv")
eap <- read.csv("../data/eap_main.csv")

# Merge data sources into single data frame
full_df <- mesi %>% full_join(nutnet) %>% full_join(eap)

# Add helper function for aggregating lnRRs across experiments
# (determines meta-analytic mean)
source("../helper_fxns/analyse_meta.R")


#####################################################################
# Investigate N fertilization responses
#####################################################################

# Subset `full_df` to include only N addition treatments
nfert_only <- full_df %>% filter(npk == "_100")
head(nfert_only)

# What traits are included?
unique(nfert_only$response)

# Select variables
use_response_n <- c("agb_n", "agb_p", "agb", "anpp", "agb_n_mass",
                    "agb_p_mass", "leaf_n_mass", "leaf_p_mass", "anet",
                    "total_biomass", "bgb", "rootshoot", "bnpp", "lma", "sla", "rmf",
                    "gsw", "amax", "leaf_wue", "fine_root_biomass", "leaf_n_area",
                    "leaf_p_area", "asat", "vcmax", "jmax", "leaf_pnue", "leaf_np",
                    "leaf_ppue", "leaf_structural_p", "leaf_metabolic_p",
                    "leaf_nucleic_p", "leaf_residual_p", "leaf_pi")

nfert_responses <- nfert_only %>%
  filter(response %in% use_response_n) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("amax", "asat"),
                        "asat", myvar))

use_vars_n <- unique(nfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
nfert_lnRR <- nfert_responses %>% 
  dplyr::select(-doi) %>%
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t),
         logr = ifelse(myvar == "sla", -logr, logr),
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_n <- purrr::map(as.list(use_vars_n),
                    ~analyse_meta(nfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_n) <- use_vars_n

df_box_n <- purrr::map_dfr(out_n, "df_box") |> 
  left_join(
    nfert_lnRR |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

#####################################################################
# Investigate P fertilization responses
#####################################################################

# Subset `full_df` to include only P addition treatments
pfert_only <- full_df %>% filter(npk == "_010")
head(pfert_only)

# What traits are included?
unique(pfert_only$response)

# Select variables
use_response_p <- c("agb_n", "agb_p", "agb", "anpp", "agb_n_mass",
                    "agb_p_mass", "leaf_n_mass", "leaf_p_mass", "anet",
                    "total_biomass", "bgb", "rootshoot", "bnpp", "lma", "rmf",
                    "gsw", "amax", "leaf_wue", "fine_root_biomass", "leaf_n_area",
                    "leaf_p_area", "asat", "vcmax", "jmax", "leaf_pnue", "leaf_np",
                    "leaf_ppue", "leaf_structural_p", "leaf_metabolic_p",
                    "leaf_nucleic_p", "leaf_residual_p", "leaf_pi")

pfert_responses <- pfert_only %>%
  filter(response %in% use_response_p) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("amax", "asat"),
                        "asat", myvar))

use_vars_p <- unique(pfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
pfert_lnRR <- pfert_responses %>% 
  dplyr::select(-doi) %>%
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t),
         logr = ifelse(myvar == "sla", -logr, logr),
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_p <- purrr::map(as.list(use_vars_p),
                    ~analyse_meta(pfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_p) <- use_vars_p

df_box_p <- purrr::map_dfr(out_p, "df_box") |> 
  left_join(
    pfert_lnRR |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

#####################################################################
# Investigate N+P fertilization responses
#####################################################################

# Subset `full_df` to include only N addition treatments
npfert_only <- full_df %>% filter(npk == "_110")
head(npfert_only)

# What traits are included?
unique(npfert_only$response)

# Select variables
use_response_np <- c("agb_n", "agb_p", "agb", "anpp", "agb_n_mass",
                    "agb_p_mass", "leaf_n_mass", "leaf_p_mass", "anet",
                    "total_biomass", "bgb", "rootshoot", "bnpp", "lma", "sla", "rmf",
                    "gsw", "amax", "leaf_wue", "fine_root_biomass", "leaf_n_area",
                    "leaf_p_area", "asat", "vcmax", "jmax", "leaf_pnue", "leaf_np",
                    "leaf_ppue", "leaf_structural_p", "leaf_metabolic_p",
                    "leaf_nucleic_p", "leaf_residual_p", "leaf_pi")

npfert_responses <- npfert_only %>%
  filter(response %in% use_response_np) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("amax", "asat"),
                        "asat", myvar))

use_vars_np <- unique(npfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
npfert_lnRR <- npfert_responses %>% 
  dplyr::select(-doi) %>%
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t),
         logr = ifelse(myvar == "sla", -logr, logr),
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_np <- purrr::map(as.list(use_vars_np),
                    ~analyse_meta(npfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_np) <- use_vars_np

df_box_np <- purrr::map_dfr(out_np, "df_box") |> 
  left_join(
    npfert_lnRR |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

##############################################################################
# Let's put together some plots
##############################################################################

# Add exp type to all data frames to merge together
df_box_np$manip_type <- "np"
df_box_p$manip_type <- "p"
df_box_n$manip_type <- "n"

npfert_lnRR$manip_type <- "np"
pfert_lnRR$manip_type <- "p"
nfert_lnRR$manip_type <- "n"

# Merge P and NP meta results
df_box_all <- df_box_n %>%
  full_join(df_box_p) %>%
  full_join(df_box_np) %>%
  mutate(manip_type = factor(manip_type, levels = c("np", "p", "n")),
         var = factor(var, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_n", "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma")))

fert_exp_responses_all <- nfert_lnRR %>%
  full_join(pfert_lnRR) %>%
  full_join(npfert_lnRR) %>%
  mutate(manip_type = factor(manip_type, levels = c("np", "p", "n")),
         myvar = factor(myvar, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_n", "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma")))
head(fert_exp_responses_all)

# Plot nutrients. Separating by trait type to avoid plot overwhelm
meta_plot_all_leaf_nutrients <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("lma", "sla", "leaf_n_mass", "leaf_n_area", "leaf_p_mass",
                             "leaf_p_area", "leaf_np")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se)) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("lma", "leaf_n_mass", "leaf_n_area", 
                                    "leaf_p_mass", "leaf_p_area", "leaf_np")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]), 
                              expression("M"["area"]))) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_leaf_nutrients

png("../plots/CNPmeta_leafNutrients_all")

# Plot leaf phosphorus fractionation. Separating by trait type to avoid plot overwhelm
meta_plot_all_phosFract <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("leaf_residual_p", 
                             "leaf_structural_p", "leaf_nucleic_p", 
                             "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se)) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("leaf_residual_p", 
                                    "leaf_structural_p", "leaf_nucleic_p", 
                                    "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf residual P",
                              "Leaf structural P",
                              "Leaf nucleic acid P",
                              "Leaf metabolic P",
                              "Leaf Pi")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_phosFract

# Plot photosynthetic traits. Separating by trait type to avoid plot overwhelm
meta_plot_all_photo <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("asat", "vcmax", "jmax",
                             "rd", "leaf_pnue", "leaf_ppue", "leaf_wue")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se)) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("asat", "vcmax", "jmax",
                                    "rd", "leaf_pnue", "leaf_ppue", "leaf_wue")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("iWUE",
                              "PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_photo

# Plot biomass traits. Separating by trait type to avoid plot overwhelm
meta_plot_all_biomass <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("rootshoot", "rmf", "bgb", "agb", "agb_n", 
                             "agb_p", "total_biomass")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se)) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("rootshoot", "rmf", "bgb", "agb", "agb_n", 
                                    "agb_p", "total_biomass")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "Belowground biomass",
                              "Aboveground P standing stock",
                              "Aboveground N standing stock",
                              "Aboveground biomass",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_biomass



png("../plots/CNPmeta_plot_all_combined_new.png", height = 12, width = 12, 
    units = "in", res = 600)
meta_plot_all_leaf_nutrients / meta_plot_all_photo / meta_plot_all_biomass +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 12, face = "bold"))
dev.off()

png("../plots/CNPmeta_plot_leafNutrients.png", height = 8, width = 12, 
    units = "in", res = 600)
meta_plot_all_leaf_nutrients
dev.off()

png("../plots/CNPmeta_plot_photo.png", height = 8, width = 12, 
    units = "in", res = 600)
meta_plot_all_photo
dev.off()

png("../plots/CNPmeta_phosFract.png", height = 4.5, width = 12, 
    units = "in", res = 600)
meta_plot_all_phosFract
dev.off()

png("../plots/CNPmeta_biomass.png", height = 8, width = 12, 
    units = "in", res = 600)
meta_plot_all_biomass
dev.off()

