# Script that explores the effect of N, P, and N+P addition on leaf- 
# and whole-plant functional traits. Meta-analysis uses data compiled
# from full-factorial N*P manipulation experiments. Whole-plant
# measurements are aggregated at the experiment level and leaf-level
# measurements are aggregated at the species level.

# Libraries
library(tidyverse)
library(metafor)
library(ggpubr)
library(forcats)
library(patchwork)
library(naniar) # to resolve NA/<NA> issue
library(orchaRd)

# Read compiled dataset
full_df <- read.csv("../data/CNP_data_compiled.csv") %>%
  replace_with_na_all(~.x == "<NA>")

# Create experiment metadata summary
experiment_summary <- full_df %>%
  dplyr::select(citation, exp:experiment_type) %>%
  distinct(citation, exp, .keep_all = TRUE)

# Species identity summary
species_summary <- full_df %>% 
  distinct(species, family, growth_form, growth_duration,
           photo_path, n_fixer, myc_assoc, .keep_all = TRUE) %>%
  dplyr::select(species:myc_assoc) %>%
  slice(-1)

# How many experiments?
unique(experiment_summary$citation)

# How many sites?
unique(experiment_summary$exp)

# How many field experiments?
unique(filter(experiment_summary, experiment_type == "field")$citation)
unique(filter(experiment_summary, experiment_type == "greenhouse")$citation)
unique(filter(experiment_summary, experiment_type == "chamber")$citation)

# How many species?
distinct(full_df, species)
distinct(full_df, family)

# Add helper function for aggregating lnRRs across experiments
# (determines meta-analytic mean)
source("../helper_fxns/analyse_meta.R")

# Add helper function for calculating interaction effect sizes
# (uses Hedge's d)
source("../helper_fxns/calc_intxn_effSize_meta.R")

#####################################################################
# Investigate N fertilization responses
#####################################################################

# Subset `full_df` to include only N addition treatments
nfert_only <- full_df %>% filter(npk == "_100")
head(nfert_only)

# What traits are included?
unique(nfert_only$response)

# Which experiments are included?
unique(nfert_only$citation)

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
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("total_leaf_area", "tla"),
                        "tla", myvar))

use_vars_n <- unique(nfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`. Also, convert SLA 
# to Marea (logr = opposite sign, stdev/error/variance remains the same)

nfert_lnRR <- nfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t))  %>%
  mutate(logr = ifelse(myvar == "sla", logr * -1, logr), 
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_n <- purrr::map(as.list(use_vars_n),
                    ~analyse_meta(nfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_n) <- use_vars_n

df_box_n_k <- data.frame(var = use_vars_n,
                          k = c("(18)", "(12)", "(125)", "(30)", "(23)", "(139)",
                                "(133)", "(85)", "(42)", "(63)", "(40)",
                                NA, "(37)", "(43)", "(115)", "(113)", "(37)", 
                                "(84)", "(79)", "(42)", "(40)", "(58)", "(59)", 
                                "(20)", "(21)", "(21)", "(20)", "(5)"),
                         sig.level = c("***", "***", "***", "***", "***", "***",
                                       "*. ", "ns ", "ns ", "ns ", "***", 
                                       NA, "***", "ns ", "** ", "*. ", "ns ", 
                                       "** ", "ns ", "ns ", "ns ", "ns ", "ns ", 
                                       "ns ", "ns ", "ns ", "ns ", "ns "))

df_box_n <- purrr::map_dfr(out_n, "df_box") |> 
  full_join(df_box_n_k) |>
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
                    "total_biomass", "bgb", "rootshoot", "bnpp", "lma", "sla", "rmf",
                    "gsw", "amax", "leaf_wue", "fine_root_biomass", "leaf_n_area",
                    "leaf_p_area", "asat", "vcmax", "jmax", "leaf_pnue", "leaf_np",
                    "leaf_ppue", "leaf_structural_p", "leaf_metabolic_p",
                    "leaf_nucleic_p", "leaf_residual_p", "leaf_pi")

pfert_responses <- pfert_only %>%
  filter(response %in% use_response_p) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("total_leaf_area", "tla"),
                        "tla", myvar))

use_vars_p <- unique(pfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
pfert_lnRR <- pfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t)) %>%
  mutate(logr = ifelse(myvar == "sla", logr * -1, logr), 
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_p <- purrr::map(as.list(use_vars_p),
                    ~analyse_meta(pfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_p) <- use_vars_p

df_box_p_k <- data.frame(var = use_vars_p,
                         k = c("(18)", "(12)", "(125)", "(30)", "(23)", "(139)",
                               "(133)", "(85)", "(42)", "(63)", "(40)", NA,
                               "(37)", "(43)", "(115)", "(113)", "(37)", "(84)", 
                               "(79)", "(42)", "(40)", "(58)", "(59)", "(20)",
                               "(21)", "(21)", "(20)", "(5)"),
                         sig.level = c(".  ", "***", "***", "***", "***", "ns ",
                                       "***", "ns ", "* ", "ns ", "*. ", NA,
                                       ".  ", ".  ", "***", "ns ", "ns ", "ns ", 
                                       "***", "ns ", "*. ", "ns ", "ns ", "ns ",
                                       "*. ", "ns ", "ns ", "**"))


df_box_p <- purrr::map_dfr(out_p, "df_box") |> 
  full_join(df_box_p_k) |>
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
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("total_leaf_area", "tla"),
                        "tla", myvar)) %>%
  filter(myvar != "r_eco" & myvar != "nee" & 
           myvar != "spad" & myvar != "tpu" &
           myvar != "rd" & myvar != "stom_lim" &
           myvar != "anet_mass")

use_vars_np <- unique(npfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
npfert_lnRR <- npfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t)) %>%
  mutate(logr = ifelse(myvar == "sla", logr * -1, logr), 
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_np <- purrr::map(as.list(use_vars_np),
                    ~analyse_meta(npfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_np) <- use_vars_np

df_box_np_k <- data.frame(var = use_vars_np,
                          k = c("(18)", "(12)", "(125)", "(30)", "(23)", "(139)",
                                "(133)", "(85)", "(42)", "(63)", "(40)", NA,
                                "(37)", "(43)", "(115)", "(113)", "(37)", "(84)", 
                                "(79)", "(42)", "(40)", "(58)", "(59)", "(20)",
                                "(21)", "(21)", "(20)", "(5)"),
                          sig.level = c("***", "***", "***", "** ", "*  ", "***",
                                        "***", "*  ", "***", "ns ", "***", NA,
                                        "** ", "*  ", "***", "** ", ".  ", "***", 
                                        "***", ".  ", "***", "ns ", "ns ", "ns ",
                                        "*  ", "ns ", "ns ", "** "))

df_box_np <- purrr::map_dfr(out_np, "df_box") |> 
  full_join(df_box_np_k) %>%
  left_join(
    npfert_lnRR |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

df_box_np$sig.level
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

# Merge N, P, and NP meta results
df_box_all <- df_box_n %>%
  full_join(df_box_p) %>%
  full_join(df_box_np) %>%
  filter(!is.na(middle)) %>%
  mutate(manip_type = factor(manip_type, levels = c("np", "p", "n")),
         var = factor(var, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_p_mass", "agb_n", "agb_n_mass",
                                   "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "gsw", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma"))) %>%
  mutate(middle_fixed = ifelse(middle >= 0,
                               sprintf(" %.3f", middle),
                               sprintf("%.3f", middle)),
         plot_label = paste0("paste('", middle_fixed, "'^'", sig.level, "')"),
         bold_yn = ifelse(sig.level %in% c("*", "**", "***"),
                          "bold", "plain"))

# Check that all vars are incorporated
unique(df_box_all$var)
## Yes

df_box_all$plot_label

# Merge N, P, and NP lnRR values
fert_exp_responses_all <- nfert_lnRR %>%
  full_join(pfert_lnRR) %>%
  full_join(npfert_lnRR) %>%
  mutate(manip_type = factor(manip_type, levels = c("np", "p", "n")),
         myvar = factor(myvar, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_p_mass", "agb_n", "agb_n_mass",
                                   "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "gsw", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma")))
unique(fert_exp_responses_all$myvar)

# Plot nonphotosynthetic leaf traits. Separating by trait type to avoid 
# plot overwhelm
meta_plot_all_leaf_nutrients <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("sla", "lma", "leaf_n_mass", "leaf_n_area", "leaf_p_mass",
                             "leaf_p_area", "leaf_np")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_var)) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("sla", "lma", "leaf_n_mass", "leaf_n_area", 
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
                              expression("M"["area"]),
                              "SLA")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 5)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = "Weight") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_leaf_nutrients

# Plot nonphotosynthetic leaf traits WITHOUT individual data points. Separating 
# by trait type to avoid plot overwhelm
meta_leaftraits_plot_noPoints <- ggplot(
  data = df_box_all %>% filter(var %in% c("lma", "leaf_n_mass", 
                                          "leaf_n_area", "leaf_p_mass", 
                                          "leaf_p_area", "leaf_np"))) +
  geom_crossbar(data = df_box_all %>%
                  filter(var %in% c("lma", "leaf_n_mass", "leaf_n_area",
                                    "leaf_p_mass", "leaf_p_area", "leaf_np")),
                aes(x = var, y = middle, ymin = ymin,
                    ymax = ymax, fill = manip_type),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.75)) +
  geom_text(aes(x = var, y = 0.8, 
                label = plot_label, 
                group = manip_type,
                fontface = bold_yn),
            position = position_dodge(width = 0.75),
            parse = TRUE, hjust = 0) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]), 
                              expression("M"["area"]))) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 5)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition") +
  coord_flip() +
  theme_bw(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
meta_leaftraits_plot_noPoints






# Plot leaf phosphorus fractionation. Separating by trait type to avoid plot overwhelm
meta_plot_all_phosFract <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("leaf_residual_p", 
                             "leaf_structural_p", "leaf_nucleic_p", 
                             "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_var)) +
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
  scale_size(range = c(0.25, 5)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = "Weight") +
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
              shape = 21, aes(size = 1/logr_var)) +
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
  scale_size(range = c(0.25, 5)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = "Weight") +
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
              shape = 21, aes(size = 1/logr_var)) +
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
       size = "Weight") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_biomass


# Write plots
# png("../plots/CNPmeta_plot_all_combined_new.png", height = 12, width = 12, 
#     units = "in", res = 600)
meta_plot_all_leaf_nutrients / meta_plot_all_photo / meta_plot_all_biomass +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 12, face = "bold"))
# dev.off()

# png("../plots/CNPmeta_plot_leafNutrients.png", height = 8, width = 12, 
#    units = "in", res = 600)
# meta_plot_all_leaf_nutrients
# dev.off()
# 
# png("../plots/CNPmeta_plot_photo.png", height = 8, width = 12, 
#     units = "in", res = 600)
# meta_plot_all_photo
# dev.off()
# 
# png("../plots/CNPmeta_phosFract.png", height = 4.5, width = 12, 
#     units = "in", res = 600)
# meta_plot_all_phosFract
# dev.off()
# 
# png("../plots/CNPmeta_biomass.png", height = 8, width = 12, 
#     units = "in", res = 600)
# meta_plot_all_biomass
# dev.off()

##############################################################################
# Some prep work for calculating interaction effect sizes (need
# to turn summary statistics into long format)
##############################################################################

# Create df that only includes control trt, prep for merge to long format
full_df_control <- full_df %>%
  mutate(unique_id = ifelse(
    is.na(species),
    str_c(citation, "XXX", exp, "XXX", response, "XXX", "none"),
    str_c(citation, "XXX", exp, "XXX", response, "XXX", species))) %>%
  group_by(unique_id) %>%
  summarize(x_c = mean(x_c),
            sd_c = mean(sd_c),
            rep_c = mean(rep_c))

# Convert df into long format to make for easy effect size calculation
full_df_trt_long <- full_df %>% 
  mutate(unique_id = ifelse(is.na(species), 
                            str_c(citation, "XXX", exp, "XXX", response, "XXX", "none"),
                            str_c(citation, "XXX", exp, "XXX", response, "XXX", species))) %>%
  dplyr::select(unique_id, fert, x_t, sd_t, rep_t) %>%
  pivot_wider(names_from = fert, values_from = c(x_t, sd_t, rep_t)) %>%
  full_join(full_df_control) %>%
  dplyr::select(unique_id,
                x_c, sd_c, rep_c,
                x_n = x_t_n, sd_n = sd_t_n, rep_n = rep_t_n,
                x_p = x_t_p, sd_p = sd_t_p, rep_p = rep_t_p,
                x_np = x_t_np, sd_np = sd_t_np, rep_np = rep_t_np)

###############################################################################
# Calculate individual, main, and interaction effect sizes
###############################################################################
CNP_effect_sizes <- data.frame(
  unique_id = full_df_trt_long$unique_id,
  calc_intxn_effSize_meta(x_a = full_df_trt_long$x_n, 
                          s_a = full_df_trt_long$sd_n, 
                          n_a = full_df_trt_long$rep_n,
                          
                          x_b = full_df_trt_long$x_p, 
                          s_b = full_df_trt_long$sd_p, 
                          n_b = full_df_trt_long$rep_p,
                          
                          x_c = full_df_trt_long$x_c, 
                          s_c = full_df_trt_long$sd_c, 
                          n_c = full_df_trt_long$rep_c,
                          
                          x_ab = full_df_trt_long$x_np, 
                          s_ab = full_df_trt_long$sd_np, 
                          n_ab = full_df_trt_long$rep_np)) %>%
  separate(unique_id, 
           into = c("citation", "exp", "response", "species"), 
           sep = "XXX", remove = TRUE) %>%
  full_join(experiment_summary, by = c("citation", "exp")) %>%
  replace_with_na_all(~.x == "none") %>%
  full_join(species_summary, by = "species") %>%
  dplyr::select(citation, exp, latitude:experiment_type, 
                species, family:myc_assoc,  response,
                
                # Individual effects
                gNi = g_a, vNi = v_a, wNi = w_a,
                gPi = g_b, vPi = v_b, wPi = w_b,
                gNPi = g_ab, vNPi = v_ab, wNPi = w_ab,
                
                # Main effects
                dNi = dA, dvNi = v_a_main, dwNi = w_a_main,
                dPi = dB, dvPi = v_b_main, dwPi = w_b_main,
                
                # Interaction effects
                dNPi = dAB, dvNPi = v_ab_int, dwNPi = w_ab_int)

###############################################################################
# Clean CNP_effect_sizes to only include relevant traits. Also tidy up to
# group similar traits without losing integrity of original trait value in
# `response` column
###############################################################################
CNP_effect_sizes_reduced <- CNP_effect_sizes %>%
  mutate(myvar = response,
         myvar = ifelse(myvar == "amax", "asat", myvar),
         myvar = ifelse(myvar == "anet", "asat", myvar),
         myvar = ifelse(myvar == "fine_root_biomass" | myvar == "bnpp", 
                      "bgb", myvar),
         myvar = ifelse(myvar == "anpp", "agb", myvar),
         myvar = ifelse(myvar == "total_leaf_area", "tla", myvar)) %>%
  filter(myvar != "r_eco" & myvar != "nee" & 
           myvar != "spad" & myvar != "tpu" &
           myvar != "rd" & myvar != "stom_lim" &
           myvar != "anet_mass" & myvar!= "lai" & myvar != "agb_np" &
           myvar != "leaf_sugar_p" & myvar != "gpp")

###############################################################################
# Use helper fxn from `analyse_meta.R` to iterate through traits and determine
# meta-analytic mean of the interaction effect size and associated 95% CI
###############################################################################

# Determine interaction effect size across multiple experiments. 
# Takes into account experiment variance and uses experiment 
# identity as grouping factor for random intercepts
use_vars_int <- unique(CNP_effect_sizes_reduced$myvar)

out_int <- purrr::map(as.list(use_vars_int),
                     ~analyse_meta_int(CNP_effect_sizes_reduced %>%
                                     rename(var = myvar), nam_target = .))
names(out_int) <- use_vars_int

df_box_int <- purrr::map_dfr(out_int, "df_box") |> 
  left_join(
    CNP_effect_sizes_reduced |> 
      group_by(myvar) |> 
      summarise(intES_min = min(dNPi), intES_max = max(dNPi)) |> 
      rename(var = myvar),
    by = "var")

# Determine main effect size across multiple experiments. 
# Takes into account experiment variance and uses experiment 
# identity as grouping factor for random intercepts

out_mainN <- purrr::map(as.list(use_vars_int),
                        ~analyse_meta_int(CNP_effect_sizes_reduced %>%
                                            rename(var = myvar), nam_target = .))
names(out_int) <- use_vars_int


###############################################################################
# Let' make some plots!
###############################################################################
# Factor variables to appear in plots in a certain order
CNP_effect_sizes_reduced <- CNP_effect_sizes_reduced %>%
  mutate(myvar = factor(myvar, 
               levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                          "agb_n", "agb", "total_biomass", 
                          "leaf_wue", "leaf_ppue", "leaf_pnue", 
                          "jmax", "vcmax", "asat", "leaf_residual_p", 
                          "leaf_structural_p", "leaf_nucleic_p", 
                          "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi", 
                          "leaf_np", "leaf_p_area", "leaf_p_mass", 
                          "leaf_n_area", "leaf_n_mass", "lma", "sla")))

# Factor variables to appear in plots in a certain order
df_box_int <- df_box_int %>%
  mutate(var = factor(var, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_n", "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma", "sla")),
         int_type = ifelse(var == "leaf_np", 
                           "synergistic",
                           ifelse(var == "agb", 
                                  "synergistic",
                                  "additive")))

meta_intPlot_leaf_nutrients <- ggplot(
  data = CNP_effect_sizes_reduced %>%
    filter(myvar %in% c("sla", "lma", "leaf_n_mass", "leaf_n_area", 
                        "leaf_p_mass", "leaf_p_area", "leaf_np")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = dwNPi)) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("sla", "lma", "leaf_n_mass", "leaf_n_area", 
                                    "leaf_p_mass", "leaf_p_area", "leaf_np")),
                aes(x = var, y = middle, ymin = ymin, 
                    ymax = ymax, fill = int_type),
                alpha = 0.6, width = 0.4,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]), 
                              expression("M"["area"]),
                              "SLA")) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1)) +
  labs(x = "", 
       y = "Interaction effect size (Hedge's d)",
       size = "Weight",
       fill = "Interaction type") +
  scale_size(range = c(0.25, 4)) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_intPlot_leaf_nutrients

# Plot leaf phosphorus fractionation. Separating by trait type to avoid plot overwhelm
meta_intPlot_phosFract <- ggplot(
  data = CNP_effect_sizes_reduced %>%
    filter(myvar %in%  c("leaf_residual_p", 
                         "leaf_structural_p", "leaf_nucleic_p", 
                         "leaf_metabolic_p", "leaf_pi")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = dwNPi)) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("leaf_residual_p", 
                                    "leaf_structural_p", "leaf_nucleic_p", 
                                    "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi")),
                aes(x = var, y = middle, ymin = ymin, 
                    ymax = ymax, fill = int_type),
                alpha = 0.6, width = 0.4,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf residual P",
                              "Leaf structural P",
                              "Leaf nucleic acid P",
                              "Leaf metabolic P",
                              "Leaf Pi")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Interaction effect size (Hedge's d)",
       size = "Weight",
       fill = "Interaction type") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_intPlot_phosFract

# Plot photosynthetic traits. Separating by trait type to avoid plot overwhelm
meta_intPlot_photo <- ggplot(
  data = CNP_effect_sizes_reduced %>%
    filter(myvar %in% c("asat", "vcmax", "jmax",
                        "rd", "leaf_pnue", "leaf_ppue", "leaf_wue")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = dwNPi)) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("asat", "vcmax", "jmax",
                                    "rd", "leaf_pnue", "leaf_ppue", "leaf_wue")),
                aes(x = var, y = middle, ymin = ymin, 
                    ymax = ymax, fill = int_type),
                alpha = 0.6, width = 0.4, 
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("iWUE",
                              "PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_size(range = c(0.25, 4)) +
  labs(x = NULL, 
       y = "Interaction effect size (Hedge's d)",
       fill = "Interaction type",
       size = "Weight") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_intPlot_photo

# Plot biomass traits. Separating by trait type to avoid plot overwhelm
meta_intPlot_biomass <- ggplot(
  data = subset(CNP_effect_sizes_reduced, 
                myvar %in% c("rootshoot", "rmf", "bgb", "agb", "agb_n", 
                             "agb_p", "total_biomass")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = dwNPi)) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("rootshoot", "rmf", "bgb", "agb", "agb_n", 
                                    "agb_p", "total_biomass")),
                aes(x = var, y = middle, ymin = ymin, 
                    ymax = ymax, fill = int_type),
                alpha = 0.6, width = 0.4,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "Belowground biomass",
                              "Aboveground P standing stock",
                              "Aboveground N standing stock",
                              "Aboveground biomass",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Interaction effect size (Hedge's d)",
       fill = "Interaction type",
       size = "Weight") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_intPlot_biomass

# png("../plots/CNPmeta_interaction_plot.png", height = 12, width = 12, 
#    units = "in", res = 600)
meta_intPlot_leaf_nutrients / meta_intPlot_photo / meta_intPlot_biomass +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 12, face = "bold"))
# dev.off()

# png("../plots/CNPmeta_intplot_leafNutrients.png", height = 8, width = 12, 
#     units = "in", res = 600)
meta_intPlot_leaf_nutrients
# dev.off()

# png("../plots/CNPmeta_intplot_photo.png", height = 8, width = 12, 
#     units = "in", res = 600)
meta_intPlot_photo
# dev.off()

# png("../plots/CNPmeta_intplot_phosFract.png", height = 4.5, width = 12, 
#     units = "in", res = 600)
meta_intPlot_phosFract
# dev.off()

# png("../plots/CNPmeta_intplot_biomass.png", height = 8, width = 12, 
#     units = "in", res = 600)
meta_intPlot_biomass
# dev.off()

###############################################################################
# lnRR trait summaries
###############################################################################
tableXX_lnRR_summary <- data.frame(
  trait = df_box_all$var,
  response = df_box_all$manip_type,
  lnRR = round(df_box_all$middle, digits = 3),
  lowerCI = round(df_box_all$ymin, digits = 3),
  upperCI = round(df_box_all$ymax, digits = 3),
  ci_range = str_c("(", 
                   round(df_box_all$ymin, digits = 3), ", ", 
                   round(df_box_all$ymax, digits = 3), ")")) %>%
  filter(trait %in% c("lma", "leaf_n_mass", "leaf_n_area",
                      "leaf_p_mass", "leaf_p_area", "leaf_np",
                      "leaf_pi", "leaf_sugar_p", "leaf_metabolic_p",
                      "leaf_nucleic_p", "leaf_structural_p", "asat", 
                      "vcmax", "jmax", "leaf_pnue", "leaf_ppue", 
                      "leaf_wue", "total_biomass", "agb", "agb_n", 
                      "agb_p", "bgb", "rmf", "rootshoot")) %>%
  mutate(trait = factor(trait, 
                        levels = c("lma", "leaf_n_mass", "leaf_n_area",
                                   "leaf_p_mass", "leaf_p_area", "leaf_np",
                                   "leaf_pi", "leaf_sugar_p", "leaf_metabolic_p",
                                   "leaf_nucleic_p", "leaf_structural_p", "asat", 
                                   "vcmax", "jmax", "leaf_pnue", "leaf_ppue", 
                                   "leaf_wue", "total_biomass", "agb", "agb_n", 
                                   "agb_p", "bgb", "rmf", "rootshoot")),
         sig = ifelse((lowerCI > 0 & upperCI > 0) |
                        (lowerCI < 0 & upperCI < 0), "***",
                      "NS")) %>%
  arrange(trait)

filter(table1_lnRR_summary, response == "np")
  
###############################################################################
# Re-run Nmass models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Nmass - full model
nfert_nmass_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + myc_assoc + photo_path + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_n_mass"))
summary(nfert_nmass_fullModel)

# Nfert: Nmass - MAT plot
nfert_nmass_mat_plot <- mod_results(nfert_nmass_fullModel, 
                              mod = "mat",
                              group = "exp", subset = TRUE
                              )$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-5, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of N"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# Nfert: Nmass - AI plot
nfert_nmass_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of N"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Nmass - PAR plot
nfert_nmass_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of N"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Nmass - full model
pfert_nmass_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_n_mass"))
summary(pfert_nmass_fullModel)

# Pfert: Nmass - MAT plot
pfert_nmass_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-5, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of N"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# Pfert: Nmass - AI plot
pfert_nmass_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of N"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Nmass - PAR plot
pfert_nmass_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of N"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))


# N+P addition effect on Nmass - full model
npfert_nmass_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "leaf_n_mass"))
summary(npfert_nmass_fullModel)

# N+P fert: Nmass - MAT plot
npfert_nmass_mat_plot <- mod_results(npfert_nmass_fullModel, 
                                     mod = "mat", 
                                     group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-5, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of N"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# N+P fert: Nmass - AI plot
npfert_nmass_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of N"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P fert: Nmass - PAR plot
npfert_nmass_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_mass"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of N"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

png("../plots/CNPmeta_Nmass_clim_moderators.png", units = "in",
    height = 14, width = 16, res = 600)
ggarrange(nfert_nmass_mat_plot, nfert_nmass_ai_plot, nfert_nmass_par_plot,
          pfert_nmass_mat_plot, pfert_nmass_ai_plot, pfert_nmass_par_plot,
          npfert_nmass_mat_plot, npfert_nmass_ai_plot, npfert_nmass_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, hjust = 0, vjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run Pmass models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Pmass - full model
nfert_pmass_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_p_mass"))
summary(nfert_pmass_fullModel)

# Nfert: Pmass - MAT plot
nfert_pmass_mat_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-1, 0.6), breaks = seq(-1, 0.6, 0.4)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# Nfert: Pmass - AI plot
nfert_pmass_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 0.6), breaks = seq(-1, 0.6, 0.4)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Pmass - PAR plot
nfert_pmass_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1, 0.6), breaks = seq(-1, 0.6, 0.4)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Pmass - full model
pfert_pmass_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_p_mass"))
summary(pfert_pmass_fullModel)

# Pfert: Pmass - MAT plot
pfert_pmass_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_mass" & experiment_type == "field"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# Pfert: Pmass - AI plot
pfert_pmass_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Pmass - PAR plot
pfert_pmass_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))


# N+P addition effect on Pmass - full model
npfert_pmass_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "leaf_p_mass"))
summary(npfert_pmass_fullModel)

# N+P fert: Pmass - MAT plot
npfert_pmass_mat_plot <- mod_results(npfert_pmass_fullModel,
                                     mod = "mat", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "magenta", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# N+P fert: Pmass - AI plot
npfert_pmass_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(500, 1100), breaks = seq(500, 1100, 200)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P fert: Pmass - PAR plot
npfert_pmass_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_Pmass_clim_moderators.png", units = "in",
    height = 14, width = 16, res = 600)
ggarrange(nfert_pmass_mat_plot, nfert_pmass_ai_plot, nfert_pmass_par_plot,
          pfert_pmass_mat_plot, pfert_pmass_ai_plot, pfert_pmass_par_plot,
          npfert_pmass_mat_plot, npfert_pmass_ai_plot, npfert_pmass_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, hjust = 0, vjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run leaf_np models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Pmass - full model
nfert_leafnp_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_np"))
summary(nfert_leafnp_fullModel)

# Nfert: Leaf N:P - MAT plot
nfert_leafnp_mat_plot <- mod_results(nfert_leafnp_fullModel,
                                     mod = "mat", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_np"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-5, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = "ln RR of leaf N:P to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Nfert: Leaf N:P - AI plot
nfert_leafnp_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_np"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = "Aridity Index",
       y = "ln RR of leaf N:P to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Leaf N:P - PAR plot
nfert_leafnp_par_plot <- mod_results(nfert_leafnp_fullModel,
                                     mod = "par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_np"),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = "ln RR of leaf N:P to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on leaf N:P - full model
pfert_leafnp_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = pfert_lnRR %>% 
                                   filter(myvar == "leaf_np"))
summary(pfert_leafnp_fullModel)

# Pfert: Leaf N:P - MAT plot
pfert_leafnp_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_np"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 1)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = "ln RR of leaf N:P to P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: Leaf N:P - AI plot
pfert_leafnp_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_np"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 1)) +
  labs(x = "Aridity Index",
       y = "ln RR of leaf N:P to P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Leaf N:P - PAR plot
pfert_leafnp_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_np"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 1)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("ln RR of leaf N:P to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))


# N+P addition effect on leaf N:P - full model
npfert_leafnp_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "leaf_np"))
summary(npfert_leafnp_fullModel)

# N+P fert: Leaf N:P - MAT plot
npfert_leafnp_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_np"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-1.5, 0.6), breaks = seq(-1.5, 0.5, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = "ln RR of leaf N:P to N+P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# N+P fert: Leaf N:P - AI plot
npfert_leafnp_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_np"),
             aes(x = ai, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1.5, 0.6), breaks = seq(-1.5, 0.5, 0.5)) +
  labs(x = "Aridity Index",
       y = "ln RR of leaf N:P to N+P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P fert: Leaf N:P - PAR plot
npfert_leafnp_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_np"),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1.5, 0.6), breaks = seq(-1.5, 0.5, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Interaction effect on Leaf N:P - full model
int_leafnp_fullModel <- rma.mv(yi = dNPi,
                               V = dvNPi,
                               W = dwNPi,
                               method = "REML", 
                               random = ~ 1 | exp, 
                               mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                               slab = exp, control = list(stepadj = 0.3), 
                               data = CNP_effect_sizes_reduced %>% 
                                 filter(myvar == "leaf_np"))
summary(int_leafnp_fullModel)

# Interaction: Pmass - MAT plot
interaction_leafnp_mat_plot <- mod_results(int_leafnp_fullModel,
                                           mod = "mat",
                                           group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(CNP_effect_sizes_reduced, myvar == "leaf_np"),
             aes(x = mat, y = dNPi, size = 1/dvNPi), alpha = 0.75) +
  geom_ribbon(aes(ymin = lowerCL, ymax = upperCL),
              alpha = 0.3, fill = "black") +
  geom_smooth(method = "loess", color = "black", linewidth = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = "Interaction effect size of leaf N:P",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Interaction: Pmass - AI plot
interaction_leafnp_ai_plot <- ggplot() +
  geom_point(data = subset(CNP_effect_sizes_reduced, myvar == "leaf_np"),
             aes(x = ai, y = dNPi, size = 1/dvNPi), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +  
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  labs(x = "Aridity Index",
       y = "Interaction effect size of leaf N:P",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Interaction: Pmass - PAR plot
interaction_leafnp_par_plot <- mod_results(int_pmass_fullModel, 
                                          mod = "par", 
                                          group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(CNP_effect_sizes_reduced, myvar == "leaf_np"),
             aes(x = par, y = dNPi, size = 1/dvNPi), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "black") +
  geom_smooth(method = "loess", linewidth = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = "Interaction effect size of leaf N:P",       
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
# png("../plots/CNPmeta_leafNP_clim_moderators.png", units = "in",
#     height = 20, width = 20, res = 600)
ggarrange(nfert_leafnp_mat_plot, nfert_leafnp_ai_plot, nfert_leafnp_par_plot,
          pfert_leafnp_mat_plot, pfert_leafnp_ai_plot, pfert_leafnp_par_plot,
          npfert_leafnp_mat_plot, npfert_leafnp_ai_plot, npfert_leafnp_par_plot,
          interaction_leafnp_mat_plot, interaction_leafnp_ai_plot, interaction_leafnp_par_plot,
          nrow = 4, ncol = 3, common.legend = TRUE, hjust = 0, vjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)",
                     "(j)", "(k)", "(l)"),
          font.label = list(size = 18))
# dev.off()

###############################################################################
# Re-run Narea models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Narea - full model
nfert_narea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_n_area"))
summary(nfert_narea_fullModel)

# Nfert: Narea - MAT plot
nfert_narea_mat_plot <- mod_results(nfert_narea_fullModel,
                                    mod = "mat", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_area"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of N"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Nfert: Narea - AI plot
nfert_narea_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of N"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Narea - PAR plot
nfert_narea_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_area"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("ln RR of N"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Narea - full model
pfert_narea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_n_area"))
summary(pfert_narea_fullModel)

# Pfert: Narea - MAT plot
pfert_narea_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_area"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of N"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: Narea - AI plot
pfert_narea_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of N"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Narea - PAR plot
pfert_narea_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_area"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("ln RR of N"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P addition effect on Narea - full model
npfert_narea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = npfert_lnRR %>% 
                                  filter(myvar == "leaf_n_area"))
summary(npfert_narea_fullModel)

# N+Pfert: Narea - MAT plot
npfert_narea_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_area"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of N"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Narea - AI plot
npfert_narea_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of N"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Narea - PAR plot
npfert_narea_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_area"),
             aes(x = par, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("ln RR of N"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_Narea_clim_moderators.png", units = "in",
    height = 14, width = 16, res = 600)
ggarrange(nfert_narea_mat_plot, nfert_narea_ai_plot, nfert_narea_par_plot,
          pfert_narea_mat_plot, pfert_narea_ai_plot, pfert_narea_par_plot,
          npfert_narea_mat_plot, npfert_narea_ai_plot, npfert_narea_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, hjust = 0, vjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run Parea models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Parea - full model
nfert_parea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_p_area"))
summary(nfert_parea_fullModel)

# Nfert: Parea - MAT plot
nfert_parea_mat_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_area"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# NPfert: Parea - AI plot
nfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Parea - PAR plot
nfert_parea_par_plot <- mod_results(nfert_parea_fullModel,
                                    mod = "par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_area"),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Parea - full model
pfert_parea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_p_area" & logr > -0.5))
summary(pfert_parea_fullModel)

# Pfert: Parea - MAT plot
pfert_parea_mat_plot <- mod_results(pfert_parea_fullModel,
                                    mod = "mat", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area" & logr > -0.5),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: Parea - AI plot
pfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Parea - PAR plot
pfert_parea_par_plot <- mod_results(pfert_parea_fullModel,
                                    mod = "par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area" & logr > -0.5),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.6, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P addition effect on Parea - full model
npfert_parea_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "leaf_p_area" & logr > -0.5))
summary(npfert_parea_fullModel)

# N+Pfert: Parea - MAT plot
npfert_parea_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area" & logr > -0.5),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Parea - AI plot
npfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Parea - PAR plot
npfert_parea_par_plot <- mod_results(npfert_parea_fullModel,
                                    mod = "par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area"),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.6, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_Parea_clim_moderators.png", units = "in",
    height = 14, width = 16, res = 600)
ggarrange(nfert_parea_mat_plot, nfert_parea_ai_plot, nfert_parea_par_plot,
          pfert_parea_mat_plot, pfert_parea_ai_plot, pfert_parea_par_plot,
          npfert_parea_mat_plot, npfert_parea_ai_plot, npfert_parea_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, hjust = 0, vjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run Marea models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Marea - full model
nfert_marea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "sla" & logr > -0.5))
summary(nfert_marea_fullModel)

# Nfert: Parea - MAT plot
ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "sla"),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# NPfert: Parea - AI plot
nfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Parea - PAR plot
nfert_parea_par_plot <- mod_results(nfert_parea_fullModel,
                                    mod = "par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_area"),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Parea - full model
pfert_parea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_p_area" & logr > -0.5))
summary(pfert_parea_fullModel)

# Pfert: Parea - MAT plot
pfert_parea_mat_plot <- mod_results(pfert_parea_fullModel,
                                    mod = "mat", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area" & logr > -0.5),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: Parea - AI plot
pfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Parea - PAR plot
pfert_parea_par_plot <- mod_results(pfert_parea_fullModel,
                                    mod = "par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area" & logr > -0.5),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.6, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P addition effect on Parea - full model
npfert_parea_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ mat + ai + par + photo_path + myc_assoc + n_fixer,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "leaf_p_area" & logr > -0.5))
summary(npfert_parea_fullModel)

# N+Pfert: Parea - MAT plot
npfert_parea_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area" & logr > -0.5),
             aes(x = mat, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("Mean annual temperature ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Parea - AI plot
npfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area"),
             aes(x = ai, y = logr, size = 1/logr_se), 
             alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = "Aridity Index",
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Parea - PAR plot
npfert_parea_par_plot <- mod_results(npfert_parea_fullModel,
                                     mod = "par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area"),
             aes(x = par, y = logr, size = 1/logr_se), alpha = 0.75) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.6, 2), breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_Parea_clim_moderators.png", units = "in",
    height = 14, width = 16, res = 600)
ggarrange(nfert_parea_mat_plot, nfert_parea_ai_plot, nfert_parea_par_plot,
          pfert_parea_mat_plot, pfert_parea_ai_plot, pfert_parea_par_plot,
          npfert_parea_mat_plot, npfert_parea_ai_plot, npfert_parea_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, hjust = 0, vjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()
