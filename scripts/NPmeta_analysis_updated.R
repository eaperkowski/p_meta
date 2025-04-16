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
library(naniar) # to resolve NA/<NA> issue

# Read data sources (MESI, NutNet, EAP manual compilation)
mesi <- read.csv("../data/mesi_main_manual.csv")
nutnet <- read.csv("../data/nutnet_main.csv")
eap <- read.csv("../data/eap_main.csv")

# Merge data sources into single data frame
full_df <- mesi %>% full_join(nutnet) %>% full_join(eap) %>%
  replace_with_na_all(~.x == "<NA>") %>%
  dplyr::select(-doi, -x_units, -se_c, -se_t) %>%
  mutate(fert = ifelse(npk == "_100", 
                       "n", ifelse(npk == "_010",
                                   "p", ifelse(npk == "_110", "np", NA))))

# Create experiment metadata summary
experiment_summary <- full_df %>%
  dplyr::select(citation, exp:experiment_type) %>%
  distinct(citation, exp, .keep_all = TRUE)

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

##############################################################################
# Some prep work for calculating nteraction effect sizes (need
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
           sep = "XXX", remove = FALSE) %>%
  full_join(experiment_summary, by = c("citation", "exp")) %>%
  dplyr::select(citation, exp, latitude:experiment_type, 
                species, response,
                
                # Individual effects
                gNi = g_a, vNi = v_a, wNi = w_a,
                gPi = g_b, vPi = v_b, wPi = w_b,
                gNPi = g_ab, vNPi = v_ab, wNPi = w_ab,
                
                # Main effects
                dNi = dA, dvNi = v_a_main, dwNi = w_a_main,
                dPi = dB, dvPi = v_b_main, dwPi = w_b_main,
                
                # Interaction effects
                dNPi = dAB, dvNPi = v_ab_int, dwNPi = w_ab_int) %>%
  replace_with_na_all(~.x == "none") 

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
           myvar != "anet_mass")
unique(CNP_effect_sizes_reduced$myvar)

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
      summarise(intES_min = min(dNPi), logr_max = max(dNPi)) |> 
      rename(var = myvar),
    by = "var")

###############################################################################
# Let' make some plots!
###############################################################################
CNP_effect_sizes_reduced <- CNP_effect_sizes_reduced %>%
  mutate(myvar = factor(myvar, 
               levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                          "agb_n", "agb", "total_biomass", 
                          "leaf_wue", "leaf_ppue", "leaf_pnue", 
                          "jmax", "vcmax", "asat", "leaf_residual_p", 
                          "leaf_structural_p", "leaf_nucleic_p", 
                          "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi", 
                          "leaf_np", "leaf_p_area", "leaf_p_mass", 
                          "leaf_n_area", "leaf_n_mass", "lma")))

df_box_int <- df_box_int %>%
  mutate(var = factor(var, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_n", "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma")))


ggplot(
  data = CNP_effect_sizes_reduced %>%
    filter(myvar %in% c("lma", "leaf_n_mass", "leaf_n_area", 
                        "leaf_p_mass", "leaf_p_area", "leaf_np")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = dwNPi)) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("lma", "leaf_n_mass", "leaf_n_area", 
                                    "leaf_p_mass", "leaf_p_area", "leaf_np")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.4, fill = "blue",
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]), 
                              expression("M"["area"]))) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1)) +
  labs(x = "", 
       y = "Interaction effect size (Hedge's d)",
       size = expression(bold("Error"^"-1"))) +
  scale_size(range = c(0.25, 4)) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))

# Plot leaf phosphorus fractionation. Separating by trait type to avoid plot overwhelm
meta_Intplot_phosFract <- ggplot(
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
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.4, fill = "blue",
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
       y = "Interaction effect size (Hedge's d)",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_Intplot_phosFract

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
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.4, fill = "blue",
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
       size = expression(bold("Error"^"-1"))) +
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
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6, fill = "blue",
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
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_intPlot_biomass


###############################################################################
# Let's play with some bivariate relationships
###############################################################################

ggplot(data = CNP_effect_sizes_reduced,
       aes(x = ))

