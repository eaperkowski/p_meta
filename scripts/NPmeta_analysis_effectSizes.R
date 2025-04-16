# Script that explores the effect of N, P, and N+P addition on leaf- 
# and whole-plant functional traits. Meta-analysis uses data compiled
# from full-factorial N*P manipulation experiments. Whole-plant
# measurements are aggregated at the experiment level and leaf-level
# measurements are aggregated at the species level.

###############################################################################
# Load libraries, data, and helper functions
###############################################################################

# Libraries
library(tidyverse)
library(ggpubr)
library(boot)
library(forcats)
library(patchwork)
library(naniar) # to resolve NA/<NA> issue

# Read data sources (MESI, NutNet, EAP manual compilation)
mesi <- read.csv("../data/mesi_main_manual.csv")
nutnet <- read.csv("../data/nutnet_main.csv")
eap <- read.csv("../data/eap_main.csv")

# Load effect size function (calculates individual, main, and interaction
# effect sizes)
source("../helper_fxns/calc_intxn_effSize_meta.R")
source("../helper_fxns/analyse_meta.R")

###############################################################################
# Merge data sources, convert to wide format for easy effect size calculation
###############################################################################
# Merge data sources into single data frame
full_df <- mesi %>% full_join(nutnet) %>% full_join(eap) %>%
  replace_with_na_all(~.x == "<NA>") %>%
  dplyr::select(-doi, -x_units, -se_c, -se_t) %>%
  mutate(fert = ifelse(npk == "_100", 
                       "n", ifelse(npk == "_010",
                                   "p", ifelse(npk == "_110", "np", NA))))

# Create experiment metadata df, prep for merge after effect size calculation
experiment_summary <- full_df %>%
  dplyr::select(citation, exp:experiment_type) %>%
  distinct(citation, exp, .keep_all = TRUE)


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
  mutate(var = response,
         var = ifelse(var == "amax", "asat", var),
         var = ifelse(var == "anet", "asat", var),
         var = ifelse(var == "fine_root_biomass" | var == "bnpp", 
                           "bgb", var),
         var = ifelse(var == "anpp", "agb", var),
         var = ifelse(var == "total_leaf_area", "tla", var)) %>%
  filter(var != "r_eco" & var != "nee" & 
           var != "spad" & var != "tpu" &
           var != "rd" & var != "stom_lim" &
           var != "anet_mass")
unique(CNP_effect_sizes_reduced$var)

###############################################################################
###############################################################################
# Let's put together some meta-regression models
############################################################################### 
###############################################################################

#############################################
# Aboveground biomass: Individual N effect
#############################################
# Model for individual N effect
agb_model_N <- rma.mv(yi = dNPi,
                      V = vNPi,
                      W = wNPi,
                      random = ~ 1 | exp,
                      slab = exp,
                      data = subset(CNP_effect_sizes_reduced, var == "leaf_residual_p"))
summary(agb_model_N)

# check for publication bias
funnel(agb_model_N)

# Start data frame
summary_results <- data.frame(trait = "agb",
                              type = "ind_n_effect",
                              predict(agb_model_N))

#############################################
# Aboveground biomass: Individual P effect
#############################################

# Model for individual N effect
agb_model <- rma.mv(yi = gNi,
                    V = vNi,
                    W = wNi,
                    random = ~ 1 | exp,
                    slab = exp,
                    data = subset(CNP_effect_sizes_reduced, var == "agb"))
summary(agb_model)

# Start data frame
summary_results <- data.frame(trait = "agb",
                              type = "ind_n_effect",
                              predict(agb_model))







