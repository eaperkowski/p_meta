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

# Helper function to calculate weighted mean
weighted_mean <- function(effect_size, weight) {
  sum(effect_size * weight) / sum(weight)
}

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
                gN = g_a, vN = v_a, wN = w_a,
                gP = g_b, vP = v_b, wP = w_b,
                gNP = g_ab, vNP = v_ab, wNP = w_ab,
                
                # Main effects
                dN = dA, dvN = v_a_main, dwN = w_a_main,
                dP = dB, dvP = v_b_main, dwP = w_b_main,
                
                # Interaction effects
                dNP = dAB, dvNP = v_ab_int, dwNP = w_ab_int) %>%
  replace_with_na_all(~.x == "none") 

###############################################################################
# From effect size, variance, and weight, calculate the weighted means for the
# individual effect size of N addition treatments. Values will be randomly
# resampled to generate upper and lower confidence intervals of effect sizes
###############################################################################
weighted_effSize_prep <- CNP_effect_sizes %>%
  mutate(response = ifelse(response == "amax", 
                      "asat", 
                      response),
         response = ifelse(response == "anet", 
                      "asat", 
                      response),
         response = ifelse(response == "fine_root_biomass" | response == "bnpp", 
                      "bgb", 
                      response),
         response = ifelse(response == "anpp", 
                           "agb", 
                           response),
         response = ifelse(response == "total_leaf_area",
                           "tla",
                           response)) %>%
  filter(response != "r_eco" & response != "nee" & 
           response != "spad" & response != "tpu" &
           response != "rd" & response != "stom_lim" &
           response != "anet_mass")
unique(weighted_effSize_prep$response)


effect_sizes %>%
  group_by(response) %>%
  summarize(
    
    # Individual effects
    wm_gN = weighted_mean(effect_size = gN, weight = wN),
    wm_gP = weighted_mean(effect_size = gP, weight = wP),
    wm_gNP = weighted_mean(effect_size = gNP, weight = wNP),
    
    # Main effects
    wm_dN = weighted_mean(effect_size = dN, weight = dwN),
    wm_dP = weighted_mean(effect_size = dP, weight = dwP),
    
    # Interaction effect
    wm_dNP = weighted_mean(effect_size = dNP, weight = dwNP)
  )
  
  
effect_sizes %>%
  group_by(response) %>%
  summarize(wm = weighted_mean(effect_size = gP, weight = wP)) 
 
effect_sizes %>%
  group_by(response) %>%
  summarize(wm = weighted_mean(effect_size = gP, weight = wP))
  
 






 effect_sizes %>%
  group_by(response) %>%
  group_modify( ~ {
    boot_result <- boot(data = ., statistic = bootstrap_wm, R = 1000)
    ci = boot.ci(boot_result)
    
    tibble(
      weighted_mean = weighted_mean(.x$gN, .x$wN),
      lower_ci = ci$percent[4],
      upper_ci = ci$percent[5]
    )
  })








