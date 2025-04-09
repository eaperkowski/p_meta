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

# Load effect size function (calculates individual, main, and interaction
# effect sizes)
source("../helper_fxns/calc_intxn_effSize_meta.R")

# Merge data sources into single data frame, then convert dataset to
# wide format to make for easy determination of individual, main, and
# interaction effect sizes
full_df <- mesi %>% full_join(nutnet) %>% full_join(eap) %>%
  replace_with_na_all(~.x == "<NA>") %>%
  dplyr::select(-doi, -x_units, -se_c, -se_t) %>%
  mutate(fert = ifelse(npk == "_100", 
                       "n", ifelse(npk == "_010",
                                   "p", ifelse(npk == "_110", "np", NA))))

head(full_df)
full_df_control <- full_df %>%
  mutate(unique_id = ifelse(is.na(species), 
                            str_c(citation, "XXX", site, "XXX", response, "XXX", "none"),
                            str_c(citation, "XXX", site, "XXX", response, "XXX", species))) %>%
  group_by(unique_id) %>%
  summarize(x_c = mean(x_c),
            sd_c = mean(sd_c),
            rep_c = mean(rep_c))

full_df_trt_long <- full_df %>% 
  mutate(unique_id = ifelse(is.na(species), 
                            str_c(citation, "XXX", site, "XXX", response, "XXX", "none"),
                            str_c(citation, "XXX", site, "XXX", response, "XXX", species))) %>%
  dplyr::select(unique_id, fert, x_t, sd_t, rep_t) %>%
  pivot_wider(names_from = fert, values_from = c(x_t, sd_t, rep_t)) %>%
  full_join(full_df_control) %>%
  dplyr::select(unique_id,
                x_c, sd_c, rep_c,
                x_n = x_t_n, sd_n = sd_t_n, rep_n = rep_t_n,
                x_p = x_t_p, sd_p = sd_t_p, rep_p = rep_t_p,
                x_np = x_t_np, sd_np = sd_t_np, rep_np = rep_t_np)
 
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
                        n_ab = full_df_trt_long$rep_np)

  





