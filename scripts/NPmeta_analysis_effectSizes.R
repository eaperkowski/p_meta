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

# Load effect size function (calculates individual, main, and interaction
# effect sizes)
source("../helper_fxns/calc_intxn_effSize_meta.R")

# Merge data sources into single data frame, then convert dataset to
# wide format to make for easy determination of individual, main, and
# interaction effect sizes
full_df <- mesi %>% full_join(nutnet) %>% full_join(eap) %>%
  dplyr::select(-doi, -x_units, -se_c, -se_t) %>%
  mutate(fert = ifelse(npk == "_100", 
                       "n", ifelse(npk == "_010",
                                   "p", ifelse(npk == "_110", "np", NA))))

head(full_df)
full_df_control <- full_df %>%
  group_by(exp, response, species) %>%
  summarize(x_c = mean(x_c))




%>%
  pivot_wider(names_from = fert,
              values_from = x_c:rep_t)

head(full_df)





