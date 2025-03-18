# Reformat Wang et al. (2019) summary table for easy load
# into compiled meta-analysis datasheet

# Library
library(tidyverse)

# Read data
wang_data_summary <- read.csv("../raw_data/Wang_2019_summary_table_from_text.csv")
head(wang_data_summary)

# Prep for easy merge into compiled datasheet
wang_data_summary_control <- wang_data_summary %>%
  filter(treatment == "C") %>%
  mutate(exp = "wang2019") %>%
  select(-treatment) %>%
  pivot_wider(names_from = trait, values_from = mean:n,
              names_glue = "{trait}_{.value}") ## Need to change order of label here to make for easy merge later
names(wang_data_summary_control)[2:21] <- str_c(names(wang_data_summary_control)[2:21], "_control")

wang_data_summary_treatment <- wang_data_summary %>%
  filter(treatment != "C") %>%
  mutate(exp = "wang2019") %>%
  pivot_wider(names_from = trait, values_from = mean:n,
              names_glue = "{trait}_{.value}")
names(wang_data_summary_treatment)[3:22] <- str_c(names(wang_data_summary_treatment)[3:22], "_trt")

# Merge control and trt data frames
wang_data_summary_control %>%
  full_join(wang_data_summary_treatment, by = "exp") %>%
  dplyr::select(exp, treatment, 
                
                pmass_mean_control, pmass_mean_trt, 
                pmass_sd_control, pmass_sd_trt, 
                pmass_se_control, pmass_se_trt, 
                pmass_n_control, pmass_n_trt,
                
                pi_mean_control, pi_mean_trt, 
                pi_sd_control, pi_sd_trt, 
                pi_se_control, pi_se_trt, 
                pi_n_control, pi_n_trt,
                
                psugar_mean_control, psugar_mean_trt, 
                psugar_sd_control, psugar_sd_trt, 
                psugar_se_control, psugar_se_trt, 
                psugar_n_control, psugar_n_trt,
                
                pnucleic_mean_control, pnucleic_mean_trt, 
                pnucleic_sd_control, pnucleic_sd_trt, 
                pnucleic_se_control, pnucleic_se_trt, 
                pnucleic_n_control, pnucleic_n_trt,
                
                presidual_mean_control, presidual_mean_trt, 
                presidual_sd_control, presidual_sd_trt, 
                presidual_se_control, presidual_se_trt, 
                presidual_n_control, presidual_n_trt) %>%
  pivot_longer(cols = pmass_mean_control:presidual_n_trt,
               names_to = c("trait", "stat", "trt"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  mutate(treatment = factor(treatment, levels = c("N1", "N2", "P", 
                                                  "N1P", "N2P")),
         trait = factor(trait, levels = c("pmass", "pi", "psugar",
                                          "pnucleic", "presidual"))) %>%
  arrange(trait, treatment) %>%
  write.csv("../data_summaries/Wang2019_summary.csv", row.names = F)

