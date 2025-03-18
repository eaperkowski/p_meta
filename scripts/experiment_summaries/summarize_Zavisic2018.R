# determine Zavisic et al. (2018) treatment mean +/- SE
# for photosynthetic data

# Libraries
library(tidyverse)

# Read data
zavisic_data <- read.csv("../raw_data/Zavisic_2018_photosynthesis_cleaned.csv")
head(zavisic_data)

# Generate treatment summary statistics
zavisic_data_summary <- zavisic_data %>%
  group_by(site, treatment) %>%
  summarize(
    
    cica_n = sum(!is.na(ci.ca)),
    cica_mean = mean(ci.ca, na.rm = TRUE),
    cica_sd = sd(ci.ca, na.rm = TRUE),
    cica_se = cica_sd / sqrt(cica_n),
    
    E_n = sum(!is.na(E)),
    E_mean = mean(E, na.rm = TRUE),
    E_sd = sd(E, na.rm = TRUE),
    E_se = E_sd / sqrt(E_n),
    
    
    gs_n = sum(!is.na(gs)),
    gs_mean = mean(gs, na.rm = TRUE),
    gs_sd = sd(gs, na.rm = TRUE),
    gs_se = gs_sd / sqrt(gs_n),
    
    asat_n = sum(!is.na(anet)),
    asat_mean = mean(anet, na.rm = TRUE),
    asat_sd = sd(anet, na.rm = TRUE),
    asat_se = asat_sd / sqrt(asat_n),
    
    iwue_n = sum(!is.na(iwue)),
    iwue_mean = mean(iwue, na.rm = TRUE),
    iwue_sd = sd(iwue, na.rm = TRUE),
    iwue_se = iwue_sd / sqrt(iwue_n))


# Prep for easy merge into compiled datasheet
zavisic_data_summary_control <- zavisic_data_summary %>%
  filter(treatment == "C") %>%
  select(-treatment)
names(zavisic_data_summary_control)[2:21] <- str_c(names(zavisic_data_summary_control)[2:21], "_control")

zavisic_data_summary_treatment <- zavisic_data_summary %>%
  filter(treatment != "C")
names(zavisic_data_summary_treatment)[3:22] <- str_c(names(zavisic_data_summary_treatment)[3:22], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
zavisic_data_summary_control %>%
  full_join(zavisic_data_summary_treatment, by = "site") %>%
  dplyr::select(site, treatment,
                
                asat_mean_control, asat_mean_trt, 
                asat_sd_control, asat_sd_trt, 
                asat_se_control, asat_se_trt, 
                asat_n_control, asat_n_trt,
                
                E_mean_control, E_mean_trt, 
                E_sd_control, E_sd_trt, 
                E_se_control, E_se_trt, 
                E_n_control, E_n_trt,
                
                cica_mean_control, cica_mean_trt, 
                cica_sd_control, cica_sd_trt, 
                cica_se_control, cica_se_trt, 
                cica_n_control, cica_n_trt,
                
                gs_mean_control, gs_mean_trt, 
                gs_sd_control, gs_sd_trt, 
                gs_se_control, gs_se_trt, 
                gs_n_control, gs_n_trt,
                
                iwue_mean_control, iwue_mean_trt, 
                iwue_sd_control, iwue_sd_trt, 
                iwue_se_control, iwue_se_trt, 
                iwue_n_control, iwue_n_trt) %>%
  pivot_longer(cols = asat_mean_control:iwue_n_trt,
               names_to = c("trait", "stat", "trt"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  mutate(trait = factor(trait, levels = c("asat", "E", "cica", "gs", "iwue"))) %>%
  arrange(site, trait) %>%
  write.csv("../data_summaries/Zavisic2018_summary.csv", row.names = F)
                

