# determine Mayor et al. 2013 treatment mean +/- SE

# Libraries
library(tidyverse)

# Calculate ci:ci from d13C isotopes
source("../../../r_functions/calc_chi.R")

# Load dataset
mayor2013 <- read.csv("../../raw_data/Mayor_2013.csv") %>%
  dplyr::select(Plot, trt = Rx, N:Leaf.15N) %>%
  mutate(Leaf.NP = Leaf.N / Leaf.P,
         cica = calc_chi_c3(leaf.d13c = Leaf.13C, year = 2011)$chi)
head(mayor2013)

# Summarize dataset
mayor2013_data_summary_plot <- mayor2013 %>%
  group_by(Plot, trt, Species) %>%
  summarize(leaf.N = mean(Leaf.N, na.rm = TRUE),
            leaf.P = mean(Leaf.P, na.rm = TRUE),
            leaf.Pi = mean(Leaf.Pi, na.rm = TRUE),
            leaf.NP = mean(Leaf.NP, na.rm = TRUE),
            cica = mean(cica, na.rm = TRUE)) %>%
  arrange(Species, Plot, trt)
mayor2013_data_summary_plot  

# Data summary by plot mean
mayor2013_data_summary <- mayor2013_data_summary_plot %>%
  group_by(trt, Species) %>%
  summarize(
    
    nmass_n = sum(!is.na(leaf.N)),
    nmass_mean = mean(leaf.N, na.rm = TRUE),
    nmass_sd = sd(leaf.N, na.rm = TRUE),
    nmass_se = nmass_sd / sqrt(nmass_n),
    
    pmass_n = sum(!is.na(leaf.P)),
    pmass_mean = mean(leaf.P, na.rm = TRUE),
    pmass_sd = sd(leaf.P, na.rm = TRUE),
    pmass_se = pmass_sd / sqrt(pmass_n),
    
    leafPi_n = sum(!is.na(leaf.Pi)),
    leafPi_mean = mean(leaf.Pi, na.rm = TRUE),
    leafPi_sd = sd(leaf.Pi, na.rm = TRUE),
    leafPi_se = leafPi_sd / sqrt(leafPi_n),
    
    leafnp_n = sum(!is.na(leaf.NP)),
    leafnp_mean = mean(leaf.NP, na.rm = TRUE),
    leafnp_sd = sd(leaf.NP, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(leafnp_n),
    
    cica_n = sum(!is.na(cica)),
    cica_mean = mean(cica, na.rm = TRUE),
    cica_sd = sd(cica, na.rm = TRUE),
    cica_se = cica_sd / sqrt(cica_n)
    
    )

# Prep for easy merge into compiled datasheet
mayor2013_data_summary_control <- mayor2013_data_summary %>%
  filter(trt == "C") %>%
  ungroup(trt) %>%
  select(-trt)
names(mayor2013_data_summary_control)[2:21] <- str_c(names(mayor2013_data_summary_control)[2:21], "_control")

mayor2013_data_summary_treatment <- mayor2013_data_summary %>%
  filter(trt != "C")
names(mayor2013_data_summary_treatment)[3:22] <- str_c(names(mayor2013_data_summary_treatment)[3:22], "_trt")


# Format into easy merge into compiled datasheet, write to .csv
mayor2013_data_summary_control %>%
  full_join(mayor2013_data_summary_treatment, by = "Species") %>%
  dplyr::select(Species, trt,
                
                nmass_mean_control, nmass_mean_trt, 
                nmass_sd_control, nmass_sd_trt, 
                nmass_se_control, nmass_se_trt, 
                nmass_n_control, nmass_n_trt,
                
                pmass_mean_control, pmass_mean_trt, 
                pmass_sd_control, pmass_sd_trt, 
                pmass_se_control, pmass_se_trt, 
                pmass_n_control, pmass_n_trt,
                
                leafPi_mean_control, leafPi_mean_trt, 
                leafPi_sd_control, leafPi_sd_trt, 
                leafPi_se_control, leafPi_se_trt, 
                leafPi_n_control, leafPi_n_trt,
                
                leafnp_mean_control, leafnp_mean_trt, 
                leafnp_sd_control, leafnp_sd_trt, 
                leafnp_se_control, leafnp_se_trt, 
                leafnp_n_control, leafnp_n_trt,
                
                cica_mean_control, cica_mean_trt, 
                cica_sd_control, cica_sd_trt, 
                cica_se_control, cica_se_trt, 
                cica_n_control, cica_n_trt
                
                
                ) %>%
  pivot_longer(cols = nmass_mean_control:cica_n_trt,
               names_to = c("trait", "stat", "treatment"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:treatment,
              values_from = value) %>%
  mutate(trait = factor(trait, levels = c("nmass", "pmass", 
                                          "leafPi", "leafnp", "cica")),
         trt = factor(trt, levels = c("N", "P", "NP")),
         Species = factor(Species, levels = c("ALSEBL", "HEISCO",
                                              "TET2PA", "OENOMA"))) %>%
  group_by(Species, trait) %>%
  filter(n_control > 1 & n_trt > 1) %>%
  filter(n_distinct(trt) == 3) %>%
  arrange(trait, Species, trt) %>%
  write.csv("../../data_summaries/species_level/Mayor2013_summarized.csv", row.names = F)







