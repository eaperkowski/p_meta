# determine Nutnet treatment mean +/- SE

# Libraries
library(tidyverse)

# Load NutNet data from Hersch-Green et al. (2024), rename some of the colnames
herschgreen_data <- read.csv("../../raw_data/Hersch-Green-Petosky_metabolicMS_EDI.csv") %>%
  dplyr::select(Site:MAP_v2, Treatment, GS.pg, Cmass = Cmg.mgleaftissue, 
                Nmass = Nmg.mgleaftissue, Pmass = Pmg.mgleaftissue, 
                Amax = Amax_.µmolCO2m..s..., E = E_.mmolH2Om..s...,
                iWUE = Instant_WUE.µmolCO2mmolH2O..., total_biomass = DryMass2_mg) %>%
  mutate(leaf_np = Nmass / Pmass)
head(herschgreen_data)

# How many sites?
unique(herschgreen_data$Site)

# Calculate site-level treatment means +/- SE
herschgreen_data_summary <- herschgreen_data %>%
  group_by(Site, Treatment, Taxa) %>%
  mutate(Treatment = factor(Treatment, levels = c("C", "N", "P", "NP"))) %>%
  summarize(
    Nmass_n = sum(!is.na(Nmass)),
    Nmass_mean = mean(Nmass, na.rm = TRUE),
    Nmass_sd = sd(Nmass, na.rm = TRUE),
    Nmass_se = Nmass_sd / sqrt(Nmass_n),
    
    Pmass_n = sum(!is.na(Pmass)),
    Pmass_mean = mean(Pmass, na.rm = TRUE),
    Pmass_sd = sd(Pmass, na.rm = TRUE),
    Pmass_se = Pmass_sd / sqrt(Pmass_n),
    
    leafnp_n = sum(!is.na(leaf_np)),
    leafnp_mean = mean(leaf_np, na.rm = TRUE),
    leafnp_sd = sd(leaf_np, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(leafnp_n),
    
    Amax_n = sum(!is.na((Amax))),
    Amax_mean = mean(Amax, na.rm = TRUE),
    Amax_sd = sd(Amax, na.rm = TRUE),
    Amax_se = Amax_sd / sqrt(Amax_n),
    
    E_n = sum(!is.na((E))),
    E_mean = mean(E, na.rm = TRUE),
    E_sd = sd(E, na.rm = TRUE),
    E_se = E_sd / sqrt(E_n),
    
    iWUE_n = sum(!is.na((iWUE))),
    iWUE_mean = mean(iWUE, na.rm = TRUE),
    iWUE_sd = sd(iWUE, na.rm = TRUE),
    iWUE_se = iWUE_sd / sqrt(iWUE_n))


# Prep for easy merge into compiled datasheet
herschgreen_data_summary_control <- herschgreen_data_summary %>%
  filter(Treatment == "C") %>%
  ungroup(Treatment) %>%
  select(-Treatment)
names(herschgreen_data_summary_control)[3:26] <- str_c(names(herschgreen_data_summary_control)[3:26], "_control")

herschgreen_data_summary_treatment <- herschgreen_data_summary %>%
  filter(Treatment != "C")
names(herschgreen_data_summary_treatment)[4:27] <- str_c(names(herschgreen_data_summary_treatment)[4:27], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
herschgreen_data_summary_control %>%
  full_join(herschgreen_data_summary_treatment, by = c("Site", "Taxa")) %>%
  dplyr::select(Site, Treatment, Taxa,
                
                Nmass_mean_control, Nmass_mean_trt,
                Nmass_sd_control, Nmass_sd_trt,
                Nmass_se_control, Nmass_se_trt,
                Nmass_n_control, Nmass_n_trt,
                
                Pmass_mean_control, Pmass_mean_trt, 
                Pmass_sd_control, Pmass_sd_trt, 
                Pmass_se_control, Pmass_se_trt, 
                Pmass_n_control, Pmass_n_trt,
                
                leafnp_mean_control, leafnp_mean_trt, 
                leafnp_sd_control, leafnp_sd_trt, 
                leafnp_se_control, leafnp_se_trt, 
                leafnp_n_control, leafnp_n_trt,
                
                Amax_mean_control, Amax_mean_trt, 
                Amax_sd_control, Amax_sd_trt, 
                Amax_se_control, Amax_se_trt, 
                Amax_n_control, Amax_n_trt,
                
                E_mean_control, E_mean_trt, 
                E_sd_control, E_sd_trt, 
                E_se_control, E_se_trt, 
                E_n_control, E_n_trt,
                
                iWUE_mean_control, iWUE_mean_trt, 
                iWUE_sd_control, iWUE_sd_trt, 
                iWUE_se_control, iWUE_se_trt, 
                iWUE_n_control, iWUE_n_trt) %>%
  pivot_longer(cols = Nmass_mean_control:iWUE_n_trt,
             names_to = c("trait", "stat", "trt"), 
             names_sep = "_",
             values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  group_by(Site, Taxa, trait) %>%
  filter(n_control > 1 & n_trt > 1) %>%
  filter(n_distinct(Treatment) == 3) %>%
  arrange(trait, Site, Taxa, Treatment) %>%
  write.csv("../../data_summaries/HerschGreen2024_summary_spp.csv", row.names = F)

  



