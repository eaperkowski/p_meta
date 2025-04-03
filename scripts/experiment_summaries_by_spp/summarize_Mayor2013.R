# determine Mayor et al. 2013 treatment mean +/- SE

# Libraries
library(tidyverse)

# Load dataset
mayor2013 <- read.csv("../../raw_data/Mayor_2013.csv") %>%
  dplyr::select(Plot, trt = Rx, N:Leaf.15N) %>%
  mutate(Leaf.NP = Leaf.N / Leaf.P)
head(mayor2013)

# Summarize dataset
mayor2013_data_summary_plot <- mayor2013 %>%
  group_by(Plot, trt, Species) %>%
  summarize(leaf.n_plotMean = mean(Leaf.N, na.rm = TRUE),
            leaf.p_plotMean = mean(Leaf.P, na.rm = TRUE),
            leaf.pi_plotMean = mean(Leaf.Pi, na.rm = TRUE),
            leaf.np_plotMean = mean(Leaf.NP, na.rm = TRUE)) %>%
  arrange(Species, Plot, trt)
mayor2013_data_summary  

# Data summary by plot mean
mayor2013_data_summary <- mayor2013_data_summary_plot %>%
  group_by(trt, Species) %>%
  summarize(
    
    nmass_n = sum(!is.na(leaf.n_plotMean)),
    nmass_mean = mean(leaf.n_plotMean, na.rm = TRUE),
    nmass_sd = sd(leaf.n_plotMean, na.rm = TRUE),
    nmass_se = nmass_sd / sqrt(nmass_n),
    
    pmass_n = sum(!is.na(leaf.p_plotMean)),
    pmass_mean = mean(leaf.p_plotMean, na.rm = TRUE),
    pmass_sd = sd(leaf.p_plotMean, na.rm = TRUE),
    pmass_se = pmass_sd / sqrt(pmass_n),
    
    leafPi_n = sum(!is.na(leaf.pi_plotMean)),
    leafPi_mean = mean(leaf.pi_plotMean, na.rm = TRUE),
    leafPi_sd = sd(leaf.pi_plotMean, na.rm = TRUE),
    leafPi_se = leafPi_sd / sqrt(leafPi_n),
    
    leafnp_n = sum(!is.na(leaf.np_plotMean)),
    leafnp_mean = mean(leaf.np_plotMean, na.rm = TRUE),
    leafnp_sd = sd(leaf.np_plotMean, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(leafnp_n)
    )

# Prep for easy merge into compiled datasheet
mayor2013_data_summary_control <- mayor2013_data_summary %>%
  filter(trt == "C") %>%
  ungroup(trt) %>%
  select(-trt)
names(mayor2013_data_summary_control)[2:17] <- str_c(names(mayor2013_data_summary_control)[2:17], "_control")

mayor2013_data_summary_treatment <- mayor2013_data_summary %>%
  filter(trt != "C")
names(mayor2013_data_summary_treatment)[3:18] <- str_c(names(mayor2013_data_summary_treatment)[3:18], "_trt")


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
                leafnp_n_control, leafnp_n_trt) %>%
  pivot_longer(cols = nmass_mean_control:leafnp_n_trt,
               names_to = c("trait", "stat", "treatment"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:treatment,
              values_from = value) %>%
  mutate(trait = factor(trait, levels = c("nmass", "pmass", 
                                          "leafPi", "leafnp")),
         trt = factor(trt, levels = c("N", "P", "NP")),
         Species = factor(Species, levels = c("ALSEBL", "HEISCO",
                                              "TET2PA", "OENOMA"))) %>%
  group_by(Species, trait) %>%
  filter(n_control > 1 & n_trt > 1) %>%
  filter(n_distinct(trt) == 3) %>%
  arrange(trait, Species, trt) %>%
  write.csv("../../data_summaries/species_level/Mayor2013_summarized.csv", row.names = F)







