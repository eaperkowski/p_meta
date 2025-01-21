# summarize Cleland et al. (2019)

# Libraries
library(tidyverse)

# Load NutNet data from Cleland et al. (2019), rename some of the colnames
cleland_data <- read.csv("../raw_data/Cleland_2019_ecosystems.csv") %>%
  filter(Kµ == 0) %>%
  rename(bgb = Rootsgperm2, rmf = rootmassfraction) %>%
  mutate(agb = (bgb/rmf) - bgb,
         root_shoot_ratio = bgb / agb,
         treatment_name = factor(treatment_name, levels = c("Control", "N", "P", "NP")))

# How many sites?
unique(cleland_data$site_code)

# Determine summary statistics for Cleland et al. (2019)
cleland_data_summary <- cleland_data %>%
  group_by(site_code, treatment_name) %>%
  summarize(
    agb_n = sum(!is.na(agb)),
    agb_mean = mean(agb, na.rm = TRUE),
    agb_sd = sd(agb, na.rm = TRUE),
    agb_se = agb_sd / sqrt(agb_n),
    
    bgb_n = sum(!is.na(bgb)),
    bgb_mean = mean(bgb, na.rm = TRUE),
    bgb_sd = sd(bgb, na.rm = TRUE),
    bgb_se = bgb_sd / sqrt(bgb_n),
    
    rmf_n = sum(!is.na(rmf)),
    rmf_mean = mean(rmf, na.rm = TRUE),
    rmf_sd = sd(rmf, na.rm = TRUE),
    rmf_se = rmf_sd / sqrt(rmf_n),
    
    rootshoot_n = sum(!is.na(root_shoot_ratio)),
    rootshoot_mean = mean(root_shoot_ratio, na.rm = TRUE),
    rootshoot_sd = sd(root_shoot_ratio, na.rm = TRUE),
    rootshoot_se = rootshoot_sd / sqrt(rootshoot_n))

# Prep for easy merge into compiled datasheet
cleland_data_summary_control <- cleland_data_summary %>%
  filter(treatment_name == "Control") %>%
  select(-treatment_name)
names(cleland_data_summary_control)[2:17] <- str_c(names(cleland_data_summary_control)[2:17], "_control")

cleland_data_summary_treatment <- cleland_data_summary %>%
  filter(treatment_name != "Control")
names(cleland_data_summary_treatment)[3:18] <- str_c(names(cleland_data_summary_treatment)[3:18], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
cleland_data_summary_control %>%
  full_join(cleland_data_summary_treatment, by = "site_code") %>%
  dplyr::select(site_code, treatment_name, 
                
                agb_mean_control, agb_mean_trt, agb_sd_control, agb_sd_trt, 
                agb_se_control, agb_se_trt, agb_n_control, agb_n_trt,
                
                bgb_mean_control, bgb_mean_trt, bgb_sd_control, bgb_sd_trt, 
                bgb_se_control, bgb_se_trt, bgb_n_control, bgb_n_trt,
                
                rmf_mean_control, rmf_mean_trt, rmf_sd_control, rmf_sd_trt, 
                rmf_se_control, rmf_se_trt, rmf_n_control, rmf_n_trt,
                
                rootshoot_mean_control, rootshoot_mean_trt, rootshoot_sd_control, 
                rootshoot_sd_trt, rootshoot_se_control, rootshoot_se_trt, 
                rootshoot_n_control, rootshoot_n_trt) %>%
  
  pivot_longer(cols = agb_mean_control:rootshoot_n_trt,
               names_to = c("trait", "stat", "trt"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  arrange(trait, site_code, treatment_name) %>%
  write.csv("../data_summaries/Cleland2019_summary.csv", row.names = F)
