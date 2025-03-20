# determine Fisher et al (2013) treatment mean +/- SE

# Libraries
library(tidyverse)

# Load fisher et al. 2013 data, then calculate summary statistics
fisher_data_summary <- read.csv("../raw_data/Fisher2013_raw.csv") %>%
  filter(fert != "") %>%
  mutate(leaf_np = round(leafn_percent / leafp_percent, digits = 2),
         site = factor(site, c("TA", "TO", "SP", "WA")),
         fert = factor(fert, c("C", "N0", "P0", "NP0"))) %>%
  group_by(site, fert) %>%
  summarize(nmass_n = length(!is.na(leafn_percent)),
            nmass_mean = mean(leafn_percent, na.rm = TRUE),
            nmass_sd = sd(leafn_percent, na.rm = TRUE),
            nmass_se = nmass_sd / sqrt(nmass_n),
            
            pmass_n = length(!is.na(leafp_percent)),
            pmass_mean = mean(leafp_percent, na.rm = TRUE),
            pmass_sd = sd(leafp_percent, na.rm = TRUE),
            pmass_se = pmass_sd / sqrt(pmass_n),
            
            leafnp_n = length(!is.na(leaf_np)),
            leafnp_mean = mean(leaf_np, na.rm = TRUE),
            leafnp_sd = sd(leaf_np, na.rm = TRUE),
            leafnp_se = leafnp_sd / sqrt(leafnp_n)
            )

# Prep for easy merge into compiled datasheet
fisher_data_summary_control <- fisher_data_summary %>%
  filter(fert == "C") %>%
  select(-fert)
names(fisher_data_summary_control)[2:13] <- str_c(names(fisher_data_summary_control)[2:13], "_control")

fisher_data_summary_treatment <- fisher_data_summary %>%
  filter(fert != "C")
names(fisher_data_summary_treatment)[3:14] <- str_c(names(fisher_data_summary_treatment)[3:14], "_trt")


fisher_data_summary_control %>%
  full_join(fisher_data_summary_treatment, by = c("site")) %>%
  dplyr::select(site, fert,
                
                nmass_mean_control, nmass_mean_trt, 
                nmass_sd_control, nmass_sd_trt, 
                nmass_se_control, nmass_se_trt, 
                nmass_n_control, nmass_n_trt,
                

                pmass_mean_control, pmass_mean_trt, 
                pmass_sd_control, pmass_sd_trt, 
                pmass_se_control, pmass_se_trt, 
                pmass_n_control, pmass_n_trt,
                

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
  mutate(trait = factor(trait, levels = c("nmass", "pmass", "leafnp"))) %>%
  arrange(site, trait, fert) %>%
  write.csv("../data_summaries/Fisher2013_summary.csv", row.names = F)

