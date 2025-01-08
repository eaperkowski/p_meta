# determine Nutnet treatment mean +/- SE

# Libraries
library(tidyverse)

# Load NutNet data from Alissar paper and Firn et al., 2019, filter to
# include only control, P, N, and NP data
nutnet_data <- read.csv("../raw_data/Firn_2019_nature.csv") %>%
  filter(trt %in% c("P", "Control", "N", "NP") & Exclose == 0 & year_trt > 1) %>%
  mutate(trt = factor(trt, levels = c("Control", "N", "P", "NP")))

head(nutnet_data)

unique(nutnet_data$site_code)

# Determine LMA, Narea, Parea based on SLA
nutnet_data <- nutnet_data %>%
  mutate(lma = (1/SLA_v2) * 1000000,
         narea = (leaf_pct_N/100) * lma,
         parea = (leaf_pct_P/100) * lma,
         leaf_np = leaf_pct_N / leaf_pct_P)

# Calculate site-level treatment means +/- SE
nutnet_data %>%
  group_by(site_code, year, trt) %>%
  summarize(
    
    n = length(site_code),
    exp_duration = mean(year_trt),
    
    lma_mean = mean(lma, na.rm = TRUE),
    lma_sd = sd(lma, na.rm = TRUE),
    lma_se = lma_sd / sqrt(n),
    
    nmass_mean = mean(leaf_pct_N, na.rm = TRUE),
    nmass_sd = sd(leaf_pct_N, na.rm = TRUE),
    nmass_se = nmass_sd / sqrt(n),
    
    narea_mean = mean(narea, na.rm = TRUE),
    narea_sd = sd(narea, na.rm = TRUE),
    narea_se = narea_sd / sqrt(n),
    
    pmass_mean = mean(leaf_pct_P, na.rm = TRUE),
    pmass_sd = sd(leaf_pct_P, na.rm = TRUE),
    pmass_se = pmass_sd / sqrt(n),
    
    parea_mean = mean(parea, na.rm = TRUE),
    parea_sd = sd(parea, na.rm = TRUE),
    parea_se = parea_sd / sqrt(n),
    
    leafnp_mean = mean(leaf_np, na.rm = TRUE),
    leafnp_sd = sd(leaf_np, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(n)
    
  ) %>%
  write.csv("../data_summaries/Firn2019_nature_summary.csv", row.names = F)


# Climate data summaries
nutnet_data %>%
  group_by(site_code) %>%
  summarize(latitude = mean(latitude), 
            longitude = mean(longitude), 
            MAT = mean(MAT), 
            MAP = mean(MAP)) %>%
  write.csv("../data_summaries/Firn2019_nature_climate_summary.csv", row.names = F)

