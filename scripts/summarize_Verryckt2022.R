# determine Verryckt et al., 2022 treatment mean +/- SE

# Libraries
library(tidyverse)

# Load Verryckt et al. (2022) data and calculate Narea and Parea based on
# SLA and leaf N concentration
verryckt_data <- read.csv("../raw_data/Verryckt_2022_mature_trees.csv") %>%
  slice(-1) %>%
  mutate(across(Vcmax:Height, \(x) as.numeric(x)),
         lma = 1 / (SLA * 0.0001),
         Narea = (N / 100) * lma,
         Parea = (P / 100) * lma,
         leaf_np = as.numeric(N) / as.numeric(P))

# Subset into sites - treating as independent data points for meta analysis
# since sites vary quite a bit spatially and have different climate patterns
verryckt_data_nou <- subset(verryckt_data, Site == "NOU")
verryckt_data_par <- subset(verryckt_data, Site == "PAR")

# First, data from NOU. Selecting only wet season to control for seasonality
# and to ensure that water availability does not limit photosynthetic traits.
# Also selecting post-fertilization measurements (3 years post fertilization)
# and full-sun branches
verryckt_data_nou_summary <- verryckt_data_nou %>%
  filter(Season == "wet" & Year == 2019 & Branch == "T") %>%
  group_by(Fertilisation, Site) %>%
  summarize(
    
    n = length(lma),
    
    SPAD_mean = mean(CCI, na.rm = TRUE),
    SPAD_sd = sd(CCI, na.rm = TRUE),
    SPAD_se = SPAD_sd / sqrt(n),
    
    lma_mean = mean(lma, na.rm = TRUE),
    lma_sd = sd(lma, na.rm =TRUE),
    lma_se = lma_sd / sqrt(n),
    
    Nmass_mean = mean(N, na.rm = TRUE),
    Nmass_sd = sd(N, na.rm = TRUE),
    Nmass_se = Nmass_sd / sqrt(n),
    
    Narea_mean = mean(Narea, na.rm = TRUE),
    Narea_sd = sd(Narea, na.rm = TRUE),
    Narea_se = Narea_sd / sqrt(n),
    
    Pmass_mean = mean(P, na.rm = TRUE),
    Pmass_sd = sd(P, na.rm = TRUE),
    Pmass_se = Pmass_sd / sqrt(n),
    
    Parea_mean = mean(Parea, na.rm = TRUE),
    Parea_sd = sd(Parea, na.rm = TRUE),
    Parea_se = Parea_sd / sqrt(n),
    
    leafnp_mean = mean(leaf_np, na.rm = TRUE),
    leafnp_sd = sd(leaf_np, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(n),
    
    Vcmax_mean = mean(Vcmax, na.rm = TRUE),
    Vcmax_sd = sd(Vcmax, na.rm = TRUE),
    Vcmax_se = Vcmax_sd / sqrt(n),
    
    Jmax_mean = mean(Jmax, na.rm = TRUE),
    Jmax_sd = sd(Jmax, na.rm = TRUE),
    Jmax_se = Jmax_sd / sqrt(n),
    
    Rd_mean = mean(Rd, na.rm = TRUE),
    Rd_sd = sd(Rd, na.rm = TRUE),
    Rd_se = Rd_sd / sqrt(n)
  )
  


# Second, data from PAR. Selecting only wet season to control for seasonality
# and to ensure that water availability does not limit photosynthetic traits 
verryckt_data_par_summary <- verryckt_data_par  %>%
  filter(Season == "wet" & Year == 2019 & Branch == "T") %>%
  group_by(Fertilisation, Site) %>%
  summarize(
    
    n = length(lma),
    
    SPAD_mean = mean(CCI, na.rm = TRUE),
    SPAD_sd = sd(CCI, na.rm = TRUE),
    SPAD_se = SPAD_sd / sqrt(n),
    
    lma_mean = mean(lma, na.rm = TRUE),
    lma_sd = sd(lma, na.rm =TRUE),
    lma_se = lma_sd / sqrt(n),
    
    Nmass_mean = mean(N, na.rm = TRUE),
    Nmass_sd = sd(N, na.rm = TRUE),
    Nmass_se = Nmass_sd / sqrt(n),
    
    Narea_mean = mean(Narea, na.rm = TRUE),
    Narea_sd = sd(Narea, na.rm = TRUE),
    Narea_se = Narea_sd / sqrt(n),
    
    Pmass_mean = mean(P, na.rm = TRUE),
    Pmass_sd = sd(P, na.rm = TRUE),
    Pmass_se = Pmass_sd / sqrt(n),
    
    Parea_mean = mean(Parea, na.rm = TRUE),
    Parea_sd = sd(Parea, na.rm = TRUE),
    Parea_se = Parea_sd / sqrt(n),
    
    leafnp_mean = mean(leaf_np, na.rm = TRUE),
    leafnp_sd = sd(leaf_np, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(n),
    
    Vcmax_mean = mean(Vcmax, na.rm = TRUE),
    Vcmax_sd = sd(Vcmax, na.rm = TRUE),
    Vcmax_se = Vcmax_sd / sqrt(n),
    
    Jmax_mean = mean(Jmax, na.rm = TRUE),
    Jmax_sd = sd(Jmax, na.rm = TRUE),
    Jmax_se = Jmax_sd / sqrt(n),
    
    Rd_mean = mean(Rd, na.rm = TRUE),
    Rd_sd = sd(Rd, na.rm = TRUE),
    Rd_se = Rd_sd / sqrt(n)
    
  )

## Merge datasets and write to .csv
verryckt_data_nou_summary %>%
  full_join(verryckt_data_par_summary) %>%
  mutate(Fertilisation = factor(Fertilisation, levels = c("C", "N", "P", "NP"))) %>%
  arrange(Site, Fertilisation) %>%
  write.csv("../data_summaries/Verryckt_2022_summary.csv", row.names = F)




