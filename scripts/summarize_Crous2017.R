# determine Crous et al. (2017) treatment mean +/- SE

# Libraries
library(tidyverse)

# Load Crous et al. (2017) data
crous_data <- read.csv("../raw_data/Crous_2017_NewPhyt.csv") %>%
  mutate(Treatm = 
           factor(Treatm, 
                  levels = c("LNLP", "HNLP", "LNHP", "HNHP")),
         pnue = PSsat / Narea,
         ppue = PSsat / Parea) 
head(crous_data)

# 
crous_data %>%
  group_by(Treatm) %>%
  summarize(
    
    n = length(Treatm),
    
    Asat_mean = mean(PSsat, na.rm = TRUE),
    Asat_sd = sd(PSsat, na.rm = TRUE),
    Asat_se = Asat_sd / sqrt(n),
    
    Vcmax_mean = mean(Vcsat, na.rm = TRUE),
    Vcmax_sd = sd(Vcsat, na.rm = TRUE),
    Vcmax_se = Vcmax_sd / sqrt(n),
    
    Jmax_mean = mean(Jsat, na.rm = TRUE),
    Jmax_sd = sd(Jsat, na.rm = TRUE),
    Jmax_se = Jmax_sd / sqrt(n),
    
    gsw_mean = mean(Gssat, na.rm = TRUE),
    gsw_sd = sd(Gssat, na.rm = TRUE),
    gsw_se = gsw_sd / sqrt(n),
    
    rd_mean = mean(Rdark, na.rm = TRUE),
    rd_sd = sd(Rdark, na.rm = TRUE),
    rd_se = rd_sd / sqrt(n),
    
    lma_mean = mean(LMA, na.rm = TRUE),
    lma_sd = sd(LMA, na.rm =TRUE),
    lma_se = lma_sd / sqrt(n),
    
    Nmass_mean = mean(Nmass, na.rm = TRUE),
    Nmass_sd = sd(Nmass, na.rm = TRUE),
    Nmass_se = Nmass_sd / sqrt(n),
    
    Narea_mean = mean(Narea, na.rm = TRUE),
    Narea_sd = sd(Narea, na.rm = TRUE),
    Narea_se = Narea_sd / sqrt(n),
    
    Pmass_mean = mean(Pmass, na.rm = TRUE),
    Pmass_sd = sd(Pmass, na.rm = TRUE),
    Pmass_se = Pmass_sd / sqrt(n),
    
    Parea_mean = mean(Parea, na.rm = TRUE),
    Parea_sd = sd(Parea, na.rm = TRUE),
    Parea_se = Parea_sd / sqrt(n),
    
    leafnp_mean = mean(N.P, na.rm = TRUE),
    leafnp_sd = sd(N.P, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(n),
    
    pnue_mean = mean(pnue, na.rm = TRUE),
    pnue_sd = sd(pnue, na.rm =TRUE),
    pnue_se = pnue_sd / sqrt(n),
    
    ppue_mean = mean(ppue, na.rm = TRUE),
    ppue_sd = sd(ppue, na.rm =TRUE),
    ppue_se = ppue_sd / sqrt(n)
    
  ) %>%
  write.csv("../data_summaries/Crous_2017_NewPhyt_summarized.csv", row.names = F)
