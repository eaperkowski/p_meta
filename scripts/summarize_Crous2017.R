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
crous_data_summary <- crous_data %>%
  group_by(Treatm) %>%
  summarize(

    Asat_n = sum(!is.na(PSsat)),
    Asat_mean = mean(PSsat, na.rm = TRUE),
    Asat_sd = sd(PSsat, na.rm = TRUE),
    Asat_se = Asat_sd / sqrt(Asat_n),
    
    Vcmax_n = sum(!is.na(Vcsat)),
    Vcmax_mean = mean(Vcsat, na.rm = TRUE),
    Vcmax_sd = sd(Vcsat, na.rm = TRUE),
    Vcmax_se = Vcmax_sd / sqrt(Vcmax_n),
    
    Jmax_n = sum(!is.na(Jsat)),
    Jmax_mean = mean(Jsat, na.rm = TRUE),
    Jmax_sd = sd(Jsat, na.rm = TRUE),
    Jmax_se = Jmax_sd / sqrt(Jmax_n),
    
    gsw_n = sum(!is.na(Gssat)),
    gsw_mean = mean(Gssat, na.rm = TRUE),
    gsw_sd = sd(Gssat, na.rm = TRUE),
    gsw_se = gsw_sd / sqrt(gsw_n),
    
    rd_n = sum(!is.na(Rdark)),
    rd_mean = mean(Rdark, na.rm = TRUE),
    rd_sd = sd(Rdark, na.rm = TRUE),
    rd_se = rd_sd / sqrt(rd_n),
    
    lma_n = sum(!is.na(LMA)),
    lma_mean = mean(LMA, na.rm = TRUE),
    lma_sd = sd(LMA, na.rm =TRUE),
    lma_se = lma_sd / sqrt(lma_n),
    
    Nmass_n = sum(!is.na(Nmass)),
    Nmass_mean = mean(Nmass, na.rm = TRUE),
    Nmass_sd = sd(Nmass, na.rm = TRUE),
    Nmass_se = Nmass_sd / sqrt(Nmass_n),
    
    Narea_n = sum(!is.na(Narea)),
    Narea_mean = mean(Narea, na.rm = TRUE),
    Narea_sd = sd(Narea, na.rm = TRUE),
    Narea_se = Narea_sd / sqrt(Narea_n),
    
    Pmass_n = sum(!is.na(Pmass)),
    Pmass_mean = mean(Pmass, na.rm = TRUE),
    Pmass_sd = sd(Pmass, na.rm = TRUE),
    Pmass_se = Pmass_sd / sqrt(Pmass_n),
    
    Parea_n = sum(!is.na(Parea)),
    Parea_mean = mean(Parea, na.rm = TRUE),
    Parea_sd = sd(Parea, na.rm = TRUE),
    Parea_se = Parea_sd / sqrt(Parea_n),
    
    leafnp_n = sum(!is.na(N.P)),
    leafnp_mean = mean(N.P, na.rm = TRUE),
    leafnp_sd = sd(N.P, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(leafnp_n),
    
    pnue_n = sum(!is.na(pnue)),
    pnue_mean = mean(pnue, na.rm = TRUE),
    pnue_sd = sd(pnue, na.rm =TRUE),
    pnue_se = pnue_sd / sqrt(pnue_n),
    
    ppue_n = sum(!is.na(ppue)),
    ppue_mean = mean(ppue, na.rm = TRUE),
    ppue_sd = sd(ppue, na.rm =TRUE),
    ppue_se = ppue_sd / sqrt(ppue_n)
    
  ) 


# Prep for easy merge into compiled datasheet
crous_data_summary_control <- crous_data_summary %>%
  filter(Treatm == "LNLP") %>%
  mutate(exp = "crous2017") %>%
  select(-Treatm)
names(crous_data_summary_control)[1:52] <- str_c(names(crous_data_summary_control)[1:52], "_control")

crous_data_summary_treatment <- crous_data_summary %>%
  filter(Treatm != "LNLP") %>%
  mutate(exp = "crous2017")
names(crous_data_summary_treatment)[2:53] <- str_c(names(crous_data_summary_treatment)[2:53], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
crous_data_summary_control %>%
  full_join(crous_data_summary_treatment, by = "exp") %>%
  dplyr::select(exp, Treatm, 
                
                Asat_mean_control, Asat_mean_trt, 
                Asat_sd_control, Asat_sd_trt, 
                Asat_se_control, Asat_se_trt, 
                Asat_n_control, Asat_n_trt,
                
                Vcmax_mean_control, Vcmax_mean_trt, 
                Vcmax_sd_control, Vcmax_sd_trt, 
                Vcmax_se_control, Vcmax_se_trt, 
                Vcmax_n_control, Vcmax_n_trt,
                
                Jmax_mean_control, Jmax_mean_trt, 
                Jmax_sd_control, Jmax_sd_trt, 
                Jmax_se_control, Jmax_se_trt, 
                Jmax_n_control, Jmax_n_trt,
                
                gsw_mean_control, gsw_mean_trt, 
                gsw_sd_control, gsw_sd_trt, 
                gsw_se_control, gsw_se_trt, 
                gsw_n_control, gsw_n_trt,
                
                rd_mean_control, rd_mean_trt, 
                rd_sd_control, rd_sd_trt, 
                rd_se_control, rd_se_trt, 
                rd_n_control, rd_n_trt,
                
                lma_mean_control, lma_mean_trt, 
                lma_sd_control, lma_sd_trt, 
                lma_se_control, lma_se_trt, 
                lma_n_control, lma_n_trt,
                
                Nmass_mean_control, Nmass_mean_trt, 
                Nmass_sd_control, Nmass_sd_trt, 
                Nmass_se_control, Nmass_se_trt, 
                Nmass_n_control, Nmass_n_trt,
                
                Narea_mean_control, Narea_mean_trt, 
                Narea_sd_control, Narea_sd_trt, 
                Narea_se_control, Narea_se_trt, 
                Narea_n_control, Narea_n_trt,
                
                Pmass_mean_control, Pmass_mean_trt, 
                Pmass_sd_control, Pmass_sd_trt, 
                Pmass_se_control, Pmass_se_trt, 
                Pmass_n_control, Pmass_n_trt,
                
                Parea_mean_control, Parea_mean_trt, 
                Parea_sd_control, Parea_sd_trt, 
                Parea_se_control, Parea_se_trt, 
                Parea_n_control, Parea_n_trt,
                
                leafnp_mean_control, leafnp_mean_trt, 
                leafnp_sd_control, leafnp_sd_trt, 
                leafnp_se_control, leafnp_se_trt, 
                leafnp_n_control, leafnp_n_trt,
                
                pnue_mean_control, pnue_mean_trt, 
                pnue_sd_control, pnue_sd_trt, 
                pnue_se_control, pnue_se_trt, 
                pnue_n_control, pnue_n_trt,
                
                ppue_mean_control, ppue_mean_trt, 
                ppue_sd_control, ppue_sd_trt, 
                ppue_se_control, ppue_se_trt, 
                ppue_n_control, ppue_n_trt) %>%
  pivot_longer(cols = Asat_mean_control:ppue_n_trt,
               names_to = c("trait", "stat", "trt"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  arrange(trait) %>%
  write.csv("../data_summaries/Cleland2019_summary.csv", row.names = F)
