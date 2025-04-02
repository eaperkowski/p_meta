# determine Verryckt et al., 2022 treatment mean +/- SE

# Libraries
library(tidyverse)

# Load Verryckt et al. (2022) data and calculate Narea and Parea based on
# SLA and leaf N concentration
verryckt_data <- read.csv("../../raw_data/Verryckt_2022_mature_trees.csv") %>%
  slice(-1) %>%
  mutate(across(Vcmax:Height, \(x) as.numeric(x)),
         lma = 1 / (SLA * 0.0001),
         Narea = (N / 100) * lma,
         Parea = (P / 100) * lma,
         leaf_np = as.numeric(N) / as.numeric(P),
         Fertilisation = factor(Fertilisation, levels = c("C", "N", "P", "NP")),
         spName_full = str_c(tolower(Genus), "_", Species))


# Selecting only wet season to control for seasonality and to ensure that 
# water availability does not limit photosynthetic traits. Also selecting 
# post-fertilization measurements (3 years post fertilization) and full-sun 
# branches
verryckt_data_summary <- verryckt_data %>%
  filter(Season == "wet" & Year == 2019 & Branch == "T") %>%
  group_by(Fertilisation, Site, spName_full) %>%
  summarize(
    SPAD_n = sum(!is.na(CCI)),
    SPAD_mean = mean(CCI, na.rm = TRUE),
    SPAD_sd = sd(CCI, na.rm = TRUE),
    SPAD_se = SPAD_sd / sqrt(SPAD_n),
    
    lma_n = sum(!is.na(lma)),
    lma_mean = mean(lma, na.rm = TRUE),
    lma_sd = sd(lma, na.rm =TRUE),
    lma_se = lma_sd / sqrt(lma_n),
    
    Nmass_n = sum(!is.na(N)),
    Nmass_mean = mean(N, na.rm = TRUE),
    Nmass_sd = sd(N, na.rm = TRUE),
    Nmass_se = Nmass_sd / sqrt(Nmass_n),
    
    Narea_n = sum(!is.na(Narea)),
    Narea_mean = mean(Narea, na.rm = TRUE),
    Narea_sd = sd(Narea, na.rm = TRUE),
    Narea_se = Narea_sd / sqrt(Narea_n),
    
    Pmass_n = sum(!is.na(P)),
    Pmass_mean = mean(P, na.rm = TRUE),
    Pmass_sd = sd(P, na.rm = TRUE),
    Pmass_se = Pmass_sd / sqrt(Pmass_n),
    
    Parea_n = sum(!is.na(Parea)),
    Parea_mean = mean(Parea, na.rm = TRUE),
    Parea_sd = sd(Parea, na.rm = TRUE),
    Parea_se = Parea_sd / sqrt(Parea_n),
    
    leafnp_n = sum(!is.na(leaf_np)),
    leafnp_mean = mean(leaf_np, na.rm = TRUE),
    leafnp_sd = sd(leaf_np, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(leafnp_n),
    
    Vcmax_n = sum(!is.na(Vcmax)),
    Vcmax_mean = mean(Vcmax, na.rm = TRUE),
    Vcmax_sd = sd(Vcmax, na.rm = TRUE),
    Vcmax_se = Vcmax_sd / sqrt(Vcmax_n),
    
    Jmax_n = sum(!is.na(Jmax)),
    Jmax_mean = mean(Jmax, na.rm = TRUE),
    Jmax_sd = sd(Jmax, na.rm = TRUE),
    Jmax_se = Jmax_sd / sqrt(Jmax_n),
    
    TPU_n = sum(!is.na(TPU)),
    TPU_mean = mean(TPU, na.rm = TRUE),
    TPU_sd = sd(TPU, na.rm = TRUE),
    TPU_se = TPU_sd / sqrt(TPU_n),
    
    Rd_n = sum(!is.na(Rd)),
    Rd_mean = mean(Rd, na.rm = TRUE),
    Rd_sd = sd(Rd, na.rm = TRUE),
    Rd_se = Rd_sd / sqrt(Rd_n)
  )
  
# Prep for easy merge into compiled datasheet
verryckt_data_summary_control <- verryckt_data_summary %>%
  ungroup(Fertilisation) %>%
  filter(Fertilisation == "C") %>%
  dplyr::select(-Fertilisation)
names(verryckt_data_summary_control)[3:46] <- str_c(names(verryckt_data_summary_control)[3:46], "_control")

verryckt_data_summary_treatment <- verryckt_data_summary %>%
  filter(Fertilisation != "C")
names(verryckt_data_summary_treatment)[4:47] <- str_c(names(verryckt_data_summary_treatment)[4:47], "_trt")


# Format into easy merge into compiled datasheet, write to .csv
verryckt_data_summary_control %>%
  full_join(verryckt_data_summary_treatment, by = c("Site", "spName_full")) %>%
  dplyr::select(Site, spName_full, Fertilisation, 
                
                SPAD_mean_control, SPAD_mean_trt,
                SPAD_sd_control, SPAD_sd_trt,
                SPAD_se_control, SPAD_se_trt,
                SPAD_n_control, SPAD_n_trt,
                
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
                
                Vcmax_mean_control, Vcmax_mean_trt,
                Vcmax_sd_control, Vcmax_sd_trt,
                Vcmax_se_control, Vcmax_se_trt,
                Vcmax_n_control, Vcmax_n_trt,
                
                Jmax_mean_control, Jmax_mean_trt,
                Jmax_sd_control, Jmax_sd_trt,
                Jmax_se_control, Jmax_se_trt,
                Jmax_n_control, Jmax_n_trt,
                
                TPU_mean_control, TPU_mean_trt,
                TPU_sd_control, TPU_sd_trt,
                TPU_se_control, TPU_se_trt,
                TPU_n_control, TPU_n_trt,
                
                Rd_mean_control, Rd_mean_trt,
                Rd_sd_control, Rd_sd_trt,
                Rd_se_control, Rd_se_trt,
                Rd_n_control, Rd_n_trt) %>%
  pivot_longer(cols = SPAD_mean_control:Rd_n_trt,
               names_to = c("trait", "stat", "trt"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  mutate(trait = factor(trait, levels = c("SPAD", "lma", "Nmass", "Narea",
                                          "Pmass", "Parea", "leafnp", "Vcmax",
                                          "Jmax", "TPU", "Rd"))) %>%
  group_by(Site, spName_full) %>%
  filter(n_control > 1 & n_trt > 1) %>%
  filter(n_distinct(Fertilisation) == 3) %>%
  arrange(Site, trait, spName_full, Fertilisation) %>%
  write.csv("../../data_summaries/species_level/Verryckt2022_summary_by_spp.csv",
            row.names = F)




