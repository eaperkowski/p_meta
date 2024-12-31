# analyze Mo et al., 2019 data

# Libraries
library(tidyverse)

# Load Mo et al. (2019) data
mo_data <- read.csv("../raw_data/Mo_2019_FE.csv")

# Group by and summarize Mo et al. (2019) dataset for easy input into 
# compiled datasheet
head(mo_data)

mo_data %>%
  
  mutate(Narea = (N/1000) * LMA,
         Parea = (P/1000) * LMA) %>%
  group_by(Treatment) %>%
  summarize(n = length(PA),
            anet_area_mean = mean(PA, na.rm = TRUE),
            anet_area_sd = sd(PA, na.rm = TRUE),
            anet_area_se = anet_area_sd / sqrt(n),
            
            anet_mass_mean = mean(PM, na.rm = TRUE),
            anet_mass_sd = sd(PM, na.rm = TRUE),
            anet_mass_se = anet_mass_sd / sqrt(n),
            
            LMA_mean = mean(LMA, na.rm = TRUE),
            LMA_sd = sd(LMA, na.rm =TRUE),
            LMA_se = LMA_sd / sqrt(n),
            
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
            
            PNUE_mean = mean(PNUE, na.rm = TRUE),
            PNUE_sd = sd(PNUE, na.rm = TRUE),
            PNUE_se = PNUE_sd / sqrt(n),
            
            PPUE_mean = mean(PPUE, na.rm = TRUE),
            PPUE_sd = sd(PPUE, na.rm = TRUE),
            PPUE_se = PPUE_sd / sqrt(n),
            
            SP_mean = mean(SP, na.rm = TRUE),
            SP_sd = sd(SP, na.rm = TRUE),
            SP_se = SP_sd / sqrt(n),

            MP_mean = mean(MP, na.rm = TRUE),
            MP_sd = sd(MP, na.rm = TRUE),
            MP_se = MP_sd / sqrt(n),
            
            NP1_mean = mean(NP1, na.rm = TRUE),
            NP1_sd = sd(NP1, na.rm = TRUE),
            NP1_se = NP1_sd / sqrt(n),
            
            RP_mean = mean(RP, na.rm = TRUE),
            RP_sd = sd(RP, na.rm = TRUE),
            RP_se = RP_sd / sqrt(n)
            ) %>%
  write.csv("../raw_data/Mo_2019_FE_summarized.csv", row.names = F)
