# analyze Mo et al., 2019 data

# Libraries
library(tidyverse)

# Load Mo et al. (2019) data
mo_data <- read.csv("../raw_data/Mo_2019_FE.csv") %>%
  mutate(Narea = (N / 1000) * LMA,
         Parea = (P / 1000) * LMA,
         leaf_np = N / P,
         Treatment = factor(Treatment, levels = c("CT", "N", "P", "NP")))
head(mo_data)

# Calculate experiment summary statistics
mo_data_summary <-  mo_data %>%
  group_by(Treatment) %>%
  summarize(anet_n = sum(!is.na(PA)),
            anet_mean = mean(PA, na.rm = TRUE),
            anet_sd = sd(PA, na.rm = TRUE),
            anet_se = anet_sd / sqrt(anet_n),
            
            lma_n = sum(!is.na(LMA)),
            lma_mean = mean(LMA, na.rm = TRUE),
            lma_sd = sd(LMA, na.rm =TRUE),
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
            
            pnue_n = sum(!is.na(PNUE)),
            pnue_mean = mean(PNUE, na.rm = TRUE),
            pnue_sd = sd(PNUE, na.rm = TRUE),
            pnue_se = pnue_sd / sqrt(pnue_n),
            
            ppue_n = sum(!is.na(PPUE)),
            ppue_mean = mean(PPUE, na.rm = TRUE),
            ppue_sd = sd(PPUE, na.rm = TRUE),
            ppue_se = ppue_sd / sqrt(ppue_n),
            
            structuralp_n = sum(!is.na(SP)),
            structuralp_mean = mean(SP, na.rm = TRUE),
            structuralp_sd = sd(SP, na.rm = TRUE),
            structuralp_se = structuralp_sd / sqrt(structuralp_n),

            metabolicp_n = sum(!is.na(MP)),
            metabolicp_mean = mean(MP, na.rm = TRUE),
            metabolicp_sd = sd(MP, na.rm = TRUE),
            metabolicp_se = metabolicp_sd / sqrt(metabolicp_n),
            
            nucleicp_n = sum(!is.na(NP1)),
            nucleicp_mean = mean(NP1, na.rm = TRUE),
            nucleicp_sd = sd(NP1, na.rm = TRUE),
            nucleicp_se = nucleicp_sd / sqrt(nucleicp_n),
            
            residualp_n = sum(!is.na(RP)),
            residualp_mean = mean(RP, na.rm = TRUE),
            residualp_sd = sd(RP, na.rm = TRUE),
            residualp_se = residualp_sd / sqrt(residualp_n)
            )

# Prep for easy merge into compiled datasheet
mo_data_summary_control <- mo_data_summary %>%
  filter(Treatment == "CT") %>%
  mutate(exp = "Mo2019") %>%
  select(-Treatment)
names(mo_data_summary_control)[1:52] <- str_c(names(mo_data_summary_control)[1:52], "_control")

mo_data_summary_treatment <- mo_data_summary %>%
  filter(Treatment != "CT") %>%
  mutate(exp = "Mo2019")
names(mo_data_summary_treatment)[2:53] <- str_c(names(mo_data_summary_treatment)[2:53], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
mo_data_summary_control %>%
  full_join(mo_data_summary_treatment, by = "exp") %>%
  dplyr::select(Treatment, 
                
                anet_mean_control, anet_mean_trt, 
                anet_sd_control, anet_sd_trt, 
                anet_se_control, anet_se_trt, 
                anet_n_control, anet_n_trt,
                
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
                ppue_n_control, ppue_n_trt,
            
                structuralp_mean_control, structuralp_mean_trt, 
                structuralp_sd_control, structuralp_sd_trt, 
                structuralp_se_control, structuralp_se_trt, 
                structuralp_n_control, structuralp_n_trt,
                
                metabolicp_mean_control, metabolicp_mean_trt, 
                metabolicp_sd_control, metabolicp_sd_trt, 
                metabolicp_se_control, metabolicp_se_trt, 
                metabolicp_n_control, metabolicp_n_trt,
                
                nucleicp_mean_control, nucleicp_mean_trt, 
                nucleicp_sd_control, nucleicp_sd_trt, 
                nucleicp_se_control, nucleicp_se_trt, 
                nucleicp_n_control, nucleicp_n_trt,
                
                residualp_mean_control, residualp_mean_trt, 
                residualp_sd_control, residualp_sd_trt, 
                residualp_se_control, residualp_se_trt, 
                residualp_n_control, residualp_n_trt) %>%
  pivot_longer(cols = anet_mean_control:residualp_n_trt,
               names_to = c("trait", "stat", "trt"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  mutate(trait = factor(trait, levels = c("anet", "lma", "Nmass", "Narea",
                                          "Pmass", "Parea", "leafnp", 
                                          "pnue", "ppue", "structuralp",
                                          "metabolicp", "nucleicp", 
                                          "residualp"))) %>%
  arrange(trait) %>%
  write.csv("../data_summaries/Mo2019_FE_summarized.csv", row.names = F)
