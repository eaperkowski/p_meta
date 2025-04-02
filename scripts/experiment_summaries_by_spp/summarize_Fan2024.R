# analyze Fan et al. (2024) data

# Libraries
library(tidyverse)

# Load Mo et al. (2019) data
fan_data <- read.csv("../../raw_data/Fan_2024_pFraction_data.csv") %>%
  slice(1:122) %>%
  mutate(leaf_np = nmass / pmass)
head(fan_data)


# Summarize Fan et al. (2024) dataset 
fan_data_summary <- fan_data %>%

  group_by(treatment) %>%
  summarize(lma_n = sum(!is.na(lma)),
            lma_mean = mean(lma, na.rm = TRUE),
            lma_sd = sd(lma, na.rm =TRUE),
            lma_se = lma_sd / sqrt(lma_n),
            
            spad_n = sum(!is.na(spad)),
            spad_mean = mean(spad, na.rm = TRUE),
            spad_sd = sd(spad, na.rm =TRUE),
            spad_se = spad_sd / sqrt(spad_n),
            
            nmass_n = sum(!is.na(nmass)),
            nmass_mean = mean(nmass, na.rm = TRUE),
            nmass_sd = sd(nmass, na.rm =TRUE),
            nmass_se = nmass_sd / sqrt(nmass_n),
            
            narea_n = sum(!is.na(narea)),
            narea_mean = mean(narea, na.rm = TRUE),
            narea_sd = sd(narea, na.rm =TRUE),
            narea_se = narea_sd / sqrt(narea_n),
            
            pmass_n = sum(!is.na(parea)),
            pmass_mean = mean(pmass, na.rm = TRUE),
            pmass_sd = sd(pmass, na.rm =TRUE),
            pmass_se = pmass_sd / sqrt(pmass_n),
            
            parea_n = sum(!is.na(parea)),
            parea_mean = mean(parea, na.rm = TRUE),
            parea_sd = sd(parea, na.rm =TRUE),
            parea_se = parea_sd / sqrt(parea_n),
            
            leafnp_n = sum(!is.na(leaf_np)),
            leafnp_mean = mean(leaf_np, na.rm = TRUE),
            leafnp_sd = sd(leaf_np, na.rm =TRUE),
            leafnp_se = leafnp_sd / sqrt(leafnp_n),
            
            rgr_n = sum(!is.na(rgr)),
            rgr_mean = mean(rgr, na.rm = TRUE),
            rgr_sd = sd(rgr, na.rm =TRUE),
            rgr_se = rgr_sd / sqrt(rgr_n),
            
            amax_n = sum(!is.na(Amax)),
            amax_mean = mean(Amax, na.rm = TRUE),
            amax_sd = sd(Amax, na.rm =TRUE),
            amax_se = amax_sd / sqrt(amax_n),
            
            vcmax_n = sum(!is.na(Vcmax)),
            vcmax_mean = mean(Vcmax, na.rm = TRUE),
            vcmax_sd = sd(Vcmax, na.rm =TRUE),
            vcmax_se = vcmax_sd / sqrt(vcmax_n),
            
            jmax_n = sum(!is.na(Jmax)),
            jmax_mean = mean(Jmax, na.rm = TRUE),
            jmax_sd = sd(Jmax, na.rm =TRUE),
            jmax_se = jmax_sd / sqrt(jmax_n),
            
            pnue_n = sum(!is.na(PNUE)),
            pnue_mean = mean(PNUE, na.rm = TRUE),
            pnue_sd = sd(PNUE, na.rm =TRUE),
            pnue_se = pnue_sd / sqrt(pnue_n),
            
            ppue_n = sum(!is.na(PPUE)),
            ppue_mean = mean(PPUE, na.rm = TRUE),
            ppue_sd = sd(PPUE, na.rm =TRUE),
            ppue_se = ppue_sd / sqrt(ppue_n),
            
            pi_n = sum(!is.na(Pi)),
            pi_mean = mean(Pi, na.rm = TRUE),
            pi_sd = sd(Pi, na.rm =TRUE),
            pi_se = pi_sd / sqrt(pi_n),
            
            esterp_n = sum(!is.na(EsterP)),
            esterp_mean = mean(EsterP, na.rm = TRUE),
            esterp_sd = sd(EsterP, na.rm =TRUE),
            esterp_se = esterp_sd / sqrt(esterp_n),
            
            nucleicp_n = sum(!is.na(Nucleic.P)),
            nucleicp_mean = mean(Nucleic.P, na.rm = TRUE),
            nucleicp_sd = sd(Nucleic.P, na.rm =TRUE),
            nucleicp_se = nucleicp_sd / sqrt(nucleicp_n),
            
            lipidp_n = sum(!is.na(Lipid.P)),
            lipidp_mean = mean(Lipid.P, na.rm = TRUE),
            lipidp_sd = sd(Lipid.P, na.rm =TRUE),
            lipidp_se = lipidp_sd / sqrt(lipidp_n),
            
            residualp_n = sum(!is.na(Residual.P)),
            residualp_mean = mean(Residual.P, na.rm = TRUE),
            residualp_sd = sd(Residual.P, na.rm =TRUE),
            residualp_se = residualp_sd / sqrt(residualp_n))


# Prep for easy merge into compiled datasheet
fan_data_summary_control <- fan_data_summary %>%
  filter(treatment == "CK") %>%
  mutate(exp = "Fan2024") %>%
  select(-treatment)
names(fan_data_summary_control)[1:72] <- str_c(names(fan_data_summary_control)[1:72], "_control")

fan_data_summary_treatment <- fan_data_summary %>%
  filter(treatment != "CK") %>%
  mutate(exp = "Fan2024")
names(fan_data_summary_treatment)[2:73] <- str_c(names(fan_data_summary_treatment)[2:73], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
fan_data_summary_control %>%
  full_join(fan_data_summary_treatment, by = "exp") %>%
  dplyr::select(exp, treatment,
                
                lma_mean_control, lma_mean_trt, 
                lma_sd_control, lma_sd_trt, 
                lma_se_control, lma_se_trt, 
                lma_n_control, lma_n_trt,
                
                nmass_mean_control, nmass_mean_trt, 
                nmass_sd_control, nmass_sd_trt, 
                nmass_se_control, nmass_se_trt, 
                nmass_n_control, nmass_n_trt,
                
                narea_mean_control, narea_mean_trt, 
                narea_sd_control, narea_sd_trt, 
                narea_se_control, narea_se_trt, 
                narea_n_control, narea_n_trt,
                
                pmass_mean_control, pmass_mean_trt, 
                pmass_sd_control, pmass_sd_trt, 
                pmass_se_control, pmass_se_trt, 
                pmass_n_control, pmass_n_trt,
                
                parea_mean_control, parea_mean_trt, 
                parea_sd_control, parea_sd_trt, 
                parea_se_control, parea_se_trt, 
                parea_n_control, parea_n_trt,
                
                leafnp_mean_control, leafnp_mean_trt, 
                leafnp_sd_control, leafnp_sd_trt, 
                leafnp_se_control, leafnp_se_trt, 
                leafnp_n_control, leafnp_n_trt,
                
                rgr_mean_control, rgr_mean_trt, 
                rgr_sd_control, rgr_sd_trt, 
                rgr_se_control, rgr_se_trt, 
                rgr_n_control, rgr_n_trt,
                
                amax_mean_control, amax_mean_trt, 
                amax_sd_control, amax_sd_trt, 
                amax_se_control, amax_se_trt, 
                amax_n_control, amax_n_trt,
                
                vcmax_mean_control, vcmax_mean_trt, 
                vcmax_sd_control, vcmax_sd_trt, 
                vcmax_se_control, vcmax_se_trt, 
                vcmax_n_control, vcmax_n_trt,
                
                jmax_mean_control, jmax_mean_trt, 
                jmax_sd_control, jmax_sd_trt, 
                jmax_se_control, jmax_se_trt, 
                jmax_n_control, jmax_n_trt,
                
                pnue_mean_control, pnue_mean_trt, 
                pnue_sd_control, pnue_sd_trt, 
                pnue_se_control, pnue_se_trt, 
                pnue_n_control, pnue_n_trt,
                
                ppue_mean_control, ppue_mean_trt, 
                ppue_sd_control, ppue_sd_trt, 
                ppue_se_control, ppue_se_trt, 
                ppue_n_control, ppue_n_trt,
                
                pi_mean_control, pi_mean_trt, 
                pi_sd_control, pi_sd_trt, 
                pi_se_control, pi_se_trt, 
                pi_n_control, pi_n_trt,
                
                esterp_mean_control, esterp_mean_trt, 
                esterp_sd_control, esterp_sd_trt, 
                esterp_se_control, esterp_se_trt, 
                esterp_n_control, esterp_n_trt,
                
                nucleicp_mean_control, nucleicp_mean_trt, 
                nucleicp_sd_control, nucleicp_sd_trt, 
                nucleicp_se_control, nucleicp_se_trt, 
                nucleicp_n_control, nucleicp_n_trt,
                
                lipidp_mean_control, lipidp_mean_trt, 
                lipidp_sd_control, lipidp_sd_trt, 
                lipidp_se_control, lipidp_se_trt, 
                lipidp_n_control, lipidp_n_trt,
                
                residualp_mean_control, residualp_mean_trt, 
                residualp_sd_control, residualp_sd_trt, 
                residualp_se_control, residualp_se_trt, 
                residualp_n_control, residualp_n_trt) %>%
  pivot_longer(cols = lma_mean_control:residualp_n_trt,
               names_to = c("trait", "stat", "trt"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  write.csv("../data_summaries/Fan2024_pFraction_data_summary.csv", row.names = F)
