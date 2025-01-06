# analyze Fan et al. (2024) data



# Libraries
library(tidyverse)

# Load Mo et al. (2019) data
fan_data <- read.csv("../raw_data/Fan_2024_pFraction_data.csv")

# Group by and summarize Fan et al. (2024) dataset for easy input into 
# compiled datasheet
head(fan_data)

fan_data %>%

  group_by(treatment) %>%
  summarize(n = length(lma),
            
            lma_mean = mean(lma, na.rm = TRUE),
            lma_sd = sd(lma, na.rm =TRUE),
            lma_se = lma_sd / sqrt(n),
            
            spad_mean = mean(spad, na.rm = TRUE),
            spad_sd = sd(spad, na.rm =TRUE),
            spad_se = spad_sd / sqrt(n),
            
            nmass_mean = mean(nmass, na.rm = TRUE),
            nmass_sd = sd(nmass, na.rm =TRUE),
            nmass_se = nmass_sd / sqrt(n),
            
            narea_mean = mean(narea, na.rm = TRUE),
            narea_sd = sd(narea, na.rm =TRUE),
            narea_se = narea_sd / sqrt(n),
            
            pmass_mean = mean(pmass, na.rm = TRUE),
            pmass_sd = sd(pmass, na.rm =TRUE),
            pmass_se = pmass_sd / sqrt(n),
            
            parea_mean = mean(parea, na.rm = TRUE),
            parea_sd = sd(parea, na.rm =TRUE),
            parea_se = parea_sd / sqrt(n),
            
            rgr_mean = mean(rgr, na.rm = TRUE),
            rgr_sd = sd(rgr, na.rm =TRUE),
            rgr_se = rgr_sd / sqrt(n),
            
            amax_mean = mean(Amax, na.rm = TRUE),
            amax_sd = sd(Amax, na.rm =TRUE),
            amax_se = amax_sd / sqrt(n),
            
            vcmax_mean = mean(Vcmax, na.rm = TRUE),
            vcmax_sd = sd(Vcmax, na.rm =TRUE),
            vcmax_se = vcmax_sd / sqrt(n),
            
            jmax_mean = mean(Jmax, na.rm = TRUE),
            jmax_sd = sd(Jmax, na.rm =TRUE),
            jmax_se = jmax_sd / sqrt(n),
            
            ppue_mean = mean(PPUE, na.rm = TRUE),
            ppue_sd = sd(PPUE, na.rm =TRUE),
            ppue_se = ppue_sd / sqrt(n),
            
            pnue_mean = mean(PNUE, na.rm = TRUE),
            pnue_sd = sd(PNUE, na.rm =TRUE),
            pnue_se = pnue_sd / sqrt(n),
            
            pi_mean = mean(Pi, na.rm = TRUE),
            pi_sd = sd(Pi, na.rm =TRUE),
            pi_se = pi_sd / sqrt(n),
            
            esterp_mean = mean(EsterP, na.rm = TRUE),
            esterp_sd = sd(EsterP, na.rm =TRUE),
            esterp_se = esterp_sd / sqrt(n),
            
            nucleicp_mean = mean(Nucleic.P, na.rm = TRUE),
            nucleicp_sd = sd(Nucleic.P, na.rm =TRUE),
            nucleicp_se = nucleicp_sd / sqrt(n),
            
            lipidp_mean = mean(Lipid.P, na.rm = TRUE),
            lipidp_sd = sd(Lipid.P, na.rm =TRUE),
            lipidp_se = lipidp_sd / sqrt(n),
            
            residualp_mean = mean(Residual.P, na.rm = TRUE),
            residualp_sd = sd(Residual.P, na.rm =TRUE),
            residualp_se = residualp_sd / sqrt(n)
            
            ) %>%
  slice(-1) %>%
  write.csv("../data_summaries/Fan2024_pFraction_data_summary.csv", row.names = F)
