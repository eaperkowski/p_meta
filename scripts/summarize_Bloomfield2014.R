# determine Bloomfield et al. (2014) treatment mean +/- SE

# Libraries
library(tidyverse)

# Load Bloomfield et al. (2014) data and change varname for treatment
# such that Control (in original datasheet) is the +P treatment and
# minus P is the Control treatment
bloomfield_data <- read.csv("../raw_data/Bloomfield_2014_FPB.csv") %>%
  mutate(Treatment = ifelse(Treatment == "Control", "P", 
                            ifelse(Treatment == "Minus P", "Control", NA)),
         Treatment = factor(Treatment, levels = c("Control", "P")))

head(bloomfield_data)

# How many unique species and growth forms?
unique(bloomfield_data$Species)
unique(bloomfield_data$Adult.stature)

# Group by treatment and timepoint, then summarize Bloomfield et al. (2014)
# for easy input into compiled datasheet
bloomfield_data %>%
  group_by(Treatment, Months.in) %>%
  summarize(
    
    n = length(Months.in),
    
    Amax_mean = mean(Amax.a, na.rm = TRUE),
    Amax_sd = sd(Amax.a, na.rm = TRUE),
    Amax_se = Amax_sd / sqrt(n),
    
    Asat_mean = mean(Asat.a, na.rm = TRUE),
    Asat_sd = sd(Asat.a, na.rm = TRUE),
    Asat_se = Asat_sd / sqrt(n),
    
    Vcmax_mean = mean(Vcmax.a, na.rm = TRUE),
    Vcmax_sd = sd(Vcmax.a, na.rm = TRUE),
    Vcmax_se = Vcmax_sd / sqrt(n),
    
    Jmax_mean = mean(Jmax.a, na.rm = TRUE),
    Jmax_sd = sd(Jmax.a, na.rm = TRUE),
    Jmax_se = Jmax_sd / sqrt(n),
    
    TPU_mean = mean(TPU.a, na.rm = TRUE),
    TPU_sd = sd(TPU.a, na.rm = TRUE),
    TPU_se = TPU_sd / sqrt(n),
    
    Rd_mean = mean(Rdark.a, na.rm = TRUE),
    Rd_sd = sd(Rdark.a, na.rm = TRUE),
    Rd_se = Rd_sd / sqrt(n),
    
    cica_mean = mean(Ci.Ca, na.rm = TRUE),
    cica_sd = sd(Ci.Ca, na.rm = TRUE),
    cica_se = cica_sd / sqrt(n),
    
    lma_mean = mean(LMA, na.rm = TRUE),
    lma_sd = sd(LMA, na.rm = TRUE),
    lma_se = lma_sd / sqrt(n),
    
    nmass_mean = mean(LeafN, na.rm = TRUE),
    nmass_sd = sd(LeafN, na.rm = TRUE),
    nmass_se = nmass_sd / sqrt(n),
    
    narea_mean = mean(N.area, na.rm = TRUE),
    narea_sd = sd(N.area, na.rm = TRUE),
    narea_se = narea_sd / sqrt(n),

    pmass_mean = mean(LeafP, na.rm = TRUE),
    pmass_sd = sd(LeafP, na.rm = TRUE),
    pmass_se = pmass_sd / sqrt(n),
    
    parea_mean = mean(P.area, na.rm = TRUE),
    parea_sd = sd(P.area, na.rm = TRUE),
    parea_se = parea_sd / sqrt(n),
    
    leafnp_mean = mean(RatioNP, na.rm = TRUE),
    leafnp_sd = sd(RatioNP, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(n),
    
    pnue_mean = mean(PNUE, na.rm = TRUE),
    pnue_sd = sd(PNUE, na.rm = TRUE),
    pnue_se = pnue_sd / sqrt(n),
    
    ppue_mean = mean(PPUE, na.rm = TRUE),
    ppue_sd = sd(PPUE, na.rm = TRUE),
    ppue_se = ppue_sd / sqrt(n)
    
  ) %>%
  write.csv("../data_summaries/Bloomfield2014_FPB_summarized.csv", row.names = F)
  