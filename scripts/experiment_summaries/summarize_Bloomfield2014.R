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
unique(bloomfield_data$Harvest)

# Group by treatment and timepoint, then summarize Bloomfield et al. (2014)
# for easy input into compiled datasheet
bloomfield_data_summary <- bloomfield_data %>%
  group_by(Treatment, Months.in) %>%
  summarize(
    
    Amax_n = sum(!is.na(Amax.a)),
    Amax_mean = mean(Amax.a, na.rm = TRUE),
    Amax_sd = sd(Amax.a, na.rm = TRUE),
    Amax_se = Amax_sd / sqrt(Amax_n),
    
    Asat_n = sum(!is.na(Asat.a)),
    Asat_mean = mean(Asat.a, na.rm = TRUE),
    Asat_sd = sd(Asat.a, na.rm = TRUE),
    Asat_se = Asat_sd / sqrt(Asat_n),
    
    Vcmax_n = sum(!is.na(Vcmax.a)),
    Vcmax_mean = mean(Vcmax.a, na.rm = TRUE),
    Vcmax_sd = sd(Vcmax.a, na.rm = TRUE),
    Vcmax_se = Vcmax_sd / sqrt(Vcmax_n),
    
    Jmax_n = sum(!is.na(Jmax.a)),
    Jmax_mean = mean(Jmax.a, na.rm = TRUE),
    Jmax_sd = sd(Jmax.a, na.rm = TRUE),
    Jmax_se = Jmax_sd / sqrt(Jmax_n),
    
    TPU_n = sum(!is.na(TPU.a)),
    TPU_mean = mean(TPU.a, na.rm = TRUE),
    TPU_sd = sd(TPU.a, na.rm = TRUE),
    TPU_se = TPU_sd / sqrt(TPU_n),
    
    Rd_n = sum(!is.na(Rdark.a)),
    Rd_mean = mean(Rdark.a, na.rm = TRUE),
    Rd_sd = sd(Rdark.a, na.rm = TRUE),
    Rd_se = Rd_sd / sqrt(Rd_n),
    
    cica_n = sum(!is.na(Ci.Ca)),
    cica_mean = mean(Ci.Ca, na.rm = TRUE),
    cica_sd = sd(Ci.Ca, na.rm = TRUE),
    cica_se = cica_sd / sqrt(cica_n),
    
    lma_n = sum(!is.na(LMA)),
    lma_mean = mean(LMA, na.rm = TRUE),
    lma_sd = sd(LMA, na.rm = TRUE),
    lma_se = lma_sd / sqrt(lma_n),
    
    nmass_n = sum(!is.na(LeafN)),
    nmass_mean = mean(LeafN, na.rm = TRUE),
    nmass_sd = sd(LeafN, na.rm = TRUE),
    nmass_se = nmass_sd / sqrt(nmass_n),
    
    narea_n = sum(!is.na(N.area)),
    narea_mean = mean(N.area, na.rm = TRUE),
    narea_sd = sd(N.area, na.rm = TRUE),
    narea_se = narea_sd / sqrt(narea_n),

    pmass_n = sum(!is.na(LeafP)),
    pmass_mean = mean(LeafP, na.rm = TRUE),
    pmass_sd = sd(LeafP, na.rm = TRUE),
    pmass_se = pmass_sd / sqrt(pmass_n),
    
    parea_n = sum(!is.na(P.area)),
    parea_mean = mean(P.area, na.rm = TRUE),
    parea_sd = sd(P.area, na.rm = TRUE),
    parea_se = parea_sd / sqrt(parea_n),
    
    leafnp_n = sum(!is.na(RatioNP)),
    leafnp_mean = mean(RatioNP, na.rm = TRUE),
    leafnp_sd = sd(RatioNP, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(leafnp_n),
    
    pnue_n = sum(!is.na(PNUE)),
    pnue_mean = mean(PNUE, na.rm = TRUE),
    pnue_sd = sd(PNUE, na.rm = TRUE),
    pnue_se = pnue_sd / sqrt(pnue_n),
    
    ppue_n = sum(!is.na(PPUE)),
    ppue_mean = mean(PPUE, na.rm = TRUE),
    ppue_sd = sd(PPUE, na.rm = TRUE),
    ppue_se = ppue_sd / sqrt(ppue_n))

# # Prep for easy merge into compiled datasheet
bloomfield_data_summary_control <- bloomfield_data_summary %>%
  ungroup(Treatment) %>%
  filter(Treatment == "Control") %>%
  mutate(exp = "bloomfield2014") %>%
  dplyr::select(-Treatment)
names(bloomfield_data_summary_control)[2:61] <- str_c(names(bloomfield_data_summary_control)[2:61], "_control")

bloomfield_data_summary_treatment <- bloomfield_data_summary %>%
  ungroup(Treatment) %>%
  filter(Treatment != "Control") %>%
  mutate(exp = "bloomfield2014")
names(bloomfield_data_summary_treatment)[3:62] <- str_c(names(bloomfield_data_summary_treatment)[3:62], "_trt")


bloomfield_data_summary_control %>%
  full_join(bloomfield_data_summary_treatment, by = c("exp", "Months.in")) %>%
  dplyr::select(exp, Treatment, Months.in, 
                
                Amax_mean_control, Amax_mean_trt, 
                Amax_sd_control, Amax_sd_trt, 
                Amax_se_control, Amax_se_trt, 
                Amax_n_control, Amax_n_trt,
              
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
                
                TPU_mean_control, TPU_mean_trt, 
                TPU_sd_control, TPU_sd_trt, 
                TPU_se_control, TPU_se_trt, 
                TPU_n_control, TPU_n_trt,
                
                Rd_mean_control, Rd_mean_trt, 
                Rd_sd_control, Rd_sd_trt, 
                Rd_se_control, Rd_se_trt, 
                Rd_n_control, Rd_n_trt,
                
                cica_mean_control, cica_mean_trt, 
                cica_sd_control, cica_sd_trt, 
                cica_se_control, cica_se_trt,
                cica_n_control, cica_n_trt,
                
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
                
                pnue_mean_control, pnue_mean_trt, 
                pnue_sd_control, pnue_sd_trt, 
                pnue_se_control, pnue_se_trt, 
                pnue_n_control, pnue_n_trt,
                
                ppue_mean_control, ppue_mean_trt, 
                ppue_sd_control, ppue_sd_trt, 
                ppue_se_control, ppue_se_trt, 
                ppue_n_control, ppue_n_trt)  %>%
  pivot_longer(cols = Amax_mean_control:ppue_n_trt,
               names_to = c("trait", "stat", "trt"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:trt,
              values_from = value) %>%
  mutate(trait = factor(trait, levels = c("Amax", "Asat", "Vcmax", "Jmax",
                                          "TPU", "Rd", "cica", "lma",
                                          "nmass", "narea", "pmass", "parea",
                                          "leafnp", "pnue", "ppue"))) %>%
  arrange(trait) %>%
  write.csv("../data_summaries/Bloomfield2014_FPB_summarized.csv", row.names = F)
  