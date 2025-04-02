# Determine Bloomfield et al. (2014) treatment mean +/- SE,
# parse by species instead of experiment

# Libraries
library(tidyverse)

# Add temp_standardize fxn to calculate Vcmax25 and Jmax25
source("../../helper_fxns/temp_standardize.R")

# Load Bloomfield et al. (2014) data and change varname for treatment
# such that Control (in original datasheet) is the +P treatment and
# minus P is the Control treatment
bloomfield_data <- read.csv("../../raw_data/Bloomfield_2014_FPB.csv") %>%
  mutate(Treatment = ifelse(Treatment == "Control", "P", 
                            ifelse(Treatment == "Minus P", "Control", NA)),
         Treatment = factor(Treatment, levels = c("Control", "P")),
         Vcmax25.a = temp_standardize(Vcmax.a, "Vcmax", 25, tLeaf = avgTleaf, tGrow = avgTleaf),
         Jmax25.a = temp_standardize(Jmax.a, "Jmax", 25, tLeaf = avgTleaf, tGrow = avgTleaf),
         Jmax25.Vcmax25 = Jmax25.a / Vcmax25.a,
         Rdark25 = temp_standardize(Rdark.a, "Rd", 25, pft = "GM", tLeaf = avgTleaf, tGrow = avgTleaf))

head(bloomfield_data)

# How many unique species and growth forms?
unique(bloomfield_data$Species)
unique(bloomfield_data$Adult.stature)
unique(bloomfield_data$Harvest)

# Group by treatment and timepoint, then summarize Bloomfield et al. (2014)
# for easy input into compiled datasheet
bloomfield_data_summary <- bloomfield_data %>%
  group_by(Treatment, Species, Months.in) %>%
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
    
    Vcmax25_n = sum(!is.na(Vcmax25.a)),
    Vcmax25_mean = mean(Vcmax25.a, na.rm = TRUE),
    Vcmax25_sd = sd(Vcmax25.a, na.rm = TRUE),
    Vcmax25_se = Vcmax25_sd / sqrt(Vcmax25_n),
    
    Jmax_n = sum(!is.na(Jmax.a)),
    Jmax_mean = mean(Jmax.a, na.rm = TRUE),
    Jmax_sd = sd(Jmax.a, na.rm = TRUE),
    Jmax_se = Jmax_sd / sqrt(Jmax_n),
    
    Jmax25_n = sum(!is.na(Jmax25.a)),
    Jmax25_mean = mean(Jmax25.a, na.rm = TRUE),
    Jmax25_sd = sd(Jmax25.a, na.rm = TRUE),
    Jmax25_se = Jmax25_sd / sqrt(Jmax25_n),
    
    Jmax25Vcmax25_n = sum(!is.na(Jmax25.Vcmax25)),
    Jmax25Vcmax25_mean = mean(Jmax25.Vcmax25, na.rm = TRUE),
    Jmax25Vcmax25_sd = sd(Jmax25.Vcmax25, na.rm = TRUE),
    Jmax25Vcmax25_se = Jmax25Vcmax25_sd / sqrt(Jmax25Vcmax25_n),
    
    TPU_n = sum(!is.na(TPU.a)),
    TPU_mean = mean(TPU.a, na.rm = TRUE),
    TPU_sd = sd(TPU.a, na.rm = TRUE),
    TPU_se = TPU_sd / sqrt(TPU_n),
    
    Rd_n = sum(!is.na(Rdark.a)),
    Rd_mean = mean(Rdark.a, na.rm = TRUE),
    Rd_sd = sd(Rdark.a, na.rm = TRUE),
    Rd_se = Rd_sd / sqrt(Rd_n),
    
    Rd25_n = sum(!is.na(Rdark25)),
    Rd25_mean = mean(Rdark25, na.rm = TRUE),
    Rd25_sd = sd(Rdark25, na.rm = TRUE),
    Rd25_se = Rd25_sd / sqrt(Rd25_n),
    
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
head(bloomfield_data_summary)


# # Prep for easy merge into compiled datasheet
bloomfield_data_summary_control <- bloomfield_data_summary %>%
  ungroup(Treatment) %>%
  filter(Treatment == "Control") %>%
  dplyr::select(-Treatment)
names(bloomfield_data_summary_control)[3:78] <- str_c(names(bloomfield_data_summary_control)[3:78], "_control")

bloomfield_data_summary_treatment <- bloomfield_data_summary %>%
  ungroup(Treatment) %>%
  filter(Treatment != "Control")
names(bloomfield_data_summary_treatment)[4:79] <- str_c(names(bloomfield_data_summary_treatment)[4:79], "_trt")


bloomfield_data_summary_control %>%
  full_join(bloomfield_data_summary_treatment, by = c("Species", "Months.in")) %>%
  ungroup(Species) %>%
  dplyr::select(spp = Species, Treatment, Months.in, 
                
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
                
                Vcmax25_mean_control, Vcmax25_mean_trt, 
                Vcmax25_sd_control, Vcmax25_sd_trt, 
                Vcmax25_se_control, Vcmax25_se_trt, 
                Vcmax25_n_control, Vcmax25_n_trt,
                
                Jmax_mean_control, Jmax_mean_trt, 
                Jmax_sd_control, Jmax_sd_trt, 
                Jmax_se_control, Jmax_se_trt, 
                Jmax_n_control, Jmax_n_trt,
                
                Jmax25_mean_control, Jmax25_mean_trt, 
                Jmax25_sd_control, Jmax25_sd_trt, 
                Jmax25_se_control, Jmax25_se_trt, 
                Jmax25_n_control, Jmax25_n_trt,
                
                Jmax25Vcmax25_mean_control, Jmax25Vcmax25_mean_trt, 
                Jmax25Vcmax25_sd_control, Jmax25Vcmax25_sd_trt, 
                Jmax25Vcmax25_se_control, Jmax25Vcmax25_se_trt, 
                Jmax25Vcmax25_n_control, Jmax25Vcmax25_n_trt,
                
                TPU_mean_control, TPU_mean_trt, 
                TPU_sd_control, TPU_sd_trt, 
                TPU_se_control, TPU_se_trt, 
                TPU_n_control, TPU_n_trt,
                
                Rd_mean_control, Rd_mean_trt, 
                Rd_sd_control, Rd_sd_trt, 
                Rd_se_control, Rd_se_trt, 
                Rd_n_control, Rd_n_trt,
                
                Rd25_mean_control, Rd25_mean_trt, 
                Rd25_sd_control, Rd25_sd_trt, 
                Rd25_se_control, Rd25_se_trt, 
                Rd25_n_control, Rd25_n_trt,
                
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
  mutate(trait = factor(trait, levels = c("Amax", "Asat", "Vcmax", "Vcmax25",
                                          "Jmax", "Jmax25", "Jmax25Vcmax25",
                                          "TPU", "Rd", "Rd25", "cica", "lma",
                                          "nmass", "narea", "pmass", "parea",
                                          "leafnp", "pnue", "ppue"))) %>%
  filter(n_control > 1 & n_trt > 1) %>%
  filter(n_distinct(Treatment) == 3) %>%
  
  
  
  arrange(trait) %>%
  write.csv("../../data_summaries/species_level/Bloomfield2014_summary_spp.csv", row.names = F)
  