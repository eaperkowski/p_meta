# summarize digitized plot data

# Libraries
library(tidyverse)

###############################################################################
# Chen et al., 2024
###############################################################################

# Read file and mutate anet to be asat (to align with other data)
chen_2024_digiData <- read.csv("../raw_data/Chen_2024_jxb_digitized.csv") %>%
  mutate(variable = ifelse(variable == "anet", "asat", variable),
         treatment = ifelse(group_id == "p40_wild", "P", "C")) %>%
  dplyr::select(trait = variable, treatment, mean:n, se)
head(chen_2024_digiData)

# Prep for easy merge into compiled datasheet
chen_2024_data_summary_control <- chen_2024_digiData %>%
  filter(treatment == "C") %>%
  mutate(exp = "Chen2024") %>%
  select(-treatment)
names(chen_2024_data_summary_control)[2:5] <- str_c(names(chen_2024_data_summary_control)[2:5], "_control")

chen_2024_data_summary_treatment <- chen_2024_digiData %>%
  filter(treatment != "C") %>%
  mutate(exp = "Chen2024")
names(chen_2024_data_summary_treatment)[3:6] <- str_c(names(chen_2024_data_summary_treatment)[3:6], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
chen_2024_data_summary_control %>%
  full_join(chen_2024_data_summary_treatment, by = c("exp", "trait")) %>%
  dplyr::select(trait, treatment, mean_control, mean_trt, 
                sd_control, sd_trt, se_control, se_trt,
                n_control, n_trt) #%>%
#  write.csv("../data_summaries/Chen2024_summary.csv", row.names = F)


###############################################################################
# Ye et al. 2022
###############################################################################

# Read file and mutate anet to be asat (to align with other data)
Ye2022_digiData <- read.csv("../raw_data/Ye_2022_jfr_digitized.csv") %>%
  separate(group_id, into = c("spp", "treatment")) %>%
  dplyr::select(trait = variable, spp, treatment, mean:n, se)
head(Ye2022_digiData)

# Prep for easy merge into compiled datasheet
Ye2022_data_summary_control <- Ye2022_digiData %>%
  filter(treatment == "control") %>%
  mutate(exp = "Ye2022") %>%
  select(-treatment)
names(Ye2022_data_summary_control)[3:6] <- str_c(names(Ye2022_data_summary_control)[3:6], "_control")

Ye2022_data_summary_treatment <- Ye2022_digiData %>%
  filter(treatment != "control") %>%
  mutate(exp = "Ye2022")
names(Ye2022_data_summary_treatment)[4:7] <- str_c(names(Ye2022_data_summary_treatment)[4:7], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
Ye2022_data_summary_control %>%
  full_join(Ye2022_data_summary_treatment, by = c("exp", "trait", "spp")) %>%
  dplyr::select(trait, spp, treatment, mean_control, mean_trt, 
                sd_control, sd_trt, se_control, se_trt,
                n_control, n_trt)# %>%
 # write.csv("../data_summaries/Ye2022_summary.csv", row.names = F)
