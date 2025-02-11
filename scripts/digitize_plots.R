# Libraries
library(metaDigitise)

# Chen et al. 2024
chen2024 <- metaDigitise(dir = "../plots_to_digitize/Chen_2024_jxb/")
## write.csv(chen2024, "../raw_data/Chen_2024_jxb_digitized.csv", row.names = F)

# Ye et al. 2022
ye2022 <- metaDigitise(dir = "../plots_to_digitize/Ye_2022_jfr/")
## write.csv(ye2022, "../raw_data/Ye_2022_jfr_digitized.csv", row.names = F)

# Black et al. 2000
black2000 <- metaDigitise(dir = "../plots_to_digitize/Black_2000_pce/")
# write.csv(black2000, "../raw_data/Black_2000_pce_digitized.csv", row.names = F)

# Hayes et al. 2022
hayes2022 <- metaDigitise(dir = "../plots_to_digitize/Hayes_2022_aob/")
# write.csv(hayes2022, "../raw_data/Hayes_2022_aob.csv", row.names = F)
