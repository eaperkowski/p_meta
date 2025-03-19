# Libraries
library(metaDigitise)
library(tidyverse)

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

# Aerts et al 2003
aerts2003 <- metaDigitise(dir = "../plots_to_digitize/Aerts_2003/")
# write.csv(aerts2003, "../raw_data/Aerts_2003_oikos.csv", row.names = F)

# Arens et al. 2008
arens2008 <- metaDigitise(dir = "../plots_to_digitize/Arens_2008/")
# write.csv(arens2008, "../raw_data/Arens_2008.csv", row.names = F)

# Aydin & Uzun 2004
aydin2004 <- metaDigitise(dir = "../plots_to_digitize/Aydin_Uzun_2004/")

# Blanke et al. 2012
blanke2012 <- metaDigitise(dir = "../plots_to_digitize/Blanke_2012/")

# Borer et al. 2013
borer2013 <- metaDigitise(dir = "../plots_to_digitize/Borer_2013/")

# Bowman et al. 1993
bowman1993 <- metaDigitise(dir = "../plots_to_digitize/Bowman_1993/")
bowman1993
