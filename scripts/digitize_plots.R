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

# Davidson et al. 2004
davidson2004 <- metaDigitise(dir = "../plots_to_digitize/Davidson_2004/")

# Dong et al. 2016
dong2016 <- metaDigitise(dir = "../plots_to_digitize/Dong_2016/")

# Eller et al. 2017
eller2017 <- metaDigitise(dir = "../plots_to_digitize/Eller_2017/")

# Falk et al. 2010
falk2010 <- metaDigitise(dir = "../plots_to_digitize/Falk_2010/")

# Fornara et al. 2013
fornara2013 <- metaDigitise(dir = "../plots_to_digitize/Fornara_2013/")

# Frost et al. 2009
frost2009 <- metaDigitise(dir = "../plots_to_digitize/Frost_2009/")

# Gough et al. 2003
gough2003 <- metaDigitise(dir = "../plots_to_digitize/Gough_2003/")

# Gusewell 2002
gusewell2002 <- metaDigitise(dir = "../plots_to_digitize/Güsewell_2002/")

# Gusewell 2003
gusewell2003 <- metaDigitise(dir = "../plots_to_digitize/Güsewell_2003/")

# Harrington 2001
harrington2001 <- metaDigitise(dir = "../plots_to_digitize/Harrington_2001/")

# Haubensak 2011
haubensak2011 <- metaDigitise(dir = "../plots_to_digitize/Haubensak_2011/")

# Huff 2015
huff2015 <- metaDigitise(dir = "../plots_to_digitize/Huff_2015/")

# Ket 2011
ket2011 <- metaDigitise(dir = "../plots_to_digitize/Ket_2011/")

# Lawrence 2001
lawrence2001 <- metaDigitise(dir = "../plots_to_digitize/Lawrence_2001/")

# Li 2011
li2011 <- metaDigitise(dir = "../plots_to_digitize/Li_2011/")

# Li 2014
li2014 <- metaDigitise(dir = "../plots_to_digitize/Li_2014/")

# Li 2015
li2015 <- metaDigitise(dir = "../plots_to_digitize/Li_2015")

# Ludwig 2001
ludwig2001 <- metaDigitise(dir = "../plots_to_digitize/Ludwig_2001/")

# Lund 2009
lund2009 <- metaDigitise(dir = "../plots_to_digitize/Lund_2009/")

# McMaster 1982
mcmaster1982 <- metaDigitise(dir = "../plots_to_digitize/McMaster_1982/")

# Ngai 2004
ngai2004 <- metaDigitise(dir = "../plots_to_digitize/Ngai_2004/")

# Ngatia 2015
ngatia2015 <- metaDigitise(dir = "../plots_to_digitize/Ngatia_2015/")

# Nielsen 2009
nielsen2009 <- metaDigitise(dir = "../plots_to_digitize/Nielsen_2009/")

nielsen2009 %>%
  separate(group_id, into = c("species", "trt")) %>%
  mutate(trt = ifelse(trt == "nnp", "np", trt),
         trt = factor(trt, levels = c("control", "n", "p", "np")),
         mean = round(mean, 1),
         sd = round(sd, 1),
         se = round(se, 1)) %>%
  arrange(species, variable, trt) %>%
  dplyr::select(variable:mean, se)

