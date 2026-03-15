# Libraries
library(tidyverse)
library(metafor)
library(ggpubr)
library(naniar) # to resolve NA/<NA> issue
library(orchaRd)

# Load meta-analysis results
meta_results_int <- read.csv("../data/CNPmeta_logr_results_int.csv") %>%
  mutate(myc_nas = ifelse(myc_assoc == "NM" | myc_assoc == "AM" |
                            myc_assoc == "NM-AM", 
                          "scavenging",
                          ifelse(myc_assoc == "EcM" | myc_assoc == "ErM" | 
                                   myc_assoc == "EcM-AM", "mining", NA)))

# Check file structure
head(meta_results_int)

##############################################################################
# Marea climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "lma" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_marea_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "lma" & 
                                    !is.na(gs_mat) & dNPi < 2))

# Summary
int_marea_clim_summary <- data.frame(trait = "marea",
                                      nut_add = "int",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(int_marea_clim)),
                                      row.names = NULL)

##############################################################################
# Nmass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_n_mass" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_nmass_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_n_mass" & 
                                    !is.na(gs_mat) & dNPi > -2))

# Summary
int_nmass_clim_summary <- data.frame(trait = "nmass",
                                      nut_add = "int",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(int_nmass_clim)),
                                      row.names = NULL)

##############################################################################
# Narea climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_n_area" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_narea_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_n_area" & 
                                    !is.na(gs_mat)))

# Summary
int_narea_clim_summary <- data.frame(trait = "narea",
                                      nut_add = "int",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(int_narea_clim)),
                                      row.names = NULL)

##############################################################################
# Pmass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_p_mass" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_pmass_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_p_mass" & 
                                    !is.na(gs_mat)))

# Summary
int_pmass_clim_summary <- data.frame(trait = "pmass",
                                     nut_add = "int",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_pmass_clim)),
                                     row.names = NULL)

##############################################################################
# Parea climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_p_area" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_parea_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_p_area" & 
                                    !is.na(gs_mat) & dNPi < 2))

# Summary
int_parea_clim_summary <- data.frame(trait = "parea",
                                     nut_add = "int",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_parea_clim)),
                                     row.names = NULL)

##############################################################################
# Leaf N:P climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_np" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_leafnp_clim <- rma.mv(yi = dNPi,
                             V = vNPi,
                             W = wNPi,
                             method = "REML", 
                             random = ~ 1 | exp, 
                             mods = ~ gs_mat + gs_ai + gs_par,
                             slab = exp, 
                             control = list(stepadj = 0.3), 
                             data = meta_results_int %>% 
                               filter(response == "leaf_np" & 
                                        !is.na(gs_mat)))

# Summary
int_leafnp_clim_summary <- data.frame(trait = "leaf_np",
                                         nut_add = "int",
                                         mod = c("intrcpt", "gs_mat",
                                                 "gs_ai", "gs_par"),
                                         coef(summary(int_leafnp_clim)),
                                         row.names = NULL)


##############################################################################
# Total biomass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "tbio_gm2" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_tbio_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "tbio_gm2" & 
                                    !is.na(gs_mat)))

# Summary
int_tbio_clim_summary <- data.frame(trait = "tbio_gm2",
                                     nut_add = "int",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_tbio_clim)),
                                     row.names = NULL)

##############################################################################
# Aboveground biomass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "anpp" & !is.na(gs_mat) & dNPi < 2 & dNPi > -2) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_anpp_clim <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "anpp" & 
                                   !is.na(gs_mat) & dNPi < 2 & dNPi > -2))

# Summary
int_anpp_clim_summary <- data.frame(trait = "anpp",
                                    nut_add = "int",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(int_anpp_clim)),
                                    row.names = NULL)

##############################################################################
# Belowground biomass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "bnpp" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_bnpp_clim <- rma.mv(yi = dNPi,
                       V = vNPi,
                       W = wNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ gs_mat + gs_ai + gs_par,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(response == "bgb" & 
                                  !is.na(gs_mat)))

# Summary
int_bnpp_clim_summary <- data.frame(trait = "bnpp",
                                   nut_add = "int",
                                   mod = c("intrcpt", "gs_mat",
                                           "gs_ai", "gs_par"),
                                   coef(summary(int_bnpp_clim)),
                                   row.names = NULL)

##############################################################################
# Root mass fraction climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "rmf" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_rmf_clim <- rma.mv(yi = dNPi,
                       V = vNPi,
                       W = wNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ gs_mat + gs_ai + gs_par,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(response == "rmf" & 
                                  !is.na(gs_mat)))

# Summary
int_rmf_clim_summary <- data.frame(trait = "rmf",
                                   nut_add = "int",
                                   mod = c("intrcpt", "gs_mat",
                                           "gs_ai", "gs_par"),
                                   coef(summary(int_rmf_clim)),
                                   row.names = NULL)

##############################################################################
# Root:shoot climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "rootshoot" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_rootshoot_clim <- rma.mv(yi = dNPi,
                             V = vNPi,
                             W = wNPi,
                             method = "REML", 
                             random = ~ 1 | exp, 
                             mods = ~ gs_mat + gs_ai + gs_par,
                             slab = exp, 
                             control = list(stepadj = 0.3), 
                             data = meta_results_int %>% 
                               filter(response == "rootshoot" & 
                                        !is.na(gs_mat)))

# Summary
int_rootshoot_clim_summary <- data.frame(trait = "rootshoot",
                                   nut_add = "int",
                                   mod = c("intrcpt", "gs_mat",
                                           "gs_ai", "gs_par"),
                                   coef(summary(int_rootshoot_clim)),
                                   row.names = NULL)

##############################################################################
# Merge climate moderators and write to .csv
##############################################################################
int_marea_clim_summary %>%
  full_join(int_nmass_clim_summary) %>%
  full_join(int_narea_clim_summary) %>%
  full_join(int_pmass_clim_summary) %>%
  full_join(int_parea_clim_summary) %>%
  full_join(int_leafnp_clim_summary) %>%
  full_join(int_tbio_clim_summary) %>%
  full_join(int_anpp_clim_summary) %>%
  full_join(int_bnpp_clim_summary) %>%
  full_join(int_rmf_clim_summary) %>%
  full_join(int_rootshoot_clim_summary) %>%
  mutate(estimate_perc = (exp(estimate) - 1) * 100,
         se_perc = (exp(se) - 1) * 100,
         ci.lb_perc = (exp(ci.lb) - 1) * 100,
         ci.ub_perc = (exp(ci.ub) - 1) * 100,
         across(estimate:ci.ub_perc, ~round(.x, 3)),
         estimate_se = str_c(estimate, "±", se),
         estimate_se_perc = str_c(estimate_perc, "±", se_perc),
         ci.range = str_c("[", ci.lb, ", ", ci.ub, "]"),
         ci.range_perc = str_c("[", ci.lb_perc, ", ", ci.ub_perc, "]")) %>%
  write.csv("../data/CNPmeta_clim_moderators_int.csv", row.names = F)

##############################################################################
# Marea photosynthetic pathway
##############################################################################

# Model
int_marea_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "lma" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_marea_photo <- data.frame(trait = "marea", 
                               nut_add = "int",
                               mod = "photo",
                               mod_results(int_marea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(int_marea_pft))[2,3],
                               p = coef(summary(int_marea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_marea_myc <- data.frame(trait = "marea", 
                             nut_add = "int",
                             mod = "myc_nas",
                             mod_results(int_marea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_marea_pft))[3,3],
                             p = coef(summary(int_marea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_marea_nfix <- data.frame(trait = "marea", 
                             nut_add = "int",
                             mod = "nfix",
                             mod_results(int_marea_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_marea_pft))[4,3],
                             p = coef(summary(int_marea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_marea_pft_results <- int_marea_photo %>% 
  rbind(int_marea_myc) %>% 
  rbind(int_marea_nfix) %>%
  mutate(k = 113,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Nmass
##############################################################################

# Model
int_nmass_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "leaf_n_mass" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_nmass_photo <- data.frame(trait = "nmass", 
                              nut_add = "int",
                              mod = "photo",
                              mod_results(int_nmass_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(int_nmass_pft))[2,3],
                              p = coef(summary(int_nmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_nmass_myc <- data.frame(trait = "nmass", 
                            nut_add = "int",
                            mod = "myc_nas",
                            mod_results(int_nmass_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_nmass_pft))[3,3],
                            p = coef(summary(int_nmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_nmass_nfix <- data.frame(trait = "nmass", 
                             nut_add = "int",
                             mod = "nfix",
                             mod_results(int_nmass_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_nmass_pft))[4,3],
                             p = coef(summary(int_nmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_nmass_pft_results <- int_nmass_photo %>% 
  rbind(int_nmass_myc) %>% 
  rbind(int_nmass_nfix) %>%
  mutate(k = 136,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)


##############################################################################
# Narea - pft
##############################################################################

# Model
int_narea_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "leaf_n_area" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_narea_photo <- data.frame(trait = "narea", 
                              nut_add = "int",
                              mod = "photo",
                              mod_results(int_narea_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(int_narea_pft))[2,3],
                              p = coef(summary(int_narea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_narea_myc <- data.frame(trait = "narea", 
                            nut_add = "int",
                            mod = "myc_nas",
                            mod_results(int_narea_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_narea_pft))[3,3],
                            p = coef(summary(int_narea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_narea_nfix <- data.frame(trait = "narea", 
                             nut_add = "int",
                             mod = "nfix",
                             mod_results(int_narea_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_narea_pft))[4,3],
                             p = coef(summary(int_narea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_narea_pft_results <- int_narea_photo %>% 
  rbind(int_narea_myc) %>% 
  rbind(int_narea_nfix) %>%
  mutate(k = 87,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Pmass - plant functional type
##############################################################################
# Model
int_pmass_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "leaf_p_mass" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_pmass_photo <- data.frame(trait = "pmass", 
                              nut_add = "int",
                              mod = "photo",
                              mod_results(int_pmass_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(int_pmass_pft))[2,3],
                              p = coef(summary(int_pmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_pmass_myc <- data.frame(trait = "pmass", 
                            nut_add = "int",
                            mod = "myc_nas",
                            mod_results(int_pmass_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_pmass_pft))[3,3],
                            p = coef(summary(int_pmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_pmass_nfix <- data.frame(trait = "pmass", 
                             nut_add = "int",
                             mod = "nfix",
                             mod_results(int_pmass_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_pmass_pft))[4,3],
                             p = coef(summary(int_pmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_pmass_pft_results <- int_pmass_photo %>% 
  rbind(int_pmass_myc) %>% 
  rbind(int_pmass_nfix) %>%
  mutate(k = 130,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

######### stopped here #############

##############################################################################
# Parea - plant functional type
##############################################################################

##############
# N addition
##############
nadd_parea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "leaf_p_area" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_parea_photo <- data.frame(trait = "parea", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_parea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_parea_pft))[2,3],
                               p = coef(summary(nadd_parea_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_parea_myc <- data.frame(trait = "parea", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_parea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_parea_pft))[3,3],
                             p = coef(summary(nadd_parea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_parea_nfix <- data.frame(trait = "parea", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_parea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_parea_pft))[4,3],
                              p = coef(summary(nadd_parea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_parea_pft_results <- nadd_parea_photo %>% 
  rbind(nadd_parea_myc) %>% 
  rbind(nadd_parea_nfix)

##############
# P addition
##############
padd_parea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "leaf_p_area" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_parea_photo <- data.frame(trait = "parea", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_parea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_parea_pft))[2,3],
                               p = coef(summary(padd_parea_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_parea_myc <- data.frame(trait = "parea", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_parea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_parea_pft))[3,3],
                             p = coef(summary(padd_parea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_parea_nfix <- data.frame(trait = "parea", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_parea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_parea_pft))[4,3],
                              p = coef(summary(padd_parea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_parea_pft_results <- padd_parea_photo %>% 
  rbind(padd_parea_myc) %>% 
  rbind(padd_parea_nfix)

################
# N+P addition
################
npadd_parea_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_parea_photo <- data.frame(trait = "parea", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_parea_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_parea_pft))[2,3],
                                p = coef(summary(npadd_parea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_parea_myc <- data.frame(trait = "parea", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_parea_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_parea_pft))[3,3],
                              p = coef(summary(npadd_parea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_parea_nfix <- data.frame(trait = "parea", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_parea_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_parea_pft))[4,3],
                               p = coef(summary(npadd_parea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_parea_pft_results <- npadd_parea_photo %>% 
  rbind(npadd_parea_myc) %>% 
  rbind(npadd_parea_nfix)

#############
# Merge Parea moderator results, with some light cleaning
#############
parea_pft_summary <- rbind(nadd_parea_pft_results, 
                           padd_parea_pft_results, 
                           npadd_parea_pft_results) %>%
  mutate(k = 82,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Leaf N:P - plant functional type
##############################################################################

##############
# N addition
##############
nadd_leafnp_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_np" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_leafnp_photo <- data.frame(trait = "leaf_np", 
                                nut_add = "n",
                                mod = "photo",
                                mod_results(nadd_leafnp_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(nadd_leafnp_pft))[2,3],
                                p = coef(summary(nadd_leafnp_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_leafnp_myc <- data.frame(trait = "leaf_np", 
                              nut_add = "n",
                              mod = "myc_nas",
                              mod_results(nadd_leafnp_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_leafnp_pft))[3,3],
                              p = coef(summary(nadd_leafnp_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_leafnp_nfix <- data.frame(trait = "leaf_np", 
                               nut_add = "n",
                               mod = "nfix",
                               mod_results(nadd_leafnp_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_leafnp_pft))[4,3],
                               p = coef(summary(nadd_leafnp_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_leafnp_pft_results <- nadd_leafnp_photo %>% 
  rbind(nadd_leafnp_myc) %>% 
  rbind(nadd_leafnp_nfix)

##############
# P addition
##############
padd_leafnp_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_np" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_leafnp_photo <- data.frame(trait = "leaf_np", 
                                nut_add = "p",
                                mod = "photo",
                                mod_results(padd_leafnp_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(padd_leafnp_pft))[2,3],
                                p = coef(summary(padd_leafnp_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_leafnp_myc <- data.frame(trait = "leaf_np", 
                              nut_add = "p",
                              mod = "myc_nas",
                              mod_results(padd_leafnp_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_leafnp_pft))[3,3],
                              p = coef(summary(padd_leafnp_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_leafnp_nfix <- data.frame(trait = "leaf_np", 
                               nut_add = "p",
                               mod = "nfix",
                               mod_results(padd_leafnp_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_leafnp_pft))[4,3],
                               p = coef(summary(padd_leafnp_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_leafnp_pft_results <- padd_leafnp_photo %>% 
  rbind(padd_leafnp_myc) %>% 
  rbind(padd_leafnp_nfix)

################
# N+P addition
################
npadd_leafnp_pft <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path + myc_nas + n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_np" & 
                                      !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_leafnp_photo <- data.frame(trait = "leaf_np", 
                                 nut_add = "np",
                                 mod = "photo",
                                 mod_results(npadd_leafnp_pft, 
                                             mod = "photo_path", 
                                             group = "exp")$mod_table,
                                 z = coef(summary(npadd_leafnp_pft))[2,3],
                                 p = coef(summary(npadd_leafnp_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_leafnp_myc <- data.frame(trait = "leaf_np", 
                               nut_add = "np",
                               mod = "myc_nas",
                               mod_results(npadd_leafnp_pft, 
                                           mod = "myc_nas", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_leafnp_pft))[3,3],
                               p = coef(summary(npadd_leafnp_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_leafnp_nfix <- data.frame(trait = "leaf_np", 
                                nut_add = "np",
                                mod = "nfix",
                                mod_results(npadd_leafnp_pft, 
                                            mod = "n_fixer", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_leafnp_pft))[4,3],
                                p = coef(summary(npadd_leafnp_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_leafnp_pft_results <- npadd_leafnp_photo %>% 
  rbind(npadd_leafnp_myc) %>% 
  rbind(npadd_leafnp_nfix)

#############
# Merge leaf N:P moderator results, with some light cleaning
#############
leafnp_pft_summary <- rbind(nadd_leafnp_pft_results, 
                            padd_leafnp_pft_results, 
                            npadd_leafnp_pft_results) %>%
  mutate(k = 115,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Asat - plant functional type
##############################################################################

##############
# N addition
##############
nadd_asat_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "asat" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_asat_photo <- data.frame(trait = "asat", 
                              nut_add = "n",
                              mod = "photo",
                              mod_results(nadd_asat_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_asat_pft))[2,3],
                              p = coef(summary(nadd_asat_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_asat_myc <- data.frame(trait = "asat", 
                            nut_add = "n",
                            mod = "myc_nas",
                            mod_results(nadd_asat_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_asat_pft))[3,3],
                            p = coef(summary(nadd_asat_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_asat_nfix <- data.frame(trait = "asat", 
                             nut_add = "n",
                             mod = "nfix",
                             mod_results(nadd_asat_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_asat_pft))[4,3],
                             p = coef(summary(nadd_asat_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_asat_pft_results <- nadd_asat_photo %>% 
  rbind(nadd_asat_myc) %>% 
  rbind(nadd_asat_nfix)

##############
# P addition
##############
padd_asat_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "asat" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_asat_photo <- data.frame(trait = "asat", 
                              nut_add = "p",
                              mod = "photo",
                              mod_results(padd_asat_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_asat_pft))[2,3],
                              p = coef(summary(padd_asat_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_asat_myc <- data.frame(trait = "asat", 
                            nut_add = "p",
                            mod = "myc_nas",
                            mod_results(padd_asat_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_asat_pft))[3,3],
                            p = coef(summary(padd_asat_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_asat_nfix <- data.frame(trait = "asat", 
                             nut_add = "p",
                             mod = "nfix",
                             mod_results(padd_asat_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_asat_pft))[4,3],
                             p = coef(summary(padd_asat_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_asat_pft_results <- padd_asat_photo %>% 
  rbind(padd_asat_myc) %>% 
  rbind(padd_asat_nfix)

################
# N+P addition
################
npadd_asat_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "asat" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_asat_photo <- data.frame(trait = "asat", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_asat_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_asat_pft))[2,3],
                               p = coef(summary(npadd_asat_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_asat_myc <- data.frame(trait = "asat", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_asat_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_asat_pft))[3,3],
                             p = coef(summary(npadd_asat_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_asat_nfix <- data.frame(trait = "asat", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_asat_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_asat_pft))[4,3],
                              p = coef(summary(npadd_asat_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_asat_pft_results <- npadd_asat_photo %>% 
  rbind(npadd_asat_myc) %>% 
  rbind(npadd_asat_nfix)

#############
# Merge Asat moderator results, with some light cleaning
#############
asat_pft_summary <- rbind(nadd_asat_pft_results, 
                          padd_asat_pft_results, 
                          npadd_asat_pft_results) %>%
  mutate(k = 93,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# gsw - plant functional type
##############################################################################

##############
# N addition
##############
nadd_gsw_pft <- rma.mv(logr, 
                       logr_var,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results %>% 
                         filter(manip_type == "n" & 
                                  myvar == "gsw" & 
                                  !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_gsw_photo <- data.frame(trait = "gsw", 
                             nut_add = "n",
                             mod = "photo",
                             mod_results(nadd_gsw_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_gsw_pft))[2,3],
                             p = coef(summary(nadd_gsw_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_gsw_myc <- data.frame(trait = "gsw", 
                           nut_add = "n",
                           mod = "myc_nas",
                           mod_results(nadd_gsw_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(nadd_gsw_pft))[3,3],
                           p = coef(summary(nadd_gsw_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_gsw_nfix <- data.frame(trait = "gsw", 
                            nut_add = "n",
                            mod = "nfix",
                            mod_results(nadd_gsw_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_gsw_pft))[4,3],
                            p = coef(summary(nadd_gsw_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_gsw_pft_results <- nadd_gsw_photo %>% 
  rbind(nadd_gsw_myc) %>% 
  rbind(nadd_gsw_nfix)

##############
# P addition
##############
padd_gsw_pft <- rma.mv(logr, 
                       logr_var,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results %>% 
                         filter(manip_type == "p" & 
                                  myvar == "gsw" & 
                                  !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_gsw_photo <- data.frame(trait = "gsw", 
                             nut_add = "p",
                             mod = "photo",
                             mod_results(padd_gsw_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_gsw_pft))[2,3],
                             p = coef(summary(padd_gsw_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_gsw_myc <- data.frame(trait = "gsw", 
                           nut_add = "p",
                           mod = "myc_nas",
                           mod_results(padd_gsw_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(padd_gsw_pft))[3,3],
                           p = coef(summary(padd_gsw_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_gsw_nfix <- data.frame(trait = "gsw", 
                            nut_add = "p",
                            mod = "nfix",
                            mod_results(padd_gsw_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_gsw_pft))[4,3],
                            p = coef(summary(padd_gsw_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_gsw_pft_results <- padd_gsw_photo %>% 
  rbind(padd_gsw_myc) %>% 
  rbind(padd_gsw_nfix)

################
# N+P addition
################
npadd_gsw_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "np" & 
                                   myvar == "gsw" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_gsw_photo <- data.frame(trait = "gsw", 
                              nut_add = "np",
                              mod = "photo",
                              mod_results(npadd_gsw_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_gsw_pft))[2,3],
                              p = coef(summary(npadd_gsw_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_gsw_myc <- data.frame(trait = "gsw", 
                            nut_add = "np",
                            mod = "myc_nas",
                            mod_results(npadd_gsw_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(npadd_gsw_pft))[3,3],
                            p = coef(summary(npadd_gsw_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_gsw_nfix <- data.frame(trait = "gsw", 
                             nut_add = "np",
                             mod = "nfix",
                             mod_results(npadd_gsw_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_gsw_pft))[4,3],
                             p = coef(summary(npadd_gsw_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_gsw_pft_results <- npadd_gsw_photo %>% 
  rbind(npadd_gsw_myc) %>% 
  rbind(npadd_gsw_nfix)

#############
# Merge gsw moderator results, with some light cleaning
#############
gsw_pft_summary <- rbind(nadd_gsw_pft_results, 
                         padd_gsw_pft_results, 
                         npadd_gsw_pft_results) %>%
  mutate(k = 47,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Rd - plant functional type
##############################################################################

##############
# N addition
##############
nadd_rd_pft <- rma.mv(logr, 
                      logr_var,
                      method = "REML", 
                      random = ~ 1 | exp, 
                      mods = ~ photo_path + myc_nas + n_fixer,
                      slab = exp, 
                      control = list(stepadj = 0.3), 
                      data = meta_results %>% 
                        filter(manip_type == "n" & 
                                 myvar == "rd" & 
                                 !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_rd_photo <- data.frame(trait = "rd", 
                            nut_add = "n",
                            mod = "photo",
                            mod_results(nadd_rd_pft, 
                                        mod = "photo_path", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_rd_pft))[2,3],
                            p = coef(summary(nadd_rd_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_rd_myc <- data.frame(trait = "rd", 
                          nut_add = "n",
                          mod = "myc_nas",
                          mod_results(nadd_rd_pft, 
                                      mod = "myc_nas", 
                                      group = "exp")$mod_table,
                          z = coef(summary(nadd_rd_pft))[3,3],
                          p = coef(summary(nadd_rd_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_rd_nfix <- data.frame(trait = "rd", 
                           nut_add = "n",
                           mod = "nfix",
                           mod_results(nadd_rd_pft, 
                                       mod = "n_fixer", 
                                       group = "exp")$mod_table,
                           z = coef(summary(nadd_rd_pft))[4,3],
                           p = coef(summary(nadd_rd_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_rd_pft_results <- nadd_rd_photo %>% 
  rbind(nadd_rd_myc) %>% 
  rbind(nadd_rd_nfix)

##############
# P addition
##############
padd_rd_pft <- rma.mv(logr, 
                      logr_var,
                      method = "REML", 
                      random = ~ 1 | exp, 
                      mods = ~ photo_path + myc_nas + n_fixer,
                      slab = exp, 
                      control = list(stepadj = 0.3), 
                      data = meta_results %>% 
                        filter(manip_type == "p" & 
                                 myvar == "rd" & 
                                 !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_rd_photo <- data.frame(trait = "rd", 
                            nut_add = "p",
                            mod = "photo",
                            mod_results(padd_rd_pft, 
                                        mod = "photo_path", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_rd_pft))[2,3],
                            p = coef(summary(padd_rd_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_rd_myc <- data.frame(trait = "rd", 
                          nut_add = "p",
                          mod = "myc_nas",
                          mod_results(padd_rd_pft, 
                                      mod = "myc_nas", 
                                      group = "exp")$mod_table,
                          z = coef(summary(padd_rd_pft))[3,3],
                          p = coef(summary(padd_rd_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_rd_nfix <- data.frame(trait = "rd", 
                           nut_add = "p",
                           mod = "nfix",
                           mod_results(padd_rd_pft, 
                                       mod = "n_fixer", 
                                       group = "exp")$mod_table,
                           z = coef(summary(padd_rd_pft))[4,3],
                           p = coef(summary(padd_rd_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_rd_pft_results <- padd_rd_photo %>% 
  rbind(padd_rd_myc) %>% 
  rbind(padd_rd_nfix)

################
# N+P addition
################
npadd_rd_pft <- rma.mv(logr, 
                       logr_var,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results %>% 
                         filter(manip_type == "np" & 
                                  myvar == "rd" & 
                                  !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_rd_photo <- data.frame(trait = "rd", 
                             nut_add = "np",
                             mod = "photo",
                             mod_results(npadd_rd_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_rd_pft))[2,3],
                             p = coef(summary(npadd_rd_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_rd_myc <- data.frame(trait = "rd", 
                           nut_add = "np",
                           mod = "myc_nas",
                           mod_results(npadd_rd_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(npadd_rd_pft))[3,3],
                           p = coef(summary(npadd_rd_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_rd_nfix <- data.frame(trait = "rd", 
                            nut_add = "np",
                            mod = "nfix",
                            mod_results(npadd_rd_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(npadd_rd_pft))[4,3],
                            p = coef(summary(npadd_rd_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_rd_pft_results <- npadd_rd_photo %>% 
  rbind(npadd_rd_myc) %>% 
  rbind(npadd_rd_nfix)

#############
# Merge Rd moderator results, with some light cleaning
#############
rd_pft_summary <- rbind(nadd_rd_pft_results, 
                        padd_rd_pft_results, 
                        npadd_rd_pft_results) %>%
  mutate(k = 32,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Vcmax - plant functional type
##############################################################################

##############
# N addition
##############
nadd_vcmax_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "vcmax" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_vcmax_photo <- data.frame(trait = "vcmax", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_vcmax_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_vcmax_pft))[2,3],
                               p = coef(summary(nadd_vcmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_vcmax_myc <- data.frame(trait = "vcmax", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_vcmax_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_vcmax_pft))[3,3],
                             p = coef(summary(nadd_vcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_vcmax_nfix <- data.frame(trait = "vcmax", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_vcmax_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_vcmax_pft))[4,3],
                              p = coef(summary(nadd_vcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_vcmax_pft_results <- nadd_vcmax_photo %>% 
  rbind(nadd_vcmax_myc) %>% 
  rbind(nadd_vcmax_nfix)

##############
# P addition
##############
padd_vcmax_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "vcmax" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_vcmax_photo <- data.frame(trait = "vcmax", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_vcmax_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_vcmax_pft))[2,3],
                               p = coef(summary(padd_vcmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_vcmax_myc <- data.frame(trait = "vcmax", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_vcmax_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_vcmax_pft))[3,3],
                             p = coef(summary(padd_vcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_vcmax_nfix <- data.frame(trait = "vcmax", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_vcmax_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_vcmax_pft))[4,3],
                              p = coef(summary(padd_vcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_vcmax_pft_results <- padd_vcmax_photo %>% 
  rbind(padd_vcmax_myc) %>% 
  rbind(padd_vcmax_nfix)

################
# N+P addition
################
npadd_vcmax_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "vcmax" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_vcmax_photo <- data.frame(trait = "vcmax", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_vcmax_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_vcmax_pft))[2,3],
                                p = coef(summary(npadd_vcmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_vcmax_myc <- data.frame(trait = "vcmax", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_vcmax_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_vcmax_pft))[3,3],
                              p = coef(summary(npadd_vcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_vcmax_nfix <- data.frame(trait = "vcmax", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_vcmax_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_vcmax_pft))[4,3],
                               p = coef(summary(npadd_vcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_vcmax_pft_results <- npadd_vcmax_photo %>% 
  rbind(npadd_vcmax_myc) %>% 
  rbind(npadd_vcmax_nfix)

#############
# Merge Vcmax moderator results, with some light cleaning
#############
vcmax_pft_summary <- rbind(nadd_vcmax_pft_results, 
                           padd_vcmax_pft_results, 
                           npadd_vcmax_pft_results) %>%
  mutate(k = 42,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Jmax - plant functional type
##############################################################################

##############
# N addition
##############
nadd_jmax_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "jmax" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_jmax_photo <- data.frame(trait = "jmax", 
                              nut_add = "n",
                              mod = "photo",
                              mod_results(nadd_jmax_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_jmax_pft))[2,3],
                              p = coef(summary(nadd_jmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_jmax_myc <- data.frame(trait = "jmax", 
                            nut_add = "n",
                            mod = "myc_nas",
                            mod_results(nadd_jmax_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_jmax_pft))[3,3],
                            p = coef(summary(nadd_jmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_jmax_nfix <- data.frame(trait = "jmax", 
                             nut_add = "n",
                             mod = "nfix",
                             mod_results(nadd_jmax_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_jmax_pft))[4,3],
                             p = coef(summary(nadd_jmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_jmax_pft_results <- nadd_jmax_photo %>% 
  rbind(nadd_jmax_myc) %>% 
  rbind(nadd_jmax_nfix)

##############
# P addition
##############
padd_jmax_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "jmax" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_jmax_photo <- data.frame(trait = "jmax", 
                              nut_add = "p",
                              mod = "photo",
                              mod_results(padd_jmax_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_jmax_pft))[2,3],
                              p = coef(summary(padd_jmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_jmax_myc <- data.frame(trait = "jmax", 
                            nut_add = "p",
                            mod = "myc_nas",
                            mod_results(padd_jmax_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_jmax_pft))[3,3],
                            p = coef(summary(padd_jmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_jmax_nfix <- data.frame(trait = "jmax", 
                             nut_add = "p",
                             mod = "nfix",
                             mod_results(padd_jmax_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_jmax_pft))[4,3],
                             p = coef(summary(padd_jmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_jmax_pft_results <- padd_jmax_photo %>% 
  rbind(padd_jmax_myc) %>% 
  rbind(padd_jmax_nfix)

################
# N+P addition
################
npadd_jmax_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "jmax" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_jmax_photo <- data.frame(trait = "jmax", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_jmax_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_jmax_pft))[2,3],
                               p = coef(summary(npadd_jmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_jmax_myc <- data.frame(trait = "jmax", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_jmax_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_jmax_pft))[3,3],
                             p = coef(summary(npadd_jmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_jmax_nfix <- data.frame(trait = "jmax", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_jmax_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_jmax_pft))[4,3],
                              p = coef(summary(npadd_jmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_jmax_pft_results <- npadd_jmax_photo %>% 
  rbind(npadd_jmax_myc) %>% 
  rbind(npadd_jmax_nfix)

#############
# Merge Jmax moderator results, with some light cleaning
#############
jmax_pft_summary <- rbind(nadd_jmax_pft_results, 
                          padd_jmax_pft_results, 
                          npadd_jmax_pft_results) %>%
  mutate(k = 40,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Jmax:Vcmax - plant functional type
##############################################################################

##############
# N addition
##############
nadd_jmaxvcmax_pft <- rma.mv(logr, 
                             logr_var,
                             method = "REML", 
                             random = ~ 1 | exp, 
                             mods = ~ photo_path + myc_nas + n_fixer,
                             slab = exp, 
                             control = list(stepadj = 0.3), 
                             data = meta_results %>% 
                               filter(manip_type == "n" & 
                                        myvar == "jmax_vcmax" & 
                                        !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_jmaxvcmax_photo <- data.frame(trait = "jmax_vcmax", 
                                   nut_add = "n",
                                   mod = "photo",
                                   mod_results(nadd_jmaxvcmax_pft, 
                                               mod = "photo_path", 
                                               group = "exp")$mod_table,
                                   z = coef(summary(nadd_jmaxvcmax_pft))[2,3],
                                   p = coef(summary(nadd_jmaxvcmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_jmaxvcmax_myc <- data.frame(trait = "jmax_vcmax", 
                                 nut_add = "n",
                                 mod = "myc_nas",
                                 mod_results(nadd_jmaxvcmax_pft, 
                                             mod = "myc_nas", 
                                             group = "exp")$mod_table,
                                 z = coef(summary(nadd_jmaxvcmax_pft))[3,3],
                                 p = coef(summary(nadd_jmaxvcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_jmaxvcmax_nfix <- data.frame(trait = "jmax_vcmax", 
                                  nut_add = "n",
                                  mod = "nfix",
                                  mod_results(nadd_jmaxvcmax_pft, 
                                              mod = "n_fixer", 
                                              group = "exp")$mod_table,
                                  z = coef(summary(nadd_jmaxvcmax_pft))[4,3],
                                  p = coef(summary(nadd_jmaxvcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_jmaxvcmax_pft_results <- nadd_jmaxvcmax_photo %>% 
  rbind(nadd_jmaxvcmax_myc) %>% 
  rbind(nadd_jmaxvcmax_nfix)

##############
# P addition
##############
padd_jmaxvcmax_pft <- rma.mv(logr, 
                             logr_var,
                             method = "REML", 
                             random = ~ 1 | exp, 
                             mods = ~ photo_path + myc_nas + n_fixer,
                             slab = exp, 
                             control = list(stepadj = 0.3), 
                             data = meta_results %>% 
                               filter(manip_type == "p" & 
                                        myvar == "jmax_vcmax" & 
                                        !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_jmaxvcmax_photo <- data.frame(trait = "jmax_vcmax", 
                                   nut_add = "p",
                                   mod = "photo",
                                   mod_results(padd_jmaxvcmax_pft, 
                                               mod = "photo_path", 
                                               group = "exp")$mod_table,
                                   z = coef(summary(padd_jmaxvcmax_pft))[2,3],
                                   p = coef(summary(padd_jmaxvcmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_jmaxvcmax_myc <- data.frame(trait = "jmax_vcmax", 
                                 nut_add = "p",
                                 mod = "myc_nas",
                                 mod_results(padd_jmaxvcmax_pft, 
                                             mod = "myc_nas", 
                                             group = "exp")$mod_table,
                                 z = coef(summary(padd_jmaxvcmax_pft))[3,3],
                                 p = coef(summary(padd_jmaxvcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_jmaxvcmax_nfix <- data.frame(trait = "jmax_vcmax", 
                                  nut_add = "p",
                                  mod = "nfix",
                                  mod_results(padd_jmaxvcmax_pft, 
                                              mod = "n_fixer", 
                                              group = "exp")$mod_table,
                                  z = coef(summary(padd_jmaxvcmax_pft))[4,3],
                                  p = coef(summary(padd_jmaxvcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_jmaxvcmax_pft_results <- padd_jmaxvcmax_photo %>% 
  rbind(padd_jmaxvcmax_myc) %>% 
  rbind(padd_jmaxvcmax_nfix)

################
# N+P addition
################
npadd_jmaxvcmax_pft <- rma.mv(logr, 
                              logr_var,
                              method = "REML", 
                              random = ~ 1 | exp, 
                              mods = ~ photo_path + myc_nas + n_fixer,
                              slab = exp, 
                              control = list(stepadj = 0.3), 
                              data = meta_results %>% 
                                filter(manip_type == "np" & 
                                         myvar == "jmax_vcmax" & 
                                         !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_jmaxvcmax_photo <- data.frame(trait = "jmax_vcmax", 
                                    nut_add = "np",
                                    mod = "photo",
                                    mod_results(npadd_jmaxvcmax_pft, 
                                                mod = "photo_path", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(npadd_jmaxvcmax_pft))[2,3],
                                    p = coef(summary(npadd_jmaxvcmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_jmaxvcmax_myc <- data.frame(trait = "jmax_vcmax", 
                                  nut_add = "np",
                                  mod = "myc_nas",
                                  mod_results(npadd_jmaxvcmax_pft, 
                                              mod = "myc_nas", 
                                              group = "exp")$mod_table,
                                  z = coef(summary(npadd_jmaxvcmax_pft))[3,3],
                                  p = coef(summary(npadd_jmaxvcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_jmaxvcmax_nfix <- data.frame(trait = "jmax_vcmax", 
                                   nut_add = "np",
                                   mod = "nfix",
                                   mod_results(npadd_jmaxvcmax_pft, 
                                               mod = "n_fixer", 
                                               group = "exp")$mod_table,
                                   z = coef(summary(npadd_jmaxvcmax_pft))[4,3],
                                   p = coef(summary(npadd_jmaxvcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_jmaxvcmax_pft_results <- npadd_jmaxvcmax_photo %>% 
  rbind(npadd_jmaxvcmax_myc) %>% 
  rbind(npadd_jmaxvcmax_nfix)

#############
# Merge Jmax:Vcmax moderator results, with some light cleaning
#############
jmaxvcmax_pft_summary <- rbind(nadd_jmaxvcmax_pft_results, 
                               padd_jmaxvcmax_pft_results, 
                               npadd_jmaxvcmax_pft_results) %>%
  mutate(k = 32,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# PNUE - plant functional type
##############################################################################

##############
# N addition
##############
nadd_pnue_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "leaf_pnue" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_pnue_photo <- data.frame(trait = "pnue", 
                              nut_add = "n",
                              mod = "photo",
                              mod_results(nadd_pnue_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_pnue_pft))[2,3],
                              p = coef(summary(nadd_pnue_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_pnue_myc <- data.frame(trait = "pnue", 
                            nut_add = "n",
                            mod = "myc_nas",
                            mod_results(nadd_pnue_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_pnue_pft))[3,3],
                            p = coef(summary(nadd_pnue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_pnue_nfix <- data.frame(trait = "pnue", 
                             nut_add = "n",
                             mod = "nfix",
                             mod_results(nadd_pnue_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_pnue_pft))[4,3],
                             p = coef(summary(nadd_pnue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_pnue_pft_results <- nadd_pnue_photo %>% 
  rbind(nadd_pnue_myc) %>% 
  rbind(nadd_pnue_nfix)

##############
# P addition
##############
padd_pnue_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "leaf_pnue" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_pnue_photo <- data.frame(trait = "pnue", 
                              nut_add = "p",
                              mod = "photo",
                              mod_results(padd_pnue_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_pnue_pft))[2,3],
                              p = coef(summary(padd_pnue_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_pnue_myc <- data.frame(trait = "pnue", 
                            nut_add = "p",
                            mod = "myc_nas",
                            mod_results(padd_pnue_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_pnue_pft))[3,3],
                            p = coef(summary(padd_pnue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_pnue_nfix <- data.frame(trait = "pnue", 
                             nut_add = "p",
                             mod = "nfix",
                             mod_results(padd_pnue_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_pnue_pft))[4,3],
                             p = coef(summary(padd_pnue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_pnue_pft_results <- padd_pnue_photo %>% 
  rbind(padd_pnue_myc) %>% 
  rbind(padd_pnue_nfix)

################
# N+P addition
################
npadd_pnue_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "leaf_pnue" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_pnue_photo <- data.frame(trait = "pnue", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_pnue_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_pnue_pft))[2,3],
                               p = coef(summary(npadd_pnue_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_pnue_myc <- data.frame(trait = "pnue", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_pnue_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_pnue_pft))[3,3],
                             p = coef(summary(npadd_pnue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_pnue_nfix <- data.frame(trait = "pnue", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_pnue_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_pnue_pft))[4,3],
                              p = coef(summary(npadd_pnue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_pnue_pft_results <- npadd_pnue_photo %>% 
  rbind(npadd_pnue_myc) %>% 
  rbind(npadd_pnue_nfix)

#############
# Merge PNUE moderator results, with some light cleaning
#############
pnue_pft_summary <- rbind(nadd_pnue_pft_results, 
                          padd_pnue_pft_results, 
                          npadd_pnue_pft_results) %>%
  mutate(k = 65,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# PPUE - plant functional type
##############################################################################

##############
# N addition
##############
nadd_ppue_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "leaf_ppue" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_ppue_photo <- data.frame(trait = "ppue", 
                              nut_add = "n",
                              mod = "photo",
                              mod_results(nadd_ppue_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_ppue_pft))[2,3],
                              p = coef(summary(nadd_ppue_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_ppue_myc <- data.frame(trait = "ppue", 
                            nut_add = "n",
                            mod = "myc_nas",
                            mod_results(nadd_ppue_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_ppue_pft))[3,3],
                            p = coef(summary(nadd_ppue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_ppue_nfix <- data.frame(trait = "ppue", 
                             nut_add = "n",
                             mod = "nfix",
                             mod_results(nadd_ppue_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_ppue_pft))[4,3],
                             p = coef(summary(nadd_ppue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_ppue_pft_results <- nadd_ppue_photo %>% 
  rbind(nadd_ppue_myc) %>% 
  rbind(nadd_ppue_nfix)

##############
# P addition
##############
padd_ppue_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "leaf_ppue" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_ppue_photo <- data.frame(trait = "ppue", 
                              nut_add = "p",
                              mod = "photo",
                              mod_results(padd_pnue_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_ppue_pft))[2,3],
                              p = coef(summary(padd_ppue_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_ppue_myc <- data.frame(trait = "ppue", 
                            nut_add = "p",
                            mod = "myc_nas",
                            mod_results(padd_ppue_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_ppue_pft))[3,3],
                            p = coef(summary(padd_ppue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_ppue_nfix <- data.frame(trait = "ppue", 
                             nut_add = "p",
                             mod = "nfix",
                             mod_results(padd_ppue_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_ppue_pft))[4,3],
                             p = coef(summary(padd_ppue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_ppue_pft_results <- padd_ppue_photo %>% 
  rbind(padd_ppue_myc) %>% 
  rbind(padd_ppue_nfix)

################
# N+P addition
################
npadd_ppue_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "leaf_ppue" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_ppue_photo <- data.frame(trait = "ppue", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_ppue_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_ppue_pft))[2,3],
                               p = coef(summary(npadd_ppue_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_ppue_myc <- data.frame(trait = "ppue", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_ppue_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_ppue_pft))[3,3],
                             p = coef(summary(npadd_ppue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_ppue_nfix <- data.frame(trait = "ppue", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_ppue_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_ppue_pft))[4,3],
                              p = coef(summary(npadd_ppue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_ppue_pft_results <- npadd_ppue_photo %>% 
  rbind(npadd_ppue_myc) %>% 
  rbind(npadd_ppue_nfix)

#############
# Merge PPUE moderator results, with some light cleaning
#############
ppue_pft_summary <- rbind(nadd_ppue_pft_results, 
                          padd_ppue_pft_results, 
                          npadd_ppue_pft_results) %>%
  mutate(k = 66,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)











##############################################################################
# Merge myc moderator results and write to .csv
##############################################################################
int_marea_myc_summary %>%
  full_join(int_nmass_myc_summary) %>%
  full_join(int_narea_myc_summary) %>%
  full_join(int_pmass_myc_summary) %>%
  full_join(int_parea_myc_summary) %>%
  full_join(int_leafnp_myc_summary) %>%
  full_join(int_asat_myc_summary) %>%
  full_join(int_vcmax_myc_summary) %>%
  full_join(int_jmax_myc_summary) %>%
  full_join(int_pnue_myc_summary) %>%
  full_join(int_ppue_myc_summary) %>%
  mutate(estimate_perc = (exp(estimate) - 1) * 100,
         lowerCL_perc = (exp(lowerCL) - 1) * 100,
         upperCL_perc = (exp(upperCL) - 1) * 100,
         across(estimate:upperCL_perc, ~round(.x, 3)),
         ci.range = str_c("[", lowerCL, ", ", upperCL, "]"),
         ci.range_perc = str_c("[", lowerCL_perc, ", ", upperCL_perc, "]")) %>%
  write.csv("../data/CNPmeta_myc_moderators_int.csv", row.names = F)





