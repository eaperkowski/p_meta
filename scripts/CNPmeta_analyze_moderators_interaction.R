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
  filter(response == "lma" & !is.na(gs_mat) & gs_ai < 3 & dNPi < 2) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_marea_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "lma" & 
                                    !is.na(gs_mat) & gs_ai < 3 & dNPi < 2))

# Summary
int_marea_clim_summary <- data.frame(trait = "marea",
                                     nut_add = "int",
                                     k = 78,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_marea_clim)),
                                     row.names = NULL)

##############################################################################
# Nmass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_n_mass" & !is.na(gs_mat) & gs_ai < 3 & dNPi > -1.9) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_nmass_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_n_mass" & 
                                    !is.na(gs_mat) & gs_ai < 3 & dNPi > -1.9))

# Summary
int_nmass_clim_summary <- data.frame(trait = "nmass",
                                     nut_add = "int",
                                     k = 101,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_nmass_clim)),
                                     row.names = NULL)

##############################################################################
# Narea climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_n_area" & !is.na(gs_mat) & gs_ai < 3) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_narea_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_n_area" & 
                                    !is.na(gs_mat) & gs_ai < 3))

# Summary
int_narea_clim_summary <- data.frame(trait = "narea",
                                     nut_add = "int",
                                     k = 50,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_narea_clim)),
                                     row.names = NULL)

##############################################################################
# Pmass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_p_mass" & !is.na(gs_mat) & gs_ai < 3) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_pmass_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_p_mass" & 
                                    !is.na(gs_mat) & gs_ai < 3))

# Summary
int_pmass_clim_summary <- data.frame(trait = "pmass",
                                     nut_add = "int",
                                     k = 97,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_pmass_clim)),
                                     row.names = NULL)

##############################################################################
# Parea climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_p_area" & !is.na(gs_mat) & gs_ai < 3 & dNPi < 2) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_parea_clim <- rma.mv(yi = dNPi,
                         V = vNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_p_area" & 
                                    !is.na(gs_mat) & gs_ai < 3 & dNPi < 2))

# Summary
int_parea_clim_summary <- data.frame(trait = "parea",
                                     nut_add = "int",
                                     k = 44,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_parea_clim)),
                                     row.names = NULL)

##############################################################################
# Leaf N:P climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_np" & !is.na(gs_mat) & gs_ai < 3 & dNPi > -2) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_leafnp_clim <- rma.mv(yi = dNPi,
                          V = vNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(response == "leaf_np" & 
                                     !is.na(gs_mat) & gs_ai < 3 & dNPi > -2))

# Summary
int_leafnp_clim_summary <- data.frame(trait = "leaf_np",
                                      nut_add = "int",
                                      k = 82,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(int_leafnp_clim)),
                                      row.names = NULL)

##############################################################################
# Total biomass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "tbio_gm2" & !is.na(gs_mat) & gs_ai < 3 & dNPi < 1.25) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_tbio_clim <- rma.mv(yi = dNPi,
                        V = vNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "tbio_gm2" & 
                                   !is.na(gs_mat) & gs_ai < 3 & dNPi < 1.25))

# Summary
int_tbio_clim_summary <- data.frame(trait = "tbio_gm2",
                                    nut_add = "int",
                                    k = 29,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(int_tbio_clim)),
                                    row.names = NULL)

##############################################################################
# Aboveground biomass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "anpp" & !is.na(gs_mat) & gs_ai < 3 & dNPi < 2 & dNPi > -1) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_anpp_clim <- rma.mv(yi = dNPi,
                        V = vNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "anpp" & 
                                   !is.na(gs_mat) & gs_ai < 3 & dNPi < 2 & dNPi > -2))

# Summary
int_anpp_clim_summary <- data.frame(trait = "anpp",
                                    nut_add = "int",
                                    k = 104,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(int_anpp_clim)),
                                    row.names = NULL)

##############################################################################
# Belowground biomass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "bnpp" & !is.na(gs_mat) & gs_ai < 3 & dNPi > -1) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_bnpp_clim <- rma.mv(yi = dNPi,
                        V = vNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "bnpp" & 
                                   !is.na(gs_mat) & gs_ai < 3 & dNPi > -1))

# Summary
int_bnpp_clim_summary <- data.frame(trait = "bnpp",
                                    nut_add = "int",
                                    k = 51,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(int_bnpp_clim)),
                                    row.names = NULL)

##############################################################################
# Root mass fraction climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "rmf" & !is.na(gs_mat) & gs_ai < 3) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_rmf_clim <- rma.mv(yi = dNPi,
                       V = vNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ gs_mat + gs_ai + gs_par,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(response == "rmf" & 
                                  !is.na(gs_mat) & gs_ai < 3))

# Summary
int_rmf_clim_summary <- data.frame(trait = "rmf",
                                   nut_add = "int",
                                   k = 34,
                                   mod = c("intrcpt", "gs_mat",
                                           "gs_ai", "gs_par"),
                                   coef(summary(int_rmf_clim)),
                                   row.names = NULL)

##############################################################################
# Root:shoot climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "rootshoot" & !is.na(gs_mat) & gs_ai < 3) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_rootshoot_clim <- rma.mv(yi = dNPi,
                             V = vNPi,
                             method = "REML", 
                             random = ~ 1 | exp, 
                             mods = ~ gs_mat + gs_ai + gs_par,
                             slab = exp, 
                             control = list(stepadj = 0.3), 
                             data = meta_results_int %>% 
                               filter(response == "rootshoot" & 
                                        !is.na(gs_mat) & gs_ai < 3))

# Summary
int_rootshoot_clim_summary <- data.frame(trait = "rootshoot",
                                         nut_add = "int",
                                         k = 36,
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
  mutate(across(estimate:ci.ub, ~round(.x, 3)),
         estimate_se = str_c(sprintf("%.3f", estimate), "±", sprintf("%.3f", se)),
         ci.range = str_c("[", sprintf("%.3f", ci.lb), ", ", sprintf("%.3f", ci.ub), "]")) %>%
  dplyr::select(trait:se, estimate_se, zval:ci.ub, ci.range) %>%
  write_excel_csv("../data/CNPmeta_clim_moderators_int.csv")

##############################################################################
# Marea photosynthetic pathway
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "lma" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_marea_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
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

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_n_mass" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_nmass_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
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

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_n_area" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_narea_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
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
  mutate(k = 86,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Pmass - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_p_mass" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_pmass_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
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
  mutate(k = 128,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Parea - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_p_area" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_parea_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "leaf_p_area" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_parea_photo <- data.frame(trait = "parea", 
                              nut_add = "int",
                              mod = "photo",
                              mod_results(int_parea_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(int_parea_pft))[2,3],
                              p = coef(summary(int_parea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_parea_myc <- data.frame(trait = "parea", 
                            nut_add = "int",
                            mod = "myc_nas",
                            mod_results(int_parea_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_parea_pft))[3,3],
                            p = coef(summary(int_parea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_parea_nfix <- data.frame(trait = "parea", 
                             nut_add = "int",
                             mod = "nfix",
                             mod_results(int_parea_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_parea_pft))[4,3],
                             p = coef(summary(int_parea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_parea_pft_results <- int_parea_photo %>% 
  rbind(int_parea_myc) %>% 
  rbind(int_parea_nfix) %>%
  mutate(k = 79,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Leaf N:P - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_np" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_leafnp_pft <- rma.mv(yi = dNPi,
                         V = vNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(response == "leaf_np" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_leafnp_photo <- data.frame(trait = "leaf_np", 
                               nut_add = "int",
                               mod = "photo",
                               mod_results(int_leafnp_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(int_leafnp_pft))[2,3],
                               p = coef(summary(int_leafnp_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_leafnp_myc <- data.frame(trait = "leaf_np", 
                             nut_add = "int",
                             mod = "myc_nas",
                             mod_results(int_leafnp_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_leafnp_pft))[3,3],
                             p = coef(summary(int_leafnp_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_leafnp_nfix <- data.frame(trait = "leaf_np", 
                              nut_add = "int",
                              mod = "nfix",
                              mod_results(int_leafnp_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(int_leafnp_pft))[4,3],
                              p = coef(summary(int_leafnp_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_leafnp_pft_results <- int_leafnp_photo %>% 
  rbind(int_leafnp_myc) %>% 
  rbind(int_leafnp_nfix) %>%
  mutate(k = 115,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Asat - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "asat" & !is.na(photo_path) & dNPi < 5) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_asat_pft <- rma.mv(yi = dNPi,
                       V = vNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(response == "asat" & 
                                  !is.na(photo_path) & dNPi < 5))

# Extract photosynthetic pathway summary statistics
int_asat_photo <- data.frame(trait = "asat", 
                             nut_add = "int",
                             mod = "photo",
                             mod_results(int_asat_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_asat_pft))[2,3],
                             p = coef(summary(int_asat_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_asat_myc <- data.frame(trait = "asat", 
                           nut_add = "int",
                           mod = "myc_nas",
                           mod_results(int_asat_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(int_asat_pft))[3,3],
                           p = coef(summary(int_asat_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_asat_nfix <- data.frame(trait = "asat", 
                            nut_add = "int",
                            mod = "nfix",
                            mod_results(int_asat_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_asat_pft))[4,3],
                            p = coef(summary(int_asat_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_asat_pft_results <- int_asat_photo %>% 
  rbind(int_asat_myc) %>% 
  rbind(int_asat_nfix) %>%
  mutate(k = 90,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# gsw - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "gsw" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_gsw_pft <- rma.mv(yi = dNPi,
                      V = vNPi,
                      method = "REML", 
                      random = ~ 1 | exp, 
                      mods = ~ photo_path + myc_nas + n_fixer,
                      slab = exp, 
                      control = list(stepadj = 0.3), 
                      data = meta_results_int %>% 
                        filter(response == "gsw" & 
                                 !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_gsw_photo <- data.frame(trait = "gsw", 
                            nut_add = "int",
                            mod = "photo",
                            mod_results(int_gsw_pft, 
                                        mod = "photo_path", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_gsw_pft))[2,3],
                            p = coef(summary(int_gsw_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_gsw_myc <- data.frame(trait = "gsw", 
                          nut_add = "int",
                          mod = "myc_nas",
                          mod_results(int_gsw_pft, 
                                      mod = "myc_nas", 
                                      group = "exp")$mod_table,
                          z = coef(summary(int_gsw_pft))[3,3],
                          p = coef(summary(int_gsw_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_gsw_nfix <- data.frame(trait = "gsw", 
                           nut_add = "int",
                           mod = "nfix",
                           mod_results(int_gsw_pft, 
                                       mod = "n_fixer", 
                                       group = "exp")$mod_table,
                           z = coef(summary(int_gsw_pft))[4,3],
                           p = coef(summary(int_gsw_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_gsw_pft_results <- int_gsw_photo %>% 
  rbind(int_gsw_myc) %>% 
  rbind(int_gsw_nfix) %>%
  mutate(k = 47,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Rd - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "rd" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_rd_pft <- rma.mv(yi = dNPi,
                     V = vNPi,
                     method = "REML", 
                     random = ~ 1 | exp, 
                     mods = ~ photo_path + myc_nas + n_fixer,
                     slab = exp, 
                     control = list(stepadj = 0.3), 
                     data = meta_results_int %>% 
                       filter(response == "rd" & 
                                !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_rd_photo <- data.frame(trait = "rd", 
                           nut_add = "int",
                           mod = "photo",
                           mod_results(int_rd_pft, 
                                       mod = "photo_path", 
                                       group = "exp")$mod_table,
                           z = coef(summary(int_rd_pft))[2,3],
                           p = coef(summary(int_rd_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_rd_myc <- data.frame(trait = "rd", 
                         nut_add = "int",
                         mod = "myc_nas",
                         mod_results(int_rd_pft, 
                                     mod = "myc_nas", 
                                     group = "exp")$mod_table,
                         z = coef(summary(int_rd_pft))[3,3],
                         p = coef(summary(int_rd_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_rd_nfix <- data.frame(trait = "rd", 
                          nut_add = "int",
                          mod = "nfix",
                          mod_results(int_rd_pft, 
                                      mod = "n_fixer", 
                                      group = "exp")$mod_table,
                          z = coef(summary(int_rd_pft))[4,3],
                          p = coef(summary(int_rd_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_rd_pft_results <- int_rd_photo %>% 
  rbind(int_rd_myc) %>% 
  rbind(int_rd_nfix) %>%
  mutate(k = 32,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Vcmax - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "vcmax" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_vcmax_pft <- rma.mv(yi = dNPi,
                        V = vNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(response == "vcmax" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_vcmax_photo <- data.frame(trait = "vcmax", 
                              nut_add = "int",
                              mod = "photo",
                              mod_results(int_vcmax_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(int_vcmax_pft))[2,3],
                              p = coef(summary(int_vcmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_vcmax_myc <- data.frame(trait = "vcmax", 
                            nut_add = "int",
                            mod = "myc_nas",
                            mod_results(int_vcmax_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_vcmax_pft))[3,3],
                            p = coef(summary(int_vcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_vcmax_nfix <- data.frame(trait = "vcmax", 
                             nut_add = "int",
                             mod = "nfix",
                             mod_results(int_vcmax_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_vcmax_pft))[4,3],
                             p = coef(summary(int_vcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_vcmax_pft_results <- int_vcmax_photo %>% 
  rbind(int_vcmax_myc) %>% 
  rbind(int_vcmax_nfix) %>%
  mutate(k = 42,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Jmax - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "jmax" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_jmax_pft <- rma.mv(yi = dNPi,
                       V = vNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(response == "jmax" & 
                                  !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_jmax_photo <- data.frame(trait = "jmax", 
                             nut_add = "int",
                             mod = "photo",
                             mod_results(int_jmax_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_jmax_pft))[2,3],
                             p = coef(summary(int_jmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_jmax_myc <- data.frame(trait = "jmax", 
                           nut_add = "int",
                           mod = "myc_nas",
                           mod_results(int_vcmax_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(int_jmax_pft))[3,3],
                           p = coef(summary(int_jmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_jmax_nfix <- data.frame(trait = "jmax", 
                            nut_add = "int",
                            mod = "nfix",
                            mod_results(int_jmax_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_jmax_pft))[4,3],
                            p = coef(summary(int_jmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_jmax_pft_results <- int_jmax_photo %>% 
  rbind(int_jmax_myc) %>% 
  rbind(int_jmax_nfix) %>%
  mutate(k = 39,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Jmax:Vcmax - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "jmax_vcmax" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_jmaxvcmax_pft <- rma.mv(yi = dNPi,
                            V = vNPi,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path + myc_nas + n_fixer,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results_int %>% 
                              filter(response == "jmax_vcmax" & 
                                       !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_jmaxvcmax_photo <- data.frame(trait = "jmax_vcmax", 
                                  nut_add = "int",
                                  mod = "photo",
                                  mod_results(int_jmaxvcmax_pft, 
                                              mod = "photo_path", 
                                              group = "exp")$mod_table,
                                  z = coef(summary(int_jmaxvcmax_pft))[2,3],
                                  p = coef(summary(int_jmaxvcmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_jmaxvcmax_myc <- data.frame(trait = "jmax_vcmax", 
                                nut_add = "int",
                                mod = "myc_nas",
                                mod_results(int_jmaxvcmax_pft, 
                                            mod = "myc_nas", 
                                            group = "exp")$mod_table,
                                z = coef(summary(int_jmaxvcmax_pft))[3,3],
                                p = coef(summary(int_jmaxvcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_jmaxvcmax_nfix <- data.frame(trait = "jmax_vcmax", 
                                 nut_add = "int",
                                 mod = "nfix",
                                 mod_results(int_jmaxvcmax_pft, 
                                             mod = "n_fixer", 
                                             group = "exp")$mod_table,
                                 z = coef(summary(int_jmaxvcmax_pft))[4,3],
                                 p = coef(summary(int_jmaxvcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_jmaxvcmax_pft_results <- int_jmaxvcmax_photo %>% 
  rbind(int_jmaxvcmax_myc) %>% 
  rbind(int_jmaxvcmax_nfix) %>%
  mutate(k = 32,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# PNUE - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_pnue" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_pnue_pft <- rma.mv(yi = dNPi,
                       V = vNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(response == "leaf_pnue" & 
                                  !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_pnue_photo <- data.frame(trait = "pnue", 
                             nut_add = "int",
                             mod = "photo",
                             mod_results(int_pnue_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_pnue_pft))[2,3],
                             p = coef(summary(int_pnue_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_pnue_myc <- data.frame(trait = "pnue", 
                           nut_add = "int",
                           mod = "myc_nas",
                           mod_results(int_pnue_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(int_pnue_pft))[3,3],
                           p = coef(summary(int_pnue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_pnue_nfix <- data.frame(trait = "pnue", 
                            nut_add = "int",
                            mod = "nfix",
                            mod_results(int_pnue_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_pnue_pft))[4,3],
                            p = coef(summary(int_pnue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_pnue_pft_results <- int_pnue_photo %>% 
  rbind(int_pnue_myc) %>% 
  rbind(int_pnue_nfix) %>%
  mutate(k = 65,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# PPUE - plant functional type
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(response == "leaf_ppue" & !is.na(photo_path)) %>%
  ggplot(aes(x = photo_path, y = dNPi)) +
  geom_point()

# Model
int_ppue_pft <- rma.mv(yi = dNPi,
                       V = vNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(response == "leaf_ppue" & 
                                  !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
int_ppue_photo <- data.frame(trait = "ppue", 
                             nut_add = "int",
                             mod = "photo",
                             mod_results(int_ppue_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(int_ppue_pft))[2,3],
                             p = coef(summary(int_ppue_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
int_ppue_myc <- data.frame(trait = "ppue", 
                           nut_add = "int",
                           mod = "myc_nas",
                           mod_results(int_ppue_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(int_ppue_pft))[3,3],
                           p = coef(summary(int_ppue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
int_ppue_nfix <- data.frame(trait = "ppue", 
                            nut_add = "int",
                            mod = "nfix",
                            mod_results(int_ppue_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(int_ppue_pft))[4,3],
                            p = coef(summary(int_ppue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

#############
# Merge summary statistics into single data frame
#############

int_ppue_pft_results <- int_ppue_photo %>% 
  rbind(int_ppue_myc) %>% 
  rbind(int_ppue_nfix) %>%
  mutate(k = 64,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Merge moderator results and write to .csv
##############################################################################
int_marea_pft_results %>%
  full_join(int_nmass_pft_results) %>%
  full_join(int_narea_pft_results) %>%
  full_join(int_pmass_pft_results) %>%
  full_join(int_parea_pft_results) %>%
  full_join(int_leafnp_pft_results) %>%
  full_join(int_asat_pft_results) %>%
  full_join(int_rd_pft_results) %>%
  full_join(int_gsw_pft_results) %>%
  full_join(int_vcmax_pft_results) %>%
  full_join(int_jmax_pft_results) %>%
  full_join(int_jmaxvcmax_pft_results) %>%
  full_join(int_pnue_pft_results) %>%
  full_join(int_ppue_pft_results) %>%
  mutate(across(estimate:upperCL, ~round(as.numeric(.x), 3)),
         ci.range = str_c("[", sprintf("%.3f", lowerCL), ", ", sprintf("%.3f", upperCL), "]")) %>%
  write.csv("../data/CNPmeta_pft_moderators_int.csv", row.names = F)
