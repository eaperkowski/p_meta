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
  filter(myvar == "lma" & !is.na(gs_mat)) %>%
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
                           filter(myvar == "lma" & 
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
  filter(myvar == "leaf_n_mass" & !is.na(gs_mat)) %>%
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
                           filter(myvar == "leaf_n_mass" & 
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
  filter(myvar == "leaf_n_area" & !is.na(gs_mat)) %>%
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
                           filter(myvar == "leaf_n_area" & 
                                    !is.na(gs_mat) & dNPi < 2))

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
  filter(myvar == "leaf_p_mass" & !is.na(gs_mat)) %>%
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
                           filter(myvar == "leaf_p_mass" & 
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
  filter(myvar == "leaf_p_area" & !is.na(gs_mat)) %>%
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
                           filter(myvar == "leaf_p_area" & 
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
  filter(myvar == "leaf_np" & !is.na(gs_mat)) %>%
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
                               filter(myvar == "leaf_np" & 
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
  filter(myvar == "total_biomass" & !is.na(gs_mat)) %>%
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
                           filter(myvar == "total_biomass" & 
                                    !is.na(gs_mat)))

# Summary
int_tbio_clim_summary <- data.frame(trait = "tbio",
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
  filter(myvar == "agb" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_agb_clim <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "agb" & 
                                   !is.na(gs_mat) & dNPi < 2))

# Summary
int_agb_clim_summary <- data.frame(trait = "agb",
                                    nut_add = "int",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(int_agb_clim)),
                                    row.names = NULL)

##############################################################################
# Belowground biomass climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(myvar == "bgb" & !is.na(gs_mat)) %>%
  ggplot(aes(x = gs_mat, y = dNPi)) +
  geom_point()

# Model
int_bgb_clim <- rma.mv(yi = dNPi,
                       V = vNPi,
                       W = wNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ gs_mat + gs_ai + gs_par,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(myvar == "bgb" & 
                                  !is.na(gs_mat)))

# Summary
int_bgb_clim_summary <- data.frame(trait = "bgb",
                                   nut_add = "int",
                                   mod = c("intrcpt", "gs_mat",
                                           "gs_ai", "gs_par"),
                                   coef(summary(int_bgb_clim)),
                                   row.names = NULL)

##############################################################################
# Root mass fraction climate moderators
##############################################################################

# Visualize responses
meta_results_int %>% 
  filter(myvar == "rmf" & !is.na(gs_mat)) %>%
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
                         filter(myvar == "rmf" & 
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
  filter(myvar == "rootshoot" & !is.na(gs_mat)) %>%
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
                               filter(myvar == "rootshoot" & 
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
  full_join(int_agb_clim_summary) %>%
  full_join(int_bgb_clim_summary) %>%
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
int_marea_photo <- rma.mv(yi = dNPi,
                          V = vNPi,
                          W = wNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "lma" & 
                                     !is.na(photo_path)))

# Model summary
int_marea_photo_summary <- data.frame(trait = "marea", 
                                       nut_add = "int",
                                       mod_results(int_marea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(int_marea_photo))[2,3],
                                       p = coef(summary(int_marea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Nmass photosynthetic pathway
##############################################################################

# Model
int_nmass_photo <- rma.mv(yi = dNPi,
                          V = vNPi,
                          W = wNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "leaf_n_mass" & 
                                     !is.na(photo_path)))

# Model summary
int_nmass_photo_summary <- data.frame(trait = "nmass", 
                                      nut_add = "int",
                                      mod_results(int_nmass_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(int_nmass_photo))[2,3],
                                      p = coef(summary(int_nmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Narea photosynthetic pathway
##############################################################################

# Model
int_narea_photo <- rma.mv(yi = dNPi,
                          V = vNPi,
                          W = wNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "leaf_n_area" & 
                                     !is.na(photo_path)))

# Model summary
int_narea_photo_summary <- data.frame(trait = "narea", 
                                      nut_add = "int",
                                      mod_results(int_narea_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(int_narea_photo))[2,3],
                                      p = coef(summary(int_narea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Pmass photosynthetic pathway
##############################################################################

# Model
int_pmass_photo <- rma.mv(yi = dNPi,
                          V = vNPi,
                          W = wNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "leaf_p_mass" & 
                                     !is.na(photo_path)))

# Model summary
int_pmass_photo_summary <- data.frame(trait = "pmass", 
                                      nut_add = "int",
                                      mod_results(int_pmass_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(int_pmass_photo))[2,3],
                                      p = coef(summary(int_pmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Parea photosynthetic pathway
##############################################################################

# Model
int_parea_photo <- rma.mv(yi = dNPi,
                          V = vNPi,
                          W = wNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "leaf_p_area" & 
                                     !is.na(photo_path)))

# Model summary
int_parea_photo_summary <- data.frame(trait = "parea", 
                                      nut_add = "int",
                                      mod_results(int_parea_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(int_parea_photo))[2,3],
                                      p = coef(summary(int_parea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)


##############################################################################
# Leaf N:P photosynthetic pathway
##############################################################################

# Model
int_leafnp_photo <- rma.mv(yi = dNPi,
                          V = vNPi,
                          W = wNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "leaf_np" & 
                                     !is.na(photo_path)))


# Model summary
int_leafnp_photo_summary <- data.frame(trait = "leaf_np", 
                                      nut_add = "int",
                                      mod_results(int_leafnp_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(int_leafnp_photo))[2,3],
                                      p = coef(summary(int_leafnp_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)


##############################################################################
# Asat photosynthetic pathway
##############################################################################

# Model
int_asat_photo <- rma.mv(yi = dNPi,
                          V = vNPi,
                          W = wNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "asat" & 
                                     !is.na(photo_path)))

# Model summary
int_asat_photo_summary <- data.frame(trait = "asat", 
                                      nut_add = "int",
                                      mod_results(int_asat_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(int_asat_photo))[2,3],
                                      p = coef(summary(int_asat_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Vcmax photosynthetic pathway
##############################################################################

# Model
int_vcmax_photo <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "vcmax" & 
                                    !is.na(photo_path)))

# Model summary
int_vcmax_photo_summary <- data.frame(trait = "vcmax", 
                                     nut_add = "int",
                                     mod_results(int_vcmax_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_vcmax_photo))[2,3],
                                     p = coef(summary(int_vcmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Jmax photosynthetic pathway
##############################################################################

# Model
int_jmax_photo <- rma.mv(yi = dNPi,
                          V = vNPi,
                          W = wNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "jmax" & 
                                     !is.na(photo_path)))

# Model summary
int_jmax_photo_summary <- data.frame(trait = "jmax", 
                                      nut_add = "int",
                                      mod_results(int_jmax_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(int_jmax_photo))[2,3],
                                      p = coef(summary(int_jmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# PNUE photosynthetic pathway
##############################################################################

# Model
int_pnue_photo <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "leaf_pnue" & 
                                    !is.na(photo_path)))

# Model summary
int_pnue_photo_summary <- data.frame(trait = "pnue", 
                                     nut_add = "int",
                                     mod_results(int_pnue_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_pnue_photo))[2,3],
                                     p = coef(summary(int_pnue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# PPUE photosynthetic pathway
##############################################################################

# Model
int_ppue_photo <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "leaf_ppue" & 
                                    !is.na(photo_path)))

# Model summary
int_ppue_photo_summary <- data.frame(trait = "ppue", 
                                     nut_add = "int",
                                     mod_results(int_ppue_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_ppue_photo))[2,3],
                                     p = coef(summary(int_ppue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Total biomass photosynthetic pathway
##############################################################################

# Model
int_tbio_photo <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "total_biomass" & 
                                    !is.na(photo_path)))

# Model summary
int_tbio_photo_summary <- data.frame(trait = "total_biomass", 
                                     nut_add = "int",
                                     mod_results(int_tbio_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_tbio_photo))[2,3],
                                     p = coef(summary(int_tbio_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Aboveground biomass photosynthetic pathway
##############################################################################

# Model
int_agb_photo <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "agb" & 
                                    !is.na(photo_path)))

# Model summary
int_agb_photo_summary <- data.frame(trait = "agb", 
                                     nut_add = "int",
                                     mod_results(int_agb_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_agb_photo))[2,3],
                                     p = coef(summary(int_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)


##############################################################################
# Belowground biomass photosynthetic pathway
##############################################################################

# Model
int_bgb_photo <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "bgb" & 
                                   !is.na(photo_path)))

# Model summary
int_bgb_photo_summary <- data.frame(trait = "bgb", 
                                    nut_add = "int",
                                    mod_results(int_bgb_photo, 
                                                mod = "photo_path", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_bgb_photo))[2,3],
                                    p = coef(summary(int_bgb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Merge photo moderator results and write to .csv
##############################################################################
int_marea_photo_summary %>%
  full_join(int_nmass_photo_summary) %>%
  full_join(int_narea_photo_summary) %>%
  full_join(int_pmass_photo_summary) %>%
  full_join(int_parea_photo_summary) %>%
  full_join(int_leafnp_photo_summary) %>%
  full_join(int_asat_photo_summary) %>%
  full_join(int_vcmax_photo_summary) %>%
  full_join(int_jmax_photo_summary) %>%
  full_join(int_pnue_photo_summary) %>%
  full_join(int_ppue_photo_summary) %>%
  full_join(int_tbio_photo_summary) %>%
  full_join(int_agb_photo_summary) %>%
  full_join(int_bgb_photo_summary) %>%
  mutate(estimate_perc = (exp(estimate) - 1) * 100,
         lowerCL_perc = (exp(lowerCL) - 1) * 100,
         upperCL_perc = (exp(upperCL) - 1) * 100,
         across(estimate:upperCL_perc, ~round(.x, 3)),
         ci.range = str_c("[", lowerCL, ", ", upperCL, "]"),
         ci.range_perc = str_c("[", lowerCL_perc, ", ", upperCL_perc, "]")) %>%
  write.csv("../data/CNPmeta_photo_moderators_int.csv", row.names = F)

##############################################################################
# Marea Nfixation
##############################################################################

# Model
int_marea_nfix <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "lma" & 
                                    !is.na(n_fixer)))

# Model summary
int_marea_nfix_summary <- data.frame(trait = "marea", 
                                    nut_add = "int",
                                    mod_results(int_marea_nfix, 
                                                mod = "n_fixer", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_marea_nfix))[2,3],
                                    p = coef(summary(int_marea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Nmass Nfixation
##############################################################################

# Model
int_nmass_nfix <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "leaf_n_mass" & 
                                    !is.na(n_fixer)))

# Model summary
int_nmass_nfix_summary <- data.frame(trait = "nmass", 
                                     nut_add = "int",
                                     mod_results(int_nmass_nfix, 
                                                 mod = "n_fixer", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_nmass_nfix))[2,3],
                                     p = coef(summary(int_nmass_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Narea Nfixation
##############################################################################

# Model
int_narea_nfix <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "leaf_n_area" & 
                                    !is.na(n_fixer)))

# Model summary
int_narea_nfix_summary <- data.frame(trait = "narea", 
                                     nut_add = "int",
                                     mod_results(int_narea_nfix, 
                                                 mod = "n_fixer", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_narea_nfix))[2,3],
                                     p = coef(summary(int_narea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Pmass Nfixation
##############################################################################

# Model
int_pmass_nfix <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "leaf_p_mass" & 
                                    !is.na(n_fixer)))

# Model summary
int_pmass_nfix_summary <- data.frame(trait = "pmass", 
                                     nut_add = "int",
                                     mod_results(int_pmass_nfix, 
                                                 mod = "n_fixer", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_pmass_nfix))[2,3],
                                     p = coef(summary(int_pmass_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Parea Nfixation
##############################################################################

# Model
int_parea_nfix <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "leaf_p_area" & 
                                    !is.na(n_fixer)))

# Model summary
int_parea_nfix_summary <- data.frame(trait = "parea", 
                                     nut_add = "int",
                                     mod_results(int_parea_nfix, 
                                                 mod = "n_fixer", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_parea_nfix))[2,3],
                                     p = coef(summary(int_parea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Leaf N:P Nfixation
##############################################################################

# Model
int_leafnp_nfix <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_np" & 
                                   !is.na(n_fixer)))

# Model summary
int_leafnp_nfix_summary <- data.frame(trait = "leaf_np", 
                                    nut_add = "int",
                                    mod_results(int_leafnp_nfix, 
                                                mod = "n_fixer", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_leafnp_nfix))[2,3],
                                    p = coef(summary(int_leafnp_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Asat Nfixation
##############################################################################

# Model
int_asat_nfix <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "asat" & 
                                    !is.na(n_fixer)))

# Model summary
int_asat_nfix_summary <- data.frame(trait = "asat", 
                                     nut_add = "int",
                                     mod_results(int_asat_nfix, 
                                                 mod = "n_fixer", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_asat_nfix))[2,3],
                                     p = coef(summary(int_asat_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Vcmax Nfixation
##############################################################################

# Model
int_vcmax_nfix <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "vcmax" & 
                                    !is.na(n_fixer)))

# Model summary
int_vcmax_nfix_summary <- data.frame(trait = "vcmax", 
                                     nut_add = "int",
                                     mod_results(int_vcmax_nfix, 
                                                 mod = "n_fixer", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_vcmax_nfix))[2,3],
                                     p = coef(summary(int_vcmax_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Jmax Nfixation
##############################################################################

# Model
int_jmax_nfix <- rma.mv(yi = dNPi,
                         V = vNPi,
                         W = wNPi,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results_int %>% 
                           filter(myvar == "jmax" & 
                                    !is.na(n_fixer)))

# Model summary
int_jmax_nfix_summary <- data.frame(trait = "jmax", 
                                     nut_add = "int",
                                     mod_results(int_jmax_nfix, 
                                                 mod = "n_fixer", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(int_jmax_nfix))[2,3],
                                     p = coef(summary(int_jmax_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# PNUE Nfixation
##############################################################################

# Model
int_pnue_nfix <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_pnue" & 
                                   !is.na(n_fixer)))

# Model summary
int_pnue_nfix_summary <- data.frame(trait = "pnue", 
                                    nut_add = "int",
                                    mod_results(int_pnue_nfix, 
                                                mod = "n_fixer", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_pnue_nfix))[2,3],
                                    p = coef(summary(int_pnue_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# PPUE Nfixation
##############################################################################

# Model
int_ppue_nfix <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_ppue" & 
                                   !is.na(n_fixer)))

# Model summary
int_ppue_nfix_summary <- data.frame(trait = "ppue", 
                                    nut_add = "int",
                                    mod_results(int_ppue_nfix, 
                                                mod = "n_fixer", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_ppue_nfix))[2,3],
                                    p = coef(summary(int_ppue_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Merge nfix moderator results and write to .csv
##############################################################################
int_marea_nfix_summary %>%
  full_join(int_nmass_nfix_summary) %>%
  full_join(int_narea_nfix_summary) %>%
  full_join(int_pmass_nfix_summary) %>%
  full_join(int_parea_nfix_summary) %>%
  full_join(int_leafnp_nfix_summary) %>%
  full_join(int_asat_nfix_summary) %>%
  full_join(int_vcmax_nfix_summary) %>%
  full_join(int_jmax_nfix_summary) %>%
  full_join(int_pnue_nfix_summary) %>%
  full_join(int_ppue_nfix_summary) %>%
  mutate(estimate_perc = (exp(estimate) - 1) * 100,
         lowerCL_perc = (exp(lowerCL) - 1) * 100,
         upperCL_perc = (exp(upperCL) - 1) * 100,
         across(estimate:upperCL_perc, ~round(.x, 3)),
         ci.range = str_c("[", lowerCL, ", ", upperCL, "]"),
         ci.range_perc = str_c("[", lowerCL_perc, ", ", upperCL_perc, "]"))
  write.csv("../data/CNPmeta_nfix_moderators_int.csv", row.names = F)


##############################################################################
# Marea mycorrhizal acquisition strategy
##############################################################################

# Model
int_marea_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "lma" & 
                                   !is.na(myc_nas)))

# Model summary
int_marea_myc_summary <- data.frame(trait = "marea", 
                                    nut_add = "int",
                                    mod_results(int_marea_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_marea_myc))[2,3],
                                    p = coef(summary(int_marea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Nmass mycorrhizal acquisition strategy
##############################################################################

# Model
int_nmass_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_n_mass" & 
                                   !is.na(myc_nas)))

# Model summary
int_nmass_myc_summary <- data.frame(trait = "nmass", 
                                    nut_add = "int",
                                    mod_results(int_nmass_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_nmass_myc))[2,3],
                                    p = coef(summary(int_nmass_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Narea mycorrhizal acquisition strategy
##############################################################################

# Model
int_narea_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_n_area" & 
                                   !is.na(myc_nas)))

# Model summary
int_narea_myc_summary <- data.frame(trait = "narea", 
                                    nut_add = "int",
                                    mod_results(int_narea_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_narea_myc))[2,3],
                                    p = coef(summary(int_narea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Pmass mycorrhizal acquisition strategy
##############################################################################

# Model
int_pmass_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_p_mass" & 
                                   !is.na(myc_nas)))

# Model summary
int_pmass_myc_summary <- data.frame(trait = "pmass", 
                                    nut_add = "int",
                                    mod_results(int_pmass_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_pmass_myc))[2,3],
                                    p = coef(summary(int_pmass_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Parea mycorrhizal acquisition strategy
##############################################################################

# Model
int_parea_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_p_area" & 
                                   !is.na(myc_nas)))

# Model summary
int_parea_myc_summary <- data.frame(trait = "parea", 
                                    nut_add = "int",
                                    mod_results(int_parea_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_parea_myc))[2,3],
                                    p = coef(summary(int_parea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Leaf N:P mycorrhizal acquisition strategy
##############################################################################

# Model
int_leafnp_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_np" & 
                                   !is.na(myc_nas)))

# Model summary
int_leafnp_myc_summary <- data.frame(trait = "leafnp", 
                                    nut_add = "int",
                                    mod_results(int_leafnp_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_leafnp_myc))[2,3],
                                    p = coef(summary(int_leafnp_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Asat mycorrhizal acquisition strategy
##############################################################################

# Model
int_asat_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "asat" & 
                                   !is.na(myc_nas)))

# Model summary
int_asat_myc_summary <- data.frame(trait = "asat", 
                                    nut_add = "int",
                                    mod_results(int_asat_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_asat_myc))[2,3],
                                    p = coef(summary(int_asat_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Vcmax mycorrhizal acquisition strategy
##############################################################################

# Model
int_vcmax_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "vcmax" & 
                                   !is.na(myc_nas)))

# Model summary
int_vcmax_myc_summary <- data.frame(trait = "vcmax", 
                                    nut_add = "int",
                                    mod_results(int_vcmax_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_vcmax_myc))[2,3],
                                    p = coef(summary(int_vcmax_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# Jmax mycorrhizal acquisition strategy
##############################################################################

# Model
int_jmax_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "jmax" & 
                                   !is.na(myc_nas)))

# Model summary
int_jmax_myc_summary <- data.frame(trait = "jmax", 
                                    nut_add = "int",
                                    mod_results(int_jmax_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_jmax_myc))[2,3],
                                    p = coef(summary(int_jmax_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# PNUE mycorrhizal acquisition strategy
##############################################################################

# Model
int_pnue_myc <- rma.mv(yi = dNPi,
                        V = vNPi,
                        W = wNPi,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ myc_nas,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results_int %>% 
                          filter(myvar == "leaf_pnue" & 
                                   !is.na(myc_nas)))

# Model summary
int_pnue_myc_summary <- data.frame(trait = "pnue", 
                                    nut_add = "int",
                                    mod_results(int_pnue_myc, 
                                                mod = "myc_nas", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(int_pnue_myc))[2,3],
                                    p = coef(summary(int_pnue_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

##############################################################################
# PPUE mycorrhizal acquisition strategy
##############################################################################

# Model
int_ppue_myc <- rma.mv(yi = dNPi,
                       V = vNPi,
                       W = wNPi,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ myc_nas,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(myvar == "leaf_ppue" & 
                                  !is.na(myc_nas)))

# Model summary
int_ppue_myc_summary <- data.frame(trait = "ppue", 
                                   nut_add = "int",
                                   mod_results(int_ppue_myc, 
                                               mod = "myc_nas", 
                                               group = "exp")$mod_table,
                                   z = coef(summary(int_ppue_myc))[2,3],
                                   p = coef(summary(int_ppue_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

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





