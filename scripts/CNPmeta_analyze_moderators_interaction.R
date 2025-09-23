# Libraries
library(tidyverse)
library(metafor)
library(ggpubr)
library(naniar) # to resolve NA/<NA> issue
library(orchaRd)

# Load meta-analysis results
meta_results_int <- read.csv("../data/CNPmeta_logr_results_int.csv")

# Check file structure
head(meta_results_int)

##############################################################################
# Marea climate moderators
##############################################################################

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
                                    !is.na(gs_mat)))

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
                                    !is.na(gs_mat)))

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
int_pmass_clim_summary <- data.frame(trait = "nmass",
                                     nut_add = "int",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_pmass_clim)),
                                     row.names = NULL)

##############################################################################
# Parea climate moderators
##############################################################################

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
                                    !is.na(gs_mat)))

# Summary
int_parea_clim_summary <- data.frame(trait = "parea",
                                     nut_add = "int",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(int_parea_clim)),
                                     row.names = NULL)

##############################################################################
# Total biomass climate moderators
##############################################################################

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
                                   !is.na(gs_mat)))

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

# N addition
nadd_rmf_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "rmf" & 
                                   !is.na(gs_mat)))

nadd_rmf_clim_summary <- data.frame(trait = "rmf",
                                    nut_add = "n",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(nadd_rmf_clim)),
                                    row.names = NULL)

# P addition
padd_rmf_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "rmf" & 
                                   !is.na(gs_mat)))

padd_rmf_clim_summary <- data.frame(trait = "rmf",
                                    nut_add = "p",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(padd_rmf_clim)),
                                    row.names = NULL)

# N+P addition
npadd_rmf_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "rmf" & 
                                    !is.na(gs_mat)))

npadd_rmf_clim_summary <- data.frame(trait = "rmf",
                                     nut_add = "np",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_rmf_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
rmf_clim_summary <- rbind(nadd_rmf_clim_summary, 
                          padd_rmf_clim_summary, 
                          npadd_rmf_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Root:shoot climate moderators
##############################################################################

# N addition
nadd_rootshoot_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "rootshoot" & 
                                   !is.na(gs_mat)))

nadd_rootshoot_clim_summary <- data.frame(trait = "rootshoot",
                                    nut_add = "n",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(nadd_rootshoot_clim)),
                                    row.names = NULL)

# P addition
padd_rootshoot_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "rootshoot" & 
                                   !is.na(gs_mat)))

padd_rootshoot_clim_summary <- data.frame(trait = "rmf",
                                    nut_add = "p",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(padd_rootshoot_clim)),
                                    row.names = NULL)

# N+P addition
npadd_rootshoot_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "rootshoot" & 
                                    !is.na(gs_mat)))

npadd_rootshoot_clim_summary <- data.frame(trait = "rootshoot",
                                     nut_add = "np",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_rootshoot_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
rootshoot_clim_summary <- rbind(nadd_rootshoot_clim_summary, 
                          padd_rootshoot_clim_summary, 
                          npadd_rootshoot_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Merge climate moderators and write to .csv
##############################################################################
marea_clim_summary %>%
  full_join(nmass_clim_summary) %>%
  full_join(narea_clim_summary) %>%
  full_join(pmass_clim_summary) %>%
  full_join(parea_clim_summary) %>%
  full_join(tbio_clim_summary) %>%
  full_join(agb_clim_summary) %>%
  full_join(bgb_clim_summary) %>%
  full_join(rmf_clim_summary) %>%
  full_join(rootshoot_clim_summary) %>%
  write.csv("../data/CNPmeta_clim_moderators.csv", row.names = F)

##############################################################################
# Marea photosynthetic pathway
##############################################################################

# N addition
nadd_marea_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "lma" & 
                                     !is.na(photo_path)))

nadd_marea_photo_summary <- data.frame(trait = "marea", 
                                       nut_add = "n",
                                       mod_results(nadd_marea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_marea_photo))[2,3],
                                       p = coef(summary(nadd_marea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_marea_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "lma" & 
                                     !is.na(photo_path)))

padd_marea_photo_summary <- data.frame(trait = "marea", 
                                       nut_add = "p",
                                       mod_results(padd_marea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_marea_photo))[2,3],
                                       p = coef(summary(padd_marea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_marea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "lma" & 
                                      !is.na(photo_path)))

npadd_marea_photo_summary <- data.frame(trait = "marea", 
                                       nut_add = "np",
                                       mod_results(npadd_marea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_marea_photo))[2,3],
                                       p = coef(summary(npadd_marea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Marea moderator results, with some light cleaning
marea_photo_summary <- rbind(nadd_marea_photo_summary, 
                            padd_marea_photo_summary, 
                            npadd_marea_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Nmass photosynthetic pathway
##############################################################################

# N addition
nadd_nmass_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(photo_path)))

nadd_nmass_photo_summary <- data.frame(trait = "nmass", 
                                       nut_add = "n",
                                       mod_results(nadd_nmass_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_nmass_photo))[2,3],
                                       p = coef(summary(nadd_nmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_nmass_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(photo_path)))

padd_nmass_photo_summary <- data.frame(trait = "nmass", 
                                       nut_add = "p",
                                       mod_results(padd_nmass_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_nmass_photo))[2,3],
                                       p = coef(summary(padd_nmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_nmass_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_n_mass" & 
                                       !is.na(photo_path)))

npadd_nmass_photo_summary <- data.frame(trait = "nmass", 
                                        nut_add = "np",
                                        mod_results(npadd_nmass_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_nmass_photo))[2,3],
                                        p = coef(summary(npadd_nmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Nmass moderator results, with some light cleaning
nmass_photo_summary <- rbind(nadd_nmass_photo_summary, 
                             padd_nmass_photo_summary, 
                             npadd_nmass_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Narea photosynthetic pathway
##############################################################################

# N addition
nadd_narea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(photo_path)))

nadd_narea_photo_summary <- data.frame(trait = "narea", 
                                       nut_add = "n",
                                       mod_results(nadd_narea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_narea_photo))[2,3],
                                       p = coef(summary(nadd_narea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_narea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(photo_path)))

padd_narea_photo_summary <- data.frame(trait = "narea", 
                                       nut_add = "p",
                                       mod_results(padd_narea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_narea_photo))[2,3],
                                       p = coef(summary(padd_narea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_narea_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_n_area" & 
                                       !is.na(photo_path)))

npadd_narea_photo_summary <- data.frame(trait = "narea", 
                                        nut_add = "np",
                                        mod_results(npadd_narea_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_narea_photo))[2,3],
                                        p = coef(summary(npadd_narea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Nmass moderator results, with some light cleaning
narea_photo_summary <- rbind(nadd_narea_photo_summary, 
                             padd_narea_photo_summary, 
                             npadd_narea_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Pmass photosynthetic pathway
##############################################################################

# N addition
nadd_pmass_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(photo_path)))

nadd_pmass_photo_summary <- data.frame(trait = "pmass", 
                                       nut_add = "n",
                                       mod_results(nadd_pmass_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_pmass_photo))[2,3],
                                       p = coef(summary(nadd_pmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_pmass_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(photo_path)))

padd_pmass_photo_summary <- data.frame(trait = "pmass", 
                                       nut_add = "p",
                                       mod_results(padd_pmass_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_pmass_photo))[2,3],
                                       p = coef(summary(padd_pmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_pmass_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_p_mass" & 
                                       !is.na(photo_path)))

npadd_pmass_photo_summary <- data.frame(trait = "pmass", 
                                        nut_add = "np",
                                        mod_results(npadd_pmass_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_pmass_photo))[2,3],
                                        p = coef(summary(npadd_pmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Pmass moderator results, with some light cleaning
pmass_photo_summary <- rbind(nadd_pmass_photo_summary, 
                             padd_pmass_photo_summary, 
                             npadd_pmass_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Parea photosynthetic pathway
##############################################################################

# N addition
nadd_parea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(photo_path)))

nadd_parea_photo_summary <- data.frame(trait = "parea", 
                                       nut_add = "n",
                                       mod_results(nadd_parea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_parea_photo))[2,3],
                                       p = coef(summary(nadd_parea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_parea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(photo_path)))

padd_parea_photo_summary <- data.frame(trait = "parea", 
                                       nut_add = "p",
                                       mod_results(padd_parea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_parea_photo))[2,3],
                                       p = coef(summary(padd_parea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_parea_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_p_area" & 
                                       !is.na(photo_path)))

npadd_parea_photo_summary <- data.frame(trait = "parea", 
                                        nut_add = "np",
                                        mod_results(npadd_parea_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_parea_photo))[2,3],
                                        p = coef(summary(npadd_parea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Parea moderator results, with some light cleaning
parea_photo_summary <- rbind(nadd_parea_photo_summary, 
                             padd_parea_photo_summary, 
                             npadd_parea_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Asat photosynthetic pathway
##############################################################################

# N addition
nadd_asat_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "asat" & 
                                      !is.na(photo_path)))

nadd_asat_photo_summary <- data.frame(trait = "asat", 
                                       nut_add = "n",
                                       mod_results(nadd_asat_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_asat_photo))[2,3],
                                       p = coef(summary(nadd_asat_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_asat_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "asat" & 
                                      !is.na(photo_path)))

padd_asat_photo_summary <- data.frame(trait = "asat", 
                                       nut_add = "p",
                                       mod_results(padd_asat_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_asat_photo))[2,3],
                                       p = coef(summary(padd_asat_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_asat_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "asat" & 
                                       !is.na(photo_path)))

npadd_asat_photo_summary <- data.frame(trait = "asat", 
                                        nut_add = "np",
                                        mod_results(npadd_asat_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_asat_photo))[2,3],
                                        p = coef(summary(npadd_asat_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Asat moderator results, with some light cleaning
asat_photo_summary <- rbind(nadd_asat_photo_summary, 
                             padd_asat_photo_summary, 
                             npadd_asat_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Vcmax photosynthetic pathway
##############################################################################

# N addition
nadd_vcmax_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "vcmax" & 
                                     !is.na(photo_path)))

nadd_vcmax_photo_summary <- data.frame(trait = "vcmax", 
                                      nut_add = "n",
                                      mod_results(nadd_vcmax_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_vcmax_photo))[2,3],
                                      p = coef(summary(nadd_vcmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_vcmax_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "vcmax" & 
                                     !is.na(photo_path)))

padd_vcmax_photo_summary <- data.frame(trait = "vcmax", 
                                      nut_add = "p",
                                      mod_results(padd_vcmax_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_vcmax_photo))[2,3],
                                      p = coef(summary(padd_vcmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_vcmax_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "vcmax" & 
                                      !is.na(photo_path)))

npadd_vcmax_photo_summary <- data.frame(trait = "vcmax", 
                                       nut_add = "np",
                                       mod_results(npadd_vcmax_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_vcmax_photo))[2,3],
                                       p = coef(summary(npadd_vcmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Vcmax moderator results, with some light cleaning
vcmax_photo_summary <- rbind(nadd_vcmax_photo_summary, 
                            padd_vcmax_photo_summary, 
                            npadd_vcmax_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Jmax photosynthetic pathway
##############################################################################

# N addition
nadd_jmax_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "jmax" & 
                                      !is.na(photo_path)))

nadd_jmax_photo_summary <- data.frame(trait = "jmax", 
                                       nut_add = "n",
                                       mod_results(nadd_jmax_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_jmax_photo))[2,3],
                                       p = coef(summary(nadd_jmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_jmax_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "jmax" & 
                                      !is.na(photo_path)))

padd_jmax_photo_summary <- data.frame(trait = "jmax", 
                                       nut_add = "p",
                                       mod_results(padd_jmax_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_jmax_photo))[2,3],
                                       p = coef(summary(padd_jmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_jmax_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "jmax" & 
                                       !is.na(photo_path)))

npadd_jmax_photo_summary <- data.frame(trait = "jmax", 
                                        nut_add = "np",
                                        mod_results(npadd_jmax_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_jmax_photo))[2,3],
                                        p = coef(summary(npadd_jmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
jmax_photo_summary <- rbind(nadd_jmax_photo_summary, 
                             padd_jmax_photo_summary, 
                             npadd_jmax_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# PNUE photosynthetic pathway
##############################################################################

# N addition
nadd_pnue_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_pnue" & 
                                     !is.na(photo_path)))

nadd_pnue_photo_summary <- data.frame(trait = "pnue", 
                                      nut_add = "n",
                                      mod_results(nadd_pnue_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_pnue_photo))[2,3],
                                      p = coef(summary(nadd_pnue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_pnue_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_pnue" & 
                                     !is.na(photo_path)))

padd_pnue_photo_summary <- data.frame(trait = "pnue", 
                                      nut_add = "p",
                                      mod_results(padd_pnue_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_pnue_photo))[2,3],
                                      p = coef(summary(padd_pnue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_pnue_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_pnue" & 
                                      !is.na(photo_path)))

npadd_pnue_photo_summary <- data.frame(trait = "pnue", 
                                       nut_add = "np",
                                       mod_results(npadd_pnue_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_pnue_photo))[2,3],
                                       p = coef(summary(npadd_pnue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
pnue_photo_summary <- rbind(nadd_pnue_photo_summary, 
                            padd_pnue_photo_summary, 
                            npadd_pnue_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# PPUE photosynthetic pathway
##############################################################################

# N addition
nadd_ppue_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_ppue" & 
                                     !is.na(photo_path)))

nadd_ppue_photo_summary <- data.frame(trait = "ppue", 
                                      nut_add = "n",
                                      mod_results(nadd_ppue_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_ppue_photo))[2,3],
                                      p = coef(summary(nadd_ppue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_ppue_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_ppue" & 
                                     !is.na(photo_path)))

padd_ppue_photo_summary <- data.frame(trait = "ppue", 
                                      nut_add = "p",
                                      mod_results(padd_ppue_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_ppue_photo))[2,3],
                                      p = coef(summary(padd_ppue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_ppue_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_ppue" & 
                                      !is.na(photo_path)))

npadd_ppue_photo_summary <- data.frame(trait = "ppue", 
                                       nut_add = "np",
                                       mod_results(npadd_ppue_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_ppue_photo))[2,3],
                                       p = coef(summary(npadd_ppue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
ppue_photo_summary <- rbind(nadd_ppue_photo_summary, 
                            padd_ppue_photo_summary, 
                            npadd_ppue_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Total biomass photosynthetic pathway
##############################################################################

# N addition
nadd_tbio_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "total_biomass" & 
                                     !is.na(photo_path)))

nadd_tbio_photo_summary <- data.frame(trait = "tbio", 
                                      nut_add = "n",
                                      mod_results(nadd_tbio_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_tbio_photo))[2,3],
                                      p = coef(summary(nadd_tbio_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_tbio_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "total_biomass" & 
                                     !is.na(photo_path)))

padd_tbio_photo_summary <- data.frame(trait = "tbio", 
                                      nut_add = "p",
                                      mod_results(padd_tbio_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_tbio_photo))[2,3],
                                      p = coef(summary(padd_tbio_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_tbio_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "total_biomass" & 
                                      !is.na(photo_path)))

npadd_tbio_photo_summary <- data.frame(trait = "tbio", 
                                       nut_add = "np",
                                       mod_results(npadd_tbio_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_tbio_photo))[2,3],
                                       p = coef(summary(npadd_tbio_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge total biomass moderator results, with some light cleaning
tbio_photo_summary <- rbind(nadd_tbio_photo_summary, 
                            padd_tbio_photo_summary, 
                            npadd_tbio_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Aboveground biomass photosynthetic pathway
##############################################################################

# N addition
nadd_agb_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "agb" & 
                                     !is.na(photo_path)))

nadd_agb_photo_summary <- data.frame(trait = "agb", 
                                      nut_add = "n",
                                      mod_results(nadd_agb_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_agb_photo))[2,3],
                                      p = coef(summary(nadd_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_agb_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "agb" & 
                                     !is.na(photo_path)))

padd_agb_photo_summary <- data.frame(trait = "agb", 
                                      nut_add = "p",
                                      mod_results(padd_agb_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_agb_photo))[2,3],
                                      p = coef(summary(padd_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_agb_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "agb" & 
                                      !is.na(photo_path)))

npadd_agb_photo_summary <- data.frame(trait = "agb", 
                                       nut_add = "np",
                                       mod_results(npadd_agb_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_agb_photo))[2,3],
                                       p = coef(summary(npadd_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge aboveground biomass moderator results, with some light cleaning
agb_photo_summary <- rbind(nadd_agb_photo_summary, 
                            padd_agb_photo_summary, 
                            npadd_agb_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Belowground biomass photosynthetic pathway
##############################################################################

# N addition
nadd_bgb_photo <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "bgb" & 
                                    !is.na(photo_path)))

nadd_bgb_photo_summary <- data.frame(trait = "bgb", 
                                     nut_add = "n",
                                     mod_results(nadd_bgb_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(nadd_bgb_photo))[2,3],
                                     p = coef(summary(nadd_bgb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_bgb_photo <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "bgb" & 
                                    !is.na(photo_path)))

padd_bgb_photo_summary <- data.frame(trait = "bgb", 
                                     nut_add = "p",
                                     mod_results(padd_bgb_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(padd_bgb_photo))[2,3],
                                     p = coef(summary(padd_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_bgb_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "bgb" & 
                                     !is.na(photo_path)))

npadd_bgb_photo_summary <- data.frame(trait = "bgb", 
                                      nut_add = "np",
                                      mod_results(npadd_bgb_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(npadd_bgb_photo))[2,3],
                                      p = coef(summary(npadd_bgb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge aboveground biomass moderator results, with some light cleaning
bgb_photo_summary <- rbind(nadd_bgb_photo_summary, 
                           padd_bgb_photo_summary, 
                           npadd_bgb_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))


##############################################################################
# Merge photo moderator results and write to .csv
##############################################################################
marea_photo_summary %>%
  full_join(nmass_photo_summary) %>%
  full_join(narea_photo_summary) %>%
  full_join(pmass_photo_summary) %>%
  full_join(parea_photo_summary) %>%
  full_join(asat_photo_summary) %>%
  full_join(vcmax_photo_summary) %>%
  full_join(jmax_photo_summary) %>%
  full_join(pnue_photo_summary) %>%
  full_join(ppue_photo_summary) %>%
  full_join(tbio_photo_summary) %>%
  full_join(agb_photo_summary) %>%
  full_join(bgb_photo_summary) %>%
  write.csv("../data/CNPmeta_photo_moderators.csv", row.names = F)

##############################################################################
# Marea Nfixation
##############################################################################

# N addition
nadd_marea_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "lma" & 
                                      !is.na(n_fixer)))

nadd_marea_nfix_summary <- data.frame(trait = "marea", 
                                       nut_add = "n",
                                       mod_results(nadd_marea_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_marea_nfix))[2,3],
                                       p = coef(summary(nadd_marea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_marea_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "lma" & 
                                      !is.na(n_fixer)))

padd_marea_nfix_summary <- data.frame(trait = "marea", 
                                       nut_add = "p",
                                       mod_results(padd_marea_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_marea_nfix))[2,3],
                                       p = coef(summary(padd_marea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_marea_nfix <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ n_fixer,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "lma" & 
                                       !is.na(n_fixer)))

npadd_marea_nfix_summary <- data.frame(trait = "marea", 
                                        nut_add = "np",
                                        mod_results(npadd_marea_nfix, 
                                                    mod = "n_fixer", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_marea_nfix))[2,3],
                                        p = coef(summary(npadd_marea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Marea moderator results, with some light cleaning
marea_nfix_summary <- rbind(nadd_marea_nfix_summary, 
                             padd_marea_nfix_summary, 
                             npadd_marea_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Nmass Nfixation
##############################################################################

# N addition
nadd_nmass_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(n_fixer)))

nadd_nmass_nfix_summary <- data.frame(trait = "nmass", 
                                       nut_add = "n",
                                       mod_results(nadd_nmass_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_nmass_nfix))[2,3],
                                       p = coef(summary(nadd_nmass_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_nmass_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(n_fixer)))

padd_nmass_nfix_summary <- data.frame(trait = "nmass", 
                                       nut_add = "p",
                                       mod_results(padd_nmass_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_nmass_nfix))[2,3],
                                       p = coef(summary(padd_nmass_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_nmass_nfix <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ n_fixer,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_n_mass" & 
                                       !is.na(n_fixer)))

npadd_nmass_nfix_summary <- data.frame(trait = "nmass", 
                                        nut_add = "np",
                                        mod_results(npadd_nmass_nfix, 
                                                    mod = "n_fixer", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_nmass_nfix))[2,3],
                                        p = coef(summary(npadd_nmass_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Nmass moderator results, with some light cleaning
nmass_nfix_summary <- rbind(nadd_nmass_nfix_summary, 
                             padd_nmass_nfix_summary, 
                             npadd_nmass_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Narea Nfixation
##############################################################################

# N addition
nadd_narea_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(n_fixer)))

nadd_narea_nfix_summary <- data.frame(trait = "narea", 
                                       nut_add = "n",
                                       mod_results(nadd_narea_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_narea_nfix))[2,3],
                                       p = coef(summary(nadd_narea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_narea_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(n_fixer)))

padd_narea_nfix_summary <- data.frame(trait = "narea", 
                                       nut_add = "p",
                                       mod_results(padd_narea_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_narea_nfix))[2,3],
                                       p = coef(summary(padd_narea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_narea_nfix <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ n_fixer,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_n_area" & 
                                       !is.na(n_fixer)))

npadd_narea_nfix_summary <- data.frame(trait = "narea", 
                                        nut_add = "np",
                                        mod_results(npadd_narea_nfix, 
                                                    mod = "n_fixer", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_narea_nfix))[2,3],
                                        p = coef(summary(npadd_narea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Nmass moderator results, with some light cleaning
narea_nfix_summary <- rbind(nadd_narea_nfix_summary, 
                             padd_narea_nfix_summary, 
                             npadd_narea_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Pmass Nfixation
##############################################################################

# N addition
nadd_pmass_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(n_fixer)))

nadd_pmass_nfix_summary <- data.frame(trait = "pmass", 
                                       nut_add = "n",
                                       mod_results(nadd_pmass_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_pmass_nfix))[2,3],
                                       p = coef(summary(nadd_pmass_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_pmass_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(n_fixer)))

padd_pmass_nfix_summary <- data.frame(trait = "pmass", 
                                       nut_add = "p",
                                       mod_results(padd_pmass_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_pmass_nfix))[2,3],
                                       p = coef(summary(padd_pmass_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_pmass_nfix <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ n_fixer,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_p_mass" & 
                                       !is.na(n_fixer)))

npadd_pmass_nfix_summary <- data.frame(trait = "pmass", 
                                        nut_add = "np",
                                        mod_results(npadd_pmass_nfix, 
                                                    mod = "n_fixer", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_pmass_nfix))[2,3],
                                        p = coef(summary(npadd_pmass_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Pmass moderator results, with some light cleaning
pmass_nfix_summary <- rbind(nadd_pmass_nfix_summary, 
                             padd_pmass_nfix_summary, 
                             npadd_pmass_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Parea Nfixation
##############################################################################

# N addition
nadd_parea_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(n_fixer)))

nadd_parea_nfix_summary <- data.frame(trait = "parea", 
                                       nut_add = "n",
                                       mod_results(nadd_parea_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_parea_nfix))[2,3],
                                       p = coef(summary(nadd_parea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_parea_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(n_fixer)))

padd_parea_nfix_summary <- data.frame(trait = "parea", 
                                       nut_add = "p",
                                       mod_results(padd_parea_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_parea_nfix))[2,3],
                                       p = coef(summary(padd_parea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_parea_nfix <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ n_fixer,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_p_area" & 
                                       !is.na(n_fixer)))

npadd_parea_nfix_summary <- data.frame(trait = "parea", 
                                        nut_add = "np",
                                        mod_results(npadd_parea_nfix, 
                                                    mod = "n_fixer", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_parea_nfix))[2,3],
                                        p = coef(summary(npadd_parea_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Parea moderator results, with some light cleaning
parea_nfix_summary <- rbind(nadd_parea_nfix_summary, 
                             padd_parea_nfix_summary, 
                             npadd_parea_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Asat Nfixation
##############################################################################

# N addition
nadd_asat_nfix <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "asat" & 
                                     !is.na(n_fixer)))

nadd_asat_nfix_summary <- data.frame(trait = "asat", 
                                      nut_add = "n",
                                      mod_results(nadd_asat_nfix, 
                                                  mod = "n_fixer", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_asat_nfix))[2,3],
                                      p = coef(summary(nadd_asat_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_asat_nfix <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "asat" & 
                                     !is.na(n_fixer)))

padd_asat_nfix_summary <- data.frame(trait = "asat", 
                                      nut_add = "p",
                                      mod_results(padd_asat_nfix, 
                                                  mod = "n_fixer", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_asat_nfix))[2,3],
                                      p = coef(summary(padd_asat_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_asat_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "asat" & 
                                      !is.na(n_fixer)))

npadd_asat_nfix_summary <- data.frame(trait = "asat", 
                                       nut_add = "np",
                                       mod_results(npadd_asat_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_asat_nfix))[2,3],
                                       p = coef(summary(npadd_asat_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Asat moderator results, with some light cleaning
asat_nfix_summary <- rbind(nadd_asat_nfix_summary, 
                            padd_asat_nfix_summary, 
                            npadd_asat_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Vcmax Nfixation
##############################################################################

# N addition
nadd_vcmax_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "vcmax" & 
                                      !is.na(n_fixer)))

nadd_vcmax_nfix_summary <- data.frame(trait = "vcmax", 
                                       nut_add = "n",
                                       mod_results(nadd_vcmax_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_vcmax_nfix))[2,3],
                                       p = coef(summary(nadd_vcmax_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_vcmax_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "vcmax" & 
                                      !is.na(n_fixer)))

padd_vcmax_nfix_summary <- data.frame(trait = "vcmax", 
                                       nut_add = "p",
                                       mod_results(padd_vcmax_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_vcmax_nfix))[2,3],
                                       p = coef(summary(padd_vcmax_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_vcmax_nfix <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ n_fixer,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "vcmax" & 
                                       !is.na(n_fixer)))

npadd_vcmax_nfix_summary <- data.frame(trait = "vcmax", 
                                        nut_add = "np",
                                        mod_results(npadd_vcmax_nfix, 
                                                    mod = "n_fixer", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_vcmax_nfix))[2,3],
                                        p = coef(summary(npadd_vcmax_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Vcmax moderator results, with some light cleaning
vcmax_nfix_summary <- rbind(nadd_vcmax_nfix_summary, 
                             padd_vcmax_nfix_summary, 
                             npadd_vcmax_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Jmax Nfixation
##############################################################################

# N addition
nadd_jmax_nfix <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "jmax" & 
                                     !is.na(n_fixer)))

nadd_jmax_nfix_summary <- data.frame(trait = "jmax", 
                                      nut_add = "n",
                                      mod_results(nadd_jmax_nfix, 
                                                  mod = "n_fixer", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_jmax_nfix))[2,3],
                                      p = coef(summary(nadd_jmax_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_jmax_nfix <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "jmax" & 
                                     !is.na(n_fixer)))

padd_jmax_nfix_summary <- data.frame(trait = "jmax", 
                                      nut_add = "p",
                                      mod_results(padd_jmax_nfix, 
                                                  mod = "n_fixer", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_jmax_nfix))[2,3],
                                      p = coef(summary(padd_jmax_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_jmax_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "jmax" & 
                                      !is.na(n_fixer)))

npadd_jmax_nfix_summary <- data.frame(trait = "jmax", 
                                       nut_add = "np",
                                       mod_results(npadd_jmax_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_jmax_nfix))[2,3],
                                       p = coef(summary(npadd_jmax_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
jmax_nfix_summary <- rbind(nadd_jmax_nfix_summary, 
                            padd_jmax_nfix_summary, 
                            npadd_jmax_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# PNUE Nfixation
##############################################################################

# N addition
nadd_pnue_nfix <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_pnue" & 
                                     !is.na(n_fixer)))

nadd_pnue_nfix_summary <- data.frame(trait = "pnue", 
                                      nut_add = "n",
                                      mod_results(nadd_pnue_nfix, 
                                                  mod = "n_fixer", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_pnue_nfix))[2,3],
                                      p = coef(summary(nadd_pnue_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_pnue_nfix <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_pnue" & 
                                     !is.na(n_fixer)))

padd_pnue_nfix_summary <- data.frame(trait = "pnue", 
                                      nut_add = "p",
                                      mod_results(padd_pnue_nfix, 
                                                  mod = "n_fixer", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_pnue_nfix))[2,3],
                                      p = coef(summary(padd_pnue_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_pnue_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_pnue" & 
                                      !is.na(n_fixer)))

npadd_pnue_nfix_summary <- data.frame(trait = "pnue", 
                                       nut_add = "np",
                                       mod_results(npadd_pnue_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_pnue_nfix))[2,3],
                                       p = coef(summary(npadd_pnue_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge PNUE moderator results, with some light cleaning
pnue_nfix_summary <- rbind(nadd_pnue_nfix_summary, 
                            padd_pnue_nfix_summary, 
                            npadd_pnue_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# PPUE Nfixation
##############################################################################

# N addition
nadd_ppue_nfix <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_ppue" & 
                                     !is.na(n_fixer)))

nadd_ppue_nfix_summary <- data.frame(trait = "ppue", 
                                      nut_add = "n",
                                      mod_results(nadd_ppue_nfix, 
                                                  mod = "n_fixer", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_ppue_nfix))[2,3],
                                      p = coef(summary(nadd_ppue_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_ppue_nfix <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_ppue" & 
                                     !is.na(n_fixer)))

padd_ppue_nfix_summary <- data.frame(trait = "ppue", 
                                      nut_add = "p",
                                      mod_results(padd_ppue_nfix, 
                                                  mod = "n_fixer", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_ppue_nfix))[2,3],
                                      p = coef(summary(padd_ppue_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_ppue_nfix <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_ppue" & 
                                      !is.na(n_fixer)))

npadd_ppue_nfix_summary <- data.frame(trait = "ppue", 
                                       nut_add = "np",
                                       mod_results(npadd_ppue_nfix, 
                                                   mod = "n_fixer", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_ppue_nfix))[2,3],
                                       p = coef(summary(npadd_ppue_nfix))[2, 4]) %>%
  dplyr::select(trait, nut_add, nfix = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
ppue_nfix_summary <- rbind(nadd_ppue_nfix_summary, 
                            padd_ppue_nfix_summary, 
                            npadd_ppue_nfix_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Merge nfix moderator results and write to .csv
##############################################################################
marea_nfix_summary %>%
  full_join(nmass_nfix_summary) %>%
  full_join(narea_nfix_summary) %>%
  full_join(pmass_nfix_summary) %>%
  full_join(parea_nfix_summary) %>%
  full_join(asat_nfix_summary) %>%
  full_join(vcmax_nfix_summary) %>%
  full_join(jmax_nfix_summary) %>%
  full_join(pnue_nfix_summary) %>%
  full_join(ppue_nfix_summary) %>%
  write.csv("../data/CNPmeta_nfix_moderators.csv", row.names = F)


##############################################################################
# Marea mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_marea_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "lma" & 
                                     !is.na(myc_nas)))

nadd_marea_myc_summary <- data.frame(trait = "marea", 
                                      nut_add = "n",
                                      mod_results(nadd_marea_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_marea_myc))[2,3],
                                      p = coef(summary(nadd_marea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_marea_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "lma" & 
                                     !is.na(myc_nas)))

padd_marea_myc_summary <- data.frame(trait = "marea", 
                                      nut_add = "p",
                                      mod_results(padd_marea_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_marea_myc))[2,3],
                                      p = coef(summary(padd_marea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_marea_myc <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ myc_nas,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "lma" & 
                                      !is.na(myc_nas)))

npadd_marea_myc_summary <- data.frame(trait = "marea", 
                                       nut_add = "np",
                                       mod_results(npadd_marea_myc, 
                                                   mod = "myc_nas", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_marea_myc))[2,3],
                                       p = coef(summary(npadd_marea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Marea moderator results, with some light cleaning
marea_myc_summary <- rbind(nadd_marea_myc_summary, 
                            padd_marea_myc_summary, 
                            npadd_marea_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Nmass mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_nmass_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(myc_nas)))

nadd_nmass_myc_summary <- data.frame(trait = "nmass", 
                                      nut_add = "n",
                                      mod_results(nadd_nmass_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_nmass_myc))[2,3],
                                      p = coef(summary(nadd_nmass_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_nmass_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(myc_nas)))

padd_nmass_myc_summary <- data.frame(trait = "nmass", 
                                      nut_add = "p",
                                      mod_results(padd_nmass_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_nmass_myc))[2,3],
                                      p = coef(summary(padd_nmass_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_nmass_myc <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ myc_nas,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(myc_nas)))

npadd_nmass_myc_summary <- data.frame(trait = "nmass", 
                                       nut_add = "np",
                                       mod_results(npadd_nmass_myc, 
                                                   mod = "myc_nas", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_nmass_myc))[2,3],
                                       p = coef(summary(npadd_nmass_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Nmass moderator results, with some light cleaning
nmass_myc_summary <- rbind(nadd_nmass_myc_summary, 
                            padd_nmass_myc_summary, 
                            npadd_nmass_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Narea mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_narea_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(myc_nas)))

nadd_narea_myc_summary <- data.frame(trait = "narea", 
                                      nut_add = "n",
                                      mod_results(nadd_narea_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_narea_myc))[2,3],
                                      p = coef(summary(nadd_narea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_narea_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(myc_nas)))

padd_narea_myc_summary <- data.frame(trait = "narea", 
                                      nut_add = "p",
                                      mod_results(padd_narea_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_narea_myc))[2,3],
                                      p = coef(summary(padd_narea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_narea_myc <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ myc_nas,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(myc_nas)))

npadd_narea_myc_summary <- data.frame(trait = "narea", 
                                       nut_add = "np",
                                       mod_results(npadd_narea_myc, 
                                                   mod = "myc_nas", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_narea_myc))[2,3],
                                       p = coef(summary(npadd_narea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Nmass moderator results, with some light cleaning
narea_myc_summary <- rbind(nadd_narea_myc_summary, 
                            padd_narea_myc_summary, 
                            npadd_narea_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Pmass mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_pmass_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(myc_nas)))

nadd_pmass_myc_summary <- data.frame(trait = "pmass", 
                                      nut_add = "n",
                                      mod_results(nadd_pmass_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_pmass_myc))[2,3],
                                      p = coef(summary(nadd_pmass_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_pmass_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(myc_nas)))

padd_pmass_myc_summary <- data.frame(trait = "pmass", 
                                      nut_add = "p",
                                      mod_results(padd_pmass_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_pmass_myc))[2,3],
                                      p = coef(summary(padd_pmass_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_pmass_myc <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ myc_nas,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(myc_nas)))

npadd_pmass_myc_summary <- data.frame(trait = "pmass", 
                                       nut_add = "np",
                                       mod_results(npadd_pmass_myc, 
                                                   mod = "myc_nas", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_pmass_myc))[2,3],
                                       p = coef(summary(npadd_pmass_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Pmass moderator results, with some light cleaning
pmass_myc_summary <- rbind(nadd_pmass_myc_summary, 
                            padd_pmass_myc_summary, 
                            npadd_pmass_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Parea mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_parea_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(myc_nas)))

nadd_parea_myc_summary <- data.frame(trait = "parea", 
                                      nut_add = "n",
                                      mod_results(nadd_parea_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_parea_myc))[2,3],
                                      p = coef(summary(nadd_parea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_parea_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(myc_nas)))

padd_parea_myc_summary <- data.frame(trait = "parea", 
                                      nut_add = "p",
                                      mod_results(padd_parea_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_parea_myc))[2,3],
                                      p = coef(summary(padd_parea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_parea_myc <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ myc_nas,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(myc_nas)))

npadd_parea_myc_summary <- data.frame(trait = "parea", 
                                       nut_add = "np",
                                       mod_results(npadd_parea_myc, 
                                                   mod = "myc_nas", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_parea_myc))[2,3],
                                       p = coef(summary(npadd_parea_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Parea moderator results, with some light cleaning
parea_myc_summary <- rbind(nadd_parea_myc_summary, 
                            padd_parea_myc_summary, 
                            npadd_parea_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Asat mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_asat_myc <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ myc_nas,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "asat" & 
                                    !is.na(myc_nas)))

nadd_asat_myc_summary <- data.frame(trait = "asat", 
                                     nut_add = "n",
                                     mod_results(nadd_asat_myc, 
                                                 mod = "myc_nas", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(nadd_asat_myc))[2,3],
                                     p = coef(summary(nadd_asat_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_asat_myc <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ myc_nas,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "asat" & 
                                    !is.na(myc_nas)))

padd_asat_myc_summary <- data.frame(trait = "asat", 
                                     nut_add = "p",
                                     mod_results(padd_asat_myc, 
                                                 mod = "myc_nas", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(padd_asat_myc))[2,3],
                                     p = coef(summary(padd_asat_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_asat_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "asat" & 
                                     !is.na(myc_nas)))

npadd_asat_myc_summary <- data.frame(trait = "asat", 
                                      nut_add = "np",
                                      mod_results(npadd_asat_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(npadd_asat_myc))[2,3],
                                      p = coef(summary(npadd_asat_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Asat moderator results, with some light cleaning
asat_myc_summary <- rbind(nadd_asat_myc_summary, 
                           padd_asat_myc_summary, 
                           npadd_asat_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Vcmax mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_vcmax_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "vcmax" & 
                                     !is.na(myc_nas)))

nadd_vcmax_myc_summary <- data.frame(trait = "vcmax", 
                                      nut_add = "n",
                                      mod_results(nadd_vcmax_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_vcmax_myc))[2,3],
                                      p = coef(summary(nadd_vcmax_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_vcmax_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "vcmax" & 
                                     !is.na(myc_nas)))

padd_vcmax_myc_summary <- data.frame(trait = "vcmax", 
                                      nut_add = "p",
                                      mod_results(padd_vcmax_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_vcmax_myc))[2,3],
                                      p = coef(summary(padd_vcmax_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_vcmax_myc <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ myc_nas,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "vcmax" & 
                                      !is.na(myc_nas)))

npadd_vcmax_myc_summary <- data.frame(trait = "vcmax", 
                                       nut_add = "np",
                                       mod_results(npadd_vcmax_myc, 
                                                   mod = "myc_nas", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_vcmax_myc))[2,3],
                                       p = coef(summary(npadd_vcmax_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Vcmax moderator results, with some light cleaning
vcmax_myc_summary <- rbind(nadd_vcmax_myc_summary, 
                            padd_vcmax_myc_summary, 
                            npadd_vcmax_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Jmax mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_jmax_myc <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ myc_nas,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "jmax" & 
                                    !is.na(myc_nas)))

nadd_jmax_myc_summary <- data.frame(trait = "jmax", 
                                     nut_add = "n",
                                     mod_results(nadd_jmax_myc, 
                                                 mod = "myc_nas", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(nadd_jmax_myc))[2,3],
                                     p = coef(summary(nadd_jmax_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_jmax_myc <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ myc_nas,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "jmax" & 
                                    !is.na(myc_nas)))

padd_jmax_myc_summary <- data.frame(trait = "jmax", 
                                     nut_add = "p",
                                     mod_results(padd_jmax_myc, 
                                                 mod = "myc_nas", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(padd_jmax_myc))[2,3],
                                     p = coef(summary(padd_jmax_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_jmax_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "jmax" & 
                                     !is.na(myc_nas)))

npadd_jmax_myc_summary <- data.frame(trait = "jmax", 
                                      nut_add = "np",
                                      mod_results(npadd_jmax_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(npadd_jmax_myc))[2,3],
                                      p = coef(summary(npadd_jmax_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
jmax_myc_summary <- rbind(nadd_jmax_myc_summary, 
                           padd_jmax_myc_summary, 
                           npadd_jmax_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# PNUE mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_pnue_myc <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ myc_nas,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "leaf_pnue" & 
                                    !is.na(myc_nas)))

nadd_pnue_myc_summary <- data.frame(trait = "pnue", 
                                     nut_add = "n",
                                     mod_results(nadd_pnue_myc, 
                                                 mod = "myc_nas", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(nadd_pnue_myc))[2,3],
                                     p = coef(summary(nadd_pnue_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_pnue_myc <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ myc_nas,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "leaf_pnue" & 
                                    !is.na(myc_nas)))

padd_pnue_myc_summary <- data.frame(trait = "pnue", 
                                     nut_add = "p",
                                     mod_results(padd_pnue_myc, 
                                                 mod = "myc_nas", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(padd_pnue_myc))[2,3],
                                     p = coef(summary(padd_pnue_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_pnue_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "leaf_pnue" & 
                                     !is.na(myc_nas)))

npadd_pnue_myc_summary <- data.frame(trait = "pnue", 
                                      nut_add = "np",
                                      mod_results(npadd_pnue_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(npadd_pnue_myc))[2,3],
                                      p = coef(summary(npadd_pnue_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge PNUE moderator results, with some light cleaning
pnue_myc_summary <- rbind(nadd_pnue_myc_summary, 
                           padd_pnue_myc_summary, 
                           npadd_pnue_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# PPUE mycorrhizal acquisition strategy
##############################################################################

# N addition
nadd_ppue_myc <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ myc_nas,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "leaf_ppue" & 
                                    !is.na(myc_nas)))

nadd_ppue_myc_summary <- data.frame(trait = "ppue", 
                                     nut_add = "n",
                                     mod_results(nadd_ppue_myc, 
                                                 mod = "myc_nas", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(nadd_ppue_myc))[2,3],
                                     p = coef(summary(nadd_ppue_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_ppue_myc <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ myc_nas,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "leaf_ppue" & 
                                    !is.na(myc_nas)))

padd_ppue_myc_summary <- data.frame(trait = "ppue", 
                                     nut_add = "p",
                                     mod_results(padd_ppue_myc, 
                                                 mod = "myc_nas", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(padd_ppue_myc))[2,3],
                                     p = coef(summary(padd_ppue_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_ppue_myc <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ myc_nas,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "leaf_ppue" & 
                                     !is.na(myc_nas)))

npadd_ppue_myc_summary <- data.frame(trait = "ppue", 
                                      nut_add = "np",
                                      mod_results(npadd_ppue_myc, 
                                                  mod = "myc_nas", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(npadd_ppue_myc))[2,3],
                                      p = coef(summary(npadd_ppue_myc))[2, 4]) %>%
  dplyr::select(trait, nut_add, myc = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
ppue_myc_summary <- rbind(nadd_ppue_myc_summary, 
                           padd_ppue_myc_summary, 
                           npadd_ppue_myc_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Merge myc moderator results and write to .csv
##############################################################################
marea_myc_summary %>%
  full_join(nmass_myc_summary) %>%
  full_join(narea_myc_summary) %>%
  full_join(pmass_myc_summary) %>%
  full_join(parea_myc_summary) %>%
  full_join(asat_myc_summary) %>%
  full_join(vcmax_myc_summary) %>%
  full_join(jmax_myc_summary) %>%
  full_join(pnue_myc_summary) %>%
  full_join(ppue_myc_summary) %>%
  write.csv("../data/CNPmeta_myc_moderators.csv", row.names = F)





