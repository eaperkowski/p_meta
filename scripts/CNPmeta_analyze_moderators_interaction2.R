# Libraries
library(tidyverse)
library(metafor)
library(ggpubr)
library(naniar) # to resolve NA/<NA> issue
library(orchaRd)

# Load meta-analysis results
meta_results <- read.csv("../data/CNPmeta_logr_results.csv")
meta_results_int <- read.csv("../data/CNPmeta_logr_results_int.csv")

# Check file structure
head(meta_results)
head(meta_results_int)

##############################################################################
# Marea climate moderators
##############################################################################

# N addition
nadd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat)))

nadd_marea_clim_summary <- data.frame(trait = "marea",
                                      nut_add = "n",
                                      k = 81,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_marea_clim)),
                                      row.names = NULL)

# P addition
padd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat)))

padd_marea_clim_summary <- data.frame(trait = "marea",
                                      nut_add = "p",
                                      k = 81,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_marea_clim)),
                                      row.names = NULL)

# N+P addition
npadd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat)))

npadd_marea_clim_summary <- data.frame(trait = "marea",
                                       nut_add = "np",
                                       k = 81,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_marea_clim)),
                                       row.names = NULL)

## Merge Marea moderator results, with some light cleaning
marea_clim_summary <- rbind(nadd_marea_clim_summary, 
      padd_marea_clim_summary, 
      npadd_marea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Nmass climate moderators
##############################################################################

# N addition
nadd_nmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(gs_mat)))

nadd_nmass_clim_summary <- data.frame(trait = "nmass",
                                      nut_add = "n",
                                      k = 105,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_nmass_clim)),
                                      row.names = NULL)

# P addition
padd_nmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(gs_mat)))

padd_nmass_clim_summary <- data.frame(trait = "nmass",
                                      nut_add = "p",
                                      k = 105,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_nmass_clim)),
                                      row.names = NULL)

# N+P addition
npadd_nmass_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(gs_mat)))

npadd_nmass_clim_summary <- data.frame(trait = "nmass",
                                       nut_add = "np",
                                       k = 105,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_nmass_clim)),
                                       row.names = NULL)

## Merge Nmass moderator results, with some light cleaning
nmass_clim_summary <- rbind(nadd_nmass_clim_summary, 
                            padd_nmass_clim_summary, 
                            npadd_nmass_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Narea climate moderators
##############################################################################

# N addition
nadd_narea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(gs_mat)))

nadd_narea_clim_summary <- data.frame(trait = "narea",
                                      nut_add = "n",
                                      k = 52,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_narea_clim)),
                                      row.names = NULL)

# P addition
padd_narea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(gs_mat)))

padd_narea_clim_summary <- data.frame(trait = "narea",
                                      nut_add = "p",
                                      k = 52,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_narea_clim)),
                                      row.names = NULL)

# N+P addition
npadd_narea_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(gs_mat)))

npadd_narea_clim_summary <- data.frame(trait = "narea",
                                       nut_add = "np",
                                       k = 52,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_narea_clim)),
                                       row.names = NULL)

## Merge Narea moderator results, with some light cleaning
narea_clim_summary <- rbind(nadd_narea_clim_summary, 
                            padd_narea_clim_summary, 
                            npadd_narea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Pmass climate moderators
##############################################################################

# N addition
nadd_pmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(gs_mat)))

nadd_pmass_clim_summary <- data.frame(trait = "pmass",
                                      nut_add = "n",
                                      k = 99,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_pmass_clim)),
                                      row.names = NULL)

# P addition
padd_pmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(gs_mat)))

padd_pmass_clim_summary <- data.frame(trait = "pmass",
                                      nut_add = "p",
                                      k = 99,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_pmass_clim)),
                                      row.names = NULL)

# N+P addition
npadd_pmass_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(gs_mat)))

npadd_pmass_clim_summary <- data.frame(trait = "pmass",
                                       nut_add = "np",
                                       k = 99,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_pmass_clim)),
                                       row.names = NULL)

## Merge Pmass moderator results, with some light cleaning
pmass_clim_summary <- rbind(nadd_pmass_clim_summary, 
                            padd_pmass_clim_summary, 
                            npadd_pmass_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Parea climate moderators
##############################################################################

# N addition
nadd_parea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(gs_mat)))

nadd_parea_clim_summary <- data.frame(trait = "parea",
                                      nut_add = "n",
                                      k = 47,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_parea_clim)),
                                      row.names = NULL)

# P addition
padd_parea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(gs_mat) & logr > -1))


padd_parea_clim_summary <- data.frame(trait = "parea",
                                      nut_add = "p",
                                      k = 47,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_parea_clim)),
                                      row.names = NULL)

# N+P addition
npadd_parea_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(gs_mat)))

npadd_parea_clim_summary <- data.frame(trait = "parea",
                                       nut_add = "np",
                                       k = 47,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_parea_clim)),
                                       row.names = NULL)

## Merge Parea moderator results, with some light cleaning
parea_clim_summary <- rbind(nadd_parea_clim_summary, 
                            padd_parea_clim_summary, 
                            npadd_parea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Leaf N:P climate moderators
##############################################################################

# N addition
nadd_leafnp_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_np" & 
                                     !is.na(gs_mat)))

nadd_leafnp_clim_summary <- data.frame(trait = "leaf_np",
                                      nut_add = "n",
                                      k = 86,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_leafnp_clim)),
                                      row.names = NULL)

# P addition
padd_leafnp_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_np" & 
                                     !is.na(gs_mat)))

padd_leafnp_clim_summary <- data.frame(trait = "leaf_np",
                                      nut_add = "p",
                                      k = 86,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_leafnp_clim)),
                                      row.names = NULL)

# N+P addition
npadd_leafnp_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_np" & 
                                      !is.na(gs_mat)))

npadd_leafnp_clim_summary <- data.frame(trait = "leaf_np",
                                       nut_add = "np",
                                       k = 86,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_leafnp_clim)),
                                       row.names = NULL)

## Merge Parea moderator results, with some light cleaning
leafnp_clim_summary <- rbind(nadd_leafnp_clim_summary, 
                            padd_leafnp_clim_summary, 
                            npadd_leafnp_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Total biomass climate moderators
##############################################################################

# N addition
nadd_tbio_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "tbio_gm2" & 
                                     !is.na(gs_mat)))

nadd_tbio_clim_summary <- data.frame(trait = "tbio_gm2",
                                     nut_add = "n",
                                     k = 31,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(nadd_tbio_clim)),
                                     row.names = NULL)

# P addition
padd_tbio_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "tbio_gm2" & 
                                    !is.na(gs_mat)))

padd_tbio_clim_summary <- data.frame(trait = "tbio",
                                     nut_add = "p",
                                     k = 31,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(padd_tbio_clim)),
                                     row.names = NULL)

# N+P addition
npadd_tbio_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "tbio_gm2" & 
                                      !is.na(gs_mat)))

npadd_tbio_clim_summary <- data.frame(trait = "tbio",
                                      nut_add = "np",
                                      k = 31,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(npadd_tbio_clim)),
                                      row.names = NULL)

## Merge Parea moderator results, with some light cleaning
tbio_clim_summary <- rbind(nadd_tbio_clim_summary, 
                            padd_tbio_clim_summary, 
                            npadd_tbio_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Aboveground biomass (g/m2) climate moderators
##############################################################################

# N addition
nadd_anpp_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "anpp" & 
                                    !is.na(gs_mat)))

nadd_agb_clim_summary <- data.frame(trait = "anpp",
                                    nut_add = "n",
                                    k = 113,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(nadd_anpp_clim)),
                                    row.names = NULL)

# P addition
padd_anpp_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "anpp" & 
                                    !is.na(gs_mat)))

padd_agb_clim_summary <- data.frame(trait = "anpp",
                                    nut_add = "p",
                                    k = 113,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(padd_anpp_clim)),
                                    row.names = NULL)

# N+P addition
npadd_anpp_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "anpp" & 
                                     !is.na(gs_mat)))

npadd_agb_clim_summary <- data.frame(trait = "anpp",
                                     nut_add = "np",
                                     k = 113,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_anpp_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
agb_clim_summary <- rbind(nadd_agb_clim_summary, 
                           padd_agb_clim_summary, 
                           npadd_agb_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Belowground biomass climate moderators
##############################################################################

# N addition
nadd_bnpp_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "bnpp" & 
                                   !is.na(gs_mat)))

nadd_bnpp_clim_summary <- data.frame(trait = "bnpp",
                                     nut_add = "n",
                                     k = 53,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(nadd_bnpp_clim)),
                                     row.names = NULL)

# P addition
padd_bnpp_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "bnpp" & 
                                   !is.na(gs_mat)))

padd_bnpp_clim_summary <- data.frame(trait = "bnpp",
                                     nut_add = "p",
                                     k = 53,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(padd_bnpp_clim)),
                                     row.names = NULL)

# N+P addition
npadd_bnpp_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "bnpp" & 
                                    !is.na(gs_mat)))

npadd_bnpp_clim_summary <- data.frame(trait = "bnpp",
                                      nut_add = "np",
                                      k = 53,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(npadd_bnpp_clim)),
                                      row.names = NULL)

## Merge Parea moderator results, with some light cleaning
bgb_clim_summary <- rbind(nadd_bnpp_clim_summary, 
                          padd_bnpp_clim_summary, 
                          npadd_bnpp_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

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
                                    k = 35,
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
                                    k = 35,
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
                                     k = 35,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_rmf_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
rmf_clim_summary <- rbind(nadd_rmf_clim_summary, 
                          padd_rmf_clim_summary, 
                          npadd_rmf_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

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
                                          k = 37,
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
                                          k = 37,
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
                                           k = 37,
                                           mod = c("intrcpt", "gs_mat",
                                                   "gs_ai", "gs_par"),
                                           coef(summary(npadd_rootshoot_clim)),
                                           row.names = NULL)

## Merge Parea moderator results, with some light cleaning
rootshoot_clim_summary <- rbind(nadd_rootshoot_clim_summary, 
                          padd_rootshoot_clim_summary, 
                          npadd_rootshoot_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Merge climate moderators and write to .csv
##############################################################################
marea_clim_summary %>%
  full_join(nmass_clim_summary) %>%
  full_join(narea_clim_summary) %>%
  full_join(pmass_clim_summary) %>%
  full_join(parea_clim_summary) %>%
  full_join(leafnp_clim_summary) %>%
  full_join(tbio_clim_summary) %>%
  full_join(agb_clim_summary) %>%
  full_join(bgb_clim_summary) %>%
  full_join(rmf_clim_summary) %>%
  full_join(rootshoot_clim_summary) %>%
  write.csv("../data/CNPmeta_clim_moderators.csv", 
            row.names = F)

##############################################################################
# Marea
##############################################################################

#############
# N addition
#############
# Model 
nadd_marea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "lma" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_marea_photo <- data.frame(trait = "marea", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_marea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_marea_pft))[2,3],
                               p = coef(summary(nadd_marea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_marea_myc <- data.frame(trait = "marea", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_marea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_marea_pft))[3,3],
                             p = coef(summary(nadd_marea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_marea_nfix <- data.frame(trait = "marea", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_marea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_marea_pft))[4,3],
                              p = coef(summary(nadd_marea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_marea_pft_results <- nadd_marea_photo %>% 
  rbind(nadd_marea_myc) %>% 
  rbind(nadd_marea_nfix)

#############
# P addition
#############
# Model
padd_marea_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "lma" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_marea_photo <- data.frame(trait = "marea", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_marea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_marea_pft))[2,3],
                               p = coef(summary(padd_marea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_marea_myc <- data.frame(trait = "marea", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_marea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_marea_pft))[3,3],
                             p = coef(summary(padd_marea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_marea_nfix <- data.frame(trait = "marea", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_marea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_marea_pft))[4,3],
                              p = coef(summary(padd_marea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_marea_pft_results <- padd_marea_photo %>% 
  rbind(padd_marea_myc) %>% 
  rbind(padd_marea_nfix)

#############
# N+P addition
#############
# Model
npadd_marea_pft <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path + myc_nas + n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "lma" & 
                                      !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_marea_photo <- data.frame(trait = "marea", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_marea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_marea_pft))[2,3],
                               p = coef(summary(npadd_marea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_marea_myc <- data.frame(trait = "marea", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_marea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_marea_pft))[3,3],
                             p = coef(summary(npadd_marea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_marea_nfix <- data.frame(trait = "marea", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_marea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_marea_pft))[4,3],
                              p = coef(summary(npadd_marea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_marea_pft_results <- npadd_marea_photo %>% 
  rbind(npadd_marea_myc) %>% 
  rbind(npadd_marea_nfix)

#############
# Merge Marea moderator results, with some light cleaning
#############
marea_pft_summary <- rbind(nadd_marea_pft_results, 
                             padd_marea_pft_results, 
                             npadd_marea_pft_results) %>%
  mutate(k = 113,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Nmass
##############################################################################

#############
# N addition
#############
# Model
nadd_nmass_pft <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path + myc_nas + n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_nmass_photo <- data.frame(trait = "nmass", 
                                nut_add = "n",
                                mod = "photo",
                                mod_results(nadd_nmass_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(nadd_nmass_pft))[2,3],
                                p = coef(summary(nadd_nmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_nmass_myc <- data.frame(trait = "nmass", 
                              nut_add = "n",
                              mod = "myc_nas",
                              mod_results(nadd_nmass_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_nmass_pft))[3,3],
                              p = coef(summary(nadd_nmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_nmass_nfix <- data.frame(trait = "nmass", 
                               nut_add = "n",
                               mod = "nfix",
                               mod_results(nadd_nmass_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_nmass_pft))[4,3],
                               p = coef(summary(nadd_nmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_nmass_pft_results <- nadd_nmass_photo %>% 
  rbind(nadd_nmass_myc) %>% 
  rbind(nadd_nmass_nfix)

##############
# P addition
##############
padd_nmass_pft <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path + myc_nas + n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_nmass_photo <- data.frame(trait = "nmass", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(padd_nmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_nmass_pft))[2,3],
                               p = coef(summary(padd_nmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_nmass_myc <- data.frame(trait = "nmass", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_nmass_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_nmass_pft))[3,3],
                             p = coef(summary(padd_nmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_nmass_nfix <- data.frame(trait = "nmass", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_nmass_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_nmass_pft))[4,3],
                              p = coef(summary(padd_nmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_nmass_pft_results <- padd_nmass_photo %>% 
  rbind(padd_nmass_myc) %>% 
  rbind(padd_nmass_nfix)

################
# N+P addition
################
npadd_nmass_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_nmass_photo <- data.frame(trait = "nmass", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_nmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_nmass_pft))[2,3],
                               p = coef(summary(npadd_nmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_nmass_myc <- data.frame(trait = "nmass", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_nmass_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_nmass_pft))[3,3],
                             p = coef(summary(npadd_nmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_nmass_nfix <- data.frame(trait = "nmass", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_nmass_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_nmass_pft))[4,3],
                              p = coef(summary(npadd_nmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_nmass_pft_results <- npadd_nmass_photo %>% 
  rbind(npadd_nmass_myc) %>% 
  rbind(npadd_nmass_nfix)

#############
# Merge Nmass moderator results, with some light cleaning
#############
nmass_pft_summary <- rbind(nadd_nmass_pft_results, 
                             padd_nmass_pft_results, 
                             npadd_nmass_pft_results) %>%
  mutate(k = 136,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Narea - pft
##############################################################################

##############
# N addition
##############
nadd_narea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "leaf_n_area" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_narea_photo <- data.frame(trait = "narea", 
                             nut_add = "n",
                             mod = "photo",
                             mod_results(nadd_narea_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_narea_pft))[2,3],
                             p = coef(summary(nadd_narea_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_narea_myc <- data.frame(trait = "narea", 
                              nut_add = "n",
                              mod = "myc_nas",
                              mod_results(nadd_narea_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_narea_pft))[3,3],
                              p = coef(summary(nadd_narea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_narea_nfix <- data.frame(trait = "narea", 
                               nut_add = "n",
                               mod = "nfix",
                               mod_results(nadd_narea_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_narea_pft))[4,3],
                               p = coef(summary(nadd_narea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_narea_pft_results <- nadd_narea_photo %>% 
  rbind(nadd_narea_myc) %>% 
  rbind(nadd_narea_nfix)


##############
# P addition
##############
padd_narea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "leaf_n_area" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_narea_photo <- data.frame(trait = "narea", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_narea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_narea_pft))[2,3],
                               p = coef(summary(padd_narea_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_narea_myc <- data.frame(trait = "narea", 
                              nut_add = "p",
                              mod = "myc_nas",
                              mod_results(padd_narea_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_narea_pft))[3,3],
                              p = coef(summary(padd_narea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_narea_nfix <- data.frame(trait = "narea", 
                               nut_add = "p",
                               mod = "nfix",
                               mod_results(padd_narea_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_narea_pft))[4,3],
                               p = coef(summary(padd_narea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_narea_pft_results <- padd_narea_photo %>% 
  rbind(padd_narea_myc) %>% 
  rbind(padd_narea_nfix)

################
# N+P addition
################
npadd_narea_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_narea_photo <- data.frame(trait = "narea", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_narea_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_narea_pft))[2,3],
                                p = coef(summary(npadd_narea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_narea_myc <- data.frame(trait = "narea", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_narea_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_narea_pft))[3,3],
                              p = coef(summary(npadd_narea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_narea_nfix <- data.frame(trait = "narea", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_narea_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_narea_pft))[4,3],
                               p = coef(summary(npadd_narea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_narea_pft_results <- npadd_narea_photo %>% 
  rbind(npadd_narea_myc) %>% 
  rbind(npadd_narea_nfix)

#############
# Merge Narea moderator results, with some light cleaning
#############
narea_pft_summary <- rbind(nadd_narea_pft_results, 
                           padd_narea_pft_results, 
                           npadd_narea_pft_results) %>%
  mutate(k = 87,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)


##############################################################################
# Pmass - plant functional type
##############################################################################

##############
# N addition
##############
nadd_pmass_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "leaf_p_mass" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_pmass_photo <- data.frame(trait = "pmass", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_pmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_pmass_pft))[2,3],
                               p = coef(summary(nadd_pmass_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_pmass_myc <- data.frame(trait = "pmass", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_pmass_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_pmass_pft))[3,3],
                             p = coef(summary(nadd_pmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_pmass_nfix <- data.frame(trait = "pmass", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_pmass_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_pmass_pft))[4,3],
                              p = coef(summary(nadd_pmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_pmass_pft_results <- nadd_pmass_photo %>% 
  rbind(nadd_pmass_myc) %>% 
  rbind(nadd_pmass_nfix)


##############
# P addition
##############
padd_pmass_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "leaf_p_mass" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_pmass_photo <- data.frame(trait = "pmass", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_pmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_pmass_pft))[2,3],
                               p = coef(summary(padd_pmass_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_pmass_myc <- data.frame(trait = "pmass", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_pmass_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_pmass_pft))[3,3],
                             p = coef(summary(padd_pmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_pmass_nfix <- data.frame(trait = "pmass", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_pmass_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_pmass_pft))[4,3],
                              p = coef(summary(padd_pmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_pmass_pft_results <- padd_pmass_photo %>% 
  rbind(padd_pmass_myc) %>% 
  rbind(padd_pmass_nfix)

################
# N+P addition
################
npadd_pmass_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_pmass_photo <- data.frame(trait = "pmass", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_pmass_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_pmass_pft))[2,3],
                                p = coef(summary(npadd_pmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_pmass_myc <- data.frame(trait = "pmass", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_pmass_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_pmass_pft))[3,3],
                              p = coef(summary(npadd_pmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_pmass_nfix <- data.frame(trait = "pmass", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_pmass_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_pmass_pft))[4,3],
                               p = coef(summary(npadd_pmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_pmass_pft_results <- npadd_pmass_photo %>% 
  rbind(npadd_pmass_myc) %>% 
  rbind(npadd_pmass_nfix)

#############
# Merge Pmass moderator results, with some light cleaning
#############
pmass_pft_summary <- rbind(nadd_pmass_pft_results, 
                           padd_pmass_pft_results, 
                           npadd_pmass_pft_results) %>%
  mutate(k = 130,
         estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

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
# Merge PFT summary data frames, write to .csv
##############################################################################
marea_pft_summary %>%
  rbind(nmass_pft_summary) %>%
  rbind(narea_pft_summary) %>%
  rbind(pmass_pft_summary) %>%
  rbind(parea_pft_summary) %>%
  rbind(leafnp_pft_summary) %>%
  rbind(asat_pft_summary) %>%
  rbind(gsw_pft_summary) %>%
  rbind(rd_pft_summary) %>%
  rbind(vcmax_pft_summary) %>%
  rbind(jmax_pft_summary) %>%
  rbind(jmaxvcmax_pft_summary) %>%
  rbind(pnue_pft_summary) %>%
  rbind(ppue_pft_summary) %>%
  write.csv("../data/CNPmeta_pft_moderators.csv", 
            row.names = F)







