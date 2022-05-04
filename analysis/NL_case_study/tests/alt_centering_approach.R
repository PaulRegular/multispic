

## fix run_loo, run_retro, and vis_multispic functions given changes to core function
## One loop for analysis?

source("analysis/NL_case_study/001_helpers.R")

nlminb_loops <- 0
silent <- TRUE
r <- "3Ps"

## Limit to the most commonly caught species by region
index_spp <- multispic::index |>
    filter(region == r) |>
    pull(species) |>
    unique()
top_spp <- multispic::landings %>%
    filter(region == r, species %in% index_spp) |>
    group_by(species) |>
    summarise(total_landings = sum(landings)) |>
    arrange(-total_landings)

## Limit to top 7 (all species with data for 2J3K)
n_spp <- switch(r, `2J3K` = 7, `3LNO` = 7, `3Ps` = 7)
spp <- head(top_spp$species, n_spp)

list2env(nl_inputs_and_priors(region = r, species = spp, K_groups = ~species), envir = globalenv())

null <- multispic(inputs, species_cor = "none", temporal_cor = "none",
                  log_K_option = par_option(option = "prior",
                                            mean = mean_log_K, sd = sd_log_K),
                  log_B0_option = par_option(option = "prior",
                                             mean = mean_log_B0, sd = sd_log_B0),
                  log_r_option = par_option(option = "prior",
                                            mean = mean_log_r, sd = sd_log_r),
                  log_sd_B_option = par_option(option = "prior",
                                               mean = mean_log_sd_B, sd = sd_log_sd_B),
                  log_q_option = par_option(option = "prior",
                                            mean = mean_log_q, sd = sd_log_q),
                  log_sd_I_option = par_option(option = "prior",
                                               mean = mean_log_sd_I, sd = sd_log_sd_I),
                  logit_rho_option = par_option(option = "prior",
                                                mean = mean_logit_rho, sd = sd_logit_rho),
                  logit_phi_option = par_option(option = "prior",
                                                mean = mean_logit_phi, sd = sd_logit_phi),
                  K_betas_option = par_option(option = "prior",
                                              mean = mean_K_betas, sd = sd_K_betas),
                  pe_betas_option = par_option(option = "prior",
                                               mean = mean_pe_betas, sd = sd_pe_betas),
                  n_forecast = 1, K_groups = ~species, survey_groups = ~species_survey,
                  pe_covariates = ~0, K_covariates = ~0, silent = silent, nlminb_loops = nlminb_loops)

list2env(nl_inputs_and_priors(region = r, species = spp), envir = globalenv())

## Hypothesis: energy flow was affected by the 1991 shift, process error is affected
##             by climate, and residual variation is temporally correlated with
##             unstructured species by species correlations.
full <- multispic(inputs, species_cor = "all", temporal_cor = "ar1",
                  log_K_option = par_option(option = "prior",
                                            mean = mean_log_K, sd = sd_log_K),
                  log_B0_option = par_option(option = "prior",
                                             mean = mean_log_B0, sd = sd_log_B0),
                  log_r_option = par_option(option = "prior",
                                            mean = mean_log_r, sd = sd_log_r),
                  log_sd_B_option = par_option(option = "prior",
                                               mean = mean_log_sd_B, sd = sd_log_sd_B),
                  log_q_option = par_option(option = "prior",
                                            mean = mean_log_q, sd = sd_log_q),
                  log_sd_I_option = par_option(option = "prior",
                                               mean = mean_log_sd_I, sd = sd_log_sd_I),
                  logit_rho_option = par_option(option = "prior",
                                                mean = mean_logit_rho, sd = sd_logit_rho),
                  logit_phi_option = par_option(option = "prior",
                                                mean = mean_logit_phi, sd = sd_logit_phi),
                  K_betas_option = par_option(option = "prior",
                                              mean = mean_K_betas, sd = sd_K_betas),
                  pe_betas_option = par_option(option = "prior",
                                               mean = mean_pe_betas, sd = sd_pe_betas),
                  n_forecast = 1, K_groups = ~1, survey_groups = ~species_survey,
                  pe_covariates = ~nlci, K_covariates = ~shift, silent = silent, nlminb_loops = nlminb_loops)
# full$sd_rep
# vis_multispic(full)

## Hypothesis: energy flow was affected by the 1991 shift, process error is affected
##             by climate, and residual variation is simple noise. (i.e., temporal
##             and species correlations are explained by the shift and climate)
just_covar <- update(full, species_cor = "none", temporal_cor = "none")
# just_covar$sd_rep

## Hypothesis: energy flow was affected by the 1991 shift and process error is
##             simple noise. (i.e., temporal and species correlations are explained by
##             the shift)
just_shift <- update(just_covar, pe_covariates = ~0)
# just_shift$sd_rep

## Hypothesis: process errors are affected by climate. (i.e., temporal and species
##             correlations are explained by climate)
just_nlci <- update(just_covar, K_covariates = ~0)
# just_nlci$sd_rep

## Hypothesis: population dynamics are affected by a common carrying capacity
##             and process error is temporally correlated with unstructured
##             species by species correlations.
just_cor <- update(full, pe_covariates = ~0, K_covariates = ~0)
# just_cor$sd_rep
# vis_multispic(just_cor)

## Hypothesis: population dynamics are affected by a common carrying capacity
##             and correlation in residual variation is shared across species
##             and over time. (i.e., there is a common but unknown environmental
##             process affecting all species)
shared_cor <- update(just_cor, species_cor = "one", temporal_cor = "ar1")
# shared_cor$sd_rep
# vis_multispic(shared_cor)

## Hypothesis: population dynamics are affected by a common carrying capacity
##             and correlation in residual variation is shared across species.
##             (i.e., there is a common but unknown environmental process affecting
##             all species and the process is noisy with no temporal dependence)
just_species_cor <- update(shared_cor, temporal_cor = "none")
# just_species_cor$sd_rep

## Hypothesis: population dynamics are affected by a common carrying capacity
##             and correlation in residual variation is temporally correlated.
##             (i.e., environmental processes affect each species differently
##             but there are carry-over effects from one year to the next)
just_temporal_cor <- update(shared_cor, species_cor = "none")
# just_temporal_cor$sd_rep

## Hypothesis: population dynamics are affected by a common carrying capacity
##             and process errors are independent across time and species
no_cor <- update(shared_cor, species_cor = "none", temporal_cor = "none")
# no_cor$sd_rep

