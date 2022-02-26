
## TODO:
## - Calculate one-step ahead residuals

## - Add data from 2J3K and 3Ps eco-regions <-- done
## - Fit to each region (single-species, 5 species, 10 species, all species)
## - Try to combine regions (K_group = ~region) <-- possible, but too many correlation parameters, many of which may not make sense to estimate (e.g. cod in 3Ps vs. redfish in 2J3K)
## - Calculate species-specific leave one out scores to assess predictive ability of each model
## - Hypothesis: multispecies >> single-species inference

## - Write ecological oriented paper; aim for PNAS or Fish and Fisheries

## - Introduce break point for the K @ 1991 <-- done
## - Add NL climate index <-- done
## - Add capelin?? <-- limited by spatial and temporal coverage of survey
## - Fix Campelen index to 1 for cod? <-- would not converge, but did improve priors for catchability

## Consider:
## - make a tidy_par function and name par using factor levels
## - consider imposing a mean change in the collapse era (current fit is using K to cause the collapse)
## - consider dropping approximate uniform prior option

## Dev notes:
## - Forcing the RW structure results in unusual process errors for some species
## - During testing the fixed, random, and uniform_prior options rarely converged


library(progress) ## need to load progress to get progress handler to work...not sure why?
library(future)
plan(multisession, workers = 6)

source("analysis/001_NL_case_study_helpers.R")


## Species-specific analyses -----------------------------------------------------------------------

## Run single-species analysis first = NULL model
## Hypothesis: population dynamics are independent and governed by species-specific K

for (r in c("2J3K", "3LNO", "3Ps")) {

    ## Limit to top 7 most commonly caught species by region
    ## (exception: cod and redfish in 3Ps as there were convergence issues)
    top_spp <- multispic::landings %>%
        filter(region == r, !grepl("Grenadier", species),
               !(region == "3Ps" & species %in% c("Atlantic Cod", "Redfish spp."))) %>%
        group_by(species) %>%
        summarise(total_landings = sum(landings)) %>%
        arrange(-total_landings)

    spp <- head(top_spp$species, 7)

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
                     pe_covariates = ~0, K_covariates = ~0, silent = TRUE)

    null$loo <- run_loo(null)
    null$retro <- run_retro(null, folds = 15)

    saveRDS(null, file = paste0("analysis/exports/sp_fits_", r, ".rds"))

}





## Multispecies analysis ---------------------------------------------------------------------------

for (r in c("2J3K", "3LNO", "3Ps")) {

    message("\n", r)

    null <- readRDS(paste0("analysis/exports/sp_fits_", r, ".rds"))
    spp <- levels(null$landings$species)
    spp <- gsub(paste0("-", r), "", spp)

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
                      pe_covariates = ~nlci, K_covariates = ~shift, silent = TRUE)

    ## Hypothesis: energy flow was affected by the 1991 shift, process error is affected
    ##             by climate, and residual variation is simple noise. (i.e., temporal
    ##             and species correlations are explained by the shift and climate)
    just_covar <- update(full, species_cor = "none", temporal_cor = "none")

    ## Hypothesis: energy flow was affected by the 1991 shift and process error is
    ##             simple noise. (i.e., temporal and species correlations are explained by
    ##             the shift)
    just_shift <- update(just_covar, pe_covariates = ~0)

    ## Hypothesis: process errors are affected by climate. (i.e., temporal and species
    ##             correlations are explained by climate)
    just_nlci <- update(just_covar, K_covariates = ~0)

    ## Hypothesis: population dynamics are affected by a common carrying capacity
    ##             and process error is temporally correlated with unstructured
    ##             species by species correlations.
    just_cor <- update(full, pe_covariates = ~0, K_covariates = ~0)

    ## Hypothesis: population dynamics are affected by a common carrying capacity
    ##             and correlation in residual variation is shared across species
    ##             and over time. (i.e., there is a common but unknown environmental
    ##             process affecting all species)
    shared_cor <- update(just_cor, species_cor = "one", temporal_cor = "ar1")

    ## Hypothesis: population dynamics are affected by a common carrying capacity
    ##             and correlation in residual variation is shared across species.
    ##             (i.e., there is a common but unknown environmental process affecting
    ##             all species and the process is noisy with no temporal dependence)
    just_species_cor <- update(shared_cor, temporal_cor = "none")

    ## Hypothesis: population dynamics are affected by a common carrying capacity
    ##             and correlation in residual variation is temporally correlated.
    ##             (i.e., environmental processes affect each species differently
    ##             but there are carry-over effects from one year to the next)
    just_temporal_cor <- update(shared_cor, species_cor = "none")

    ## Hypothesis: population dynamics are affected by a common carrying capacity
    ##             and process errors are independent across time and species
    no_cor <- update(shared_cor, species_cor = "none", temporal_cor = "none")

    full$loo <- run_loo(full)
    just_covar$loo<- run_loo(just_covar)
    just_shift$loo <- run_loo(just_shift)
    just_nlci$loo  <- run_loo(just_nlci)
    just_cor$loo  <- run_loo(just_cor)
    shared_cor$loo  <- run_loo(shared_cor)
    just_species_cor$loo  <- run_loo(just_species_cor)
    just_temporal_cor$loo  <- run_loo(just_temporal_cor)
    no_cor$loo  <- run_loo(no_cor)

    full$retro <- run_retro(full, folds = 15)
    just_covar$retro<- run_retro(just_covar, folds = 15)
    just_shift$retro <- run_retro(just_shift, folds = 15)
    just_nlci$retro  <- run_retro(just_nlci, folds = 15)
    just_cor$retro  <- run_retro(just_cor, folds = 15)
    shared_cor$retro  <- run_retro(shared_cor, folds = 15)
    just_species_cor$retro  <- run_retro(just_species_cor, folds = 15)
    just_temporal_cor$retro  <- run_retro(just_temporal_cor, folds = 15)
    no_cor$retro  <- run_retro(no_cor, folds = 15)

    fits <- mget(c("full", "just_covar", "just_shift", "just_nlci", "just_cor",
                   "shared_cor", "just_species_cor", "just_temporal_cor", "no_cor"))

    saveRDS(fits, file = paste0("analysis/exports/spp_fits_", r, ".rds")) # spp = multiple species

}

