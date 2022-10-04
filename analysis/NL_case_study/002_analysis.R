
## TODO:
## - Calculate one-step ahead residuals

## - Add data from 2J3K and 3Ps eco-regions <-- done
## - Fit to each region (single-species, 5 species, 10 species, all species) <-- limited by number of species in 2J3K. Kept it simple, focusing on 7 species in each region
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
## - consider dropping approximate uniform prior option

## Dev notes:
## - Forcing the RW structure results in unusual process errors for some species
## - During testing the fixed, random, and uniform_prior options rarely converged
## - Ran into convergence issues when trying to use cumsum_nlci with 3Ps

options(warn = 1)

library(progress) ## need to load progress to get progress handler to work...not sure why?
library(future)
plan(multisession, workers = 6)

source("analysis/NL_case_study/001_helpers.R")

multispic::index |> filter(gear == "Campelen") |> pull(year) |> unique() |> length()


for (r in c("2J3K", "3LNO", "3Ps")) {

    message("\n", r)

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

    ## Run single-species analysis first = NULL model
    ## Hypothesis: population dynamics are independent and governed by species-specific K
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
                      pe_covariates = ~0, K_covariates = ~0, silent = TRUE, nlminb_loops = 0)
    # null$sd_rep
    # vis_multispic(null)


    list2env(nl_inputs_and_priors(region = r, species = spp), envir = globalenv())

    ## Hypothesis: energy flow was affected by the 1991 shift and residual variation is
    ##             temporally correlated with unstructured species by species correlations.
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
                      pe_covariates = ~0, K_covariates = ~shift, silent = TRUE, nlminb_loops = 0)
    # full$sd_rep
    # vis_multispic(full)

    ## Hypothesis: energy flow was affected by the 1991 shift and process error is
    ##             simple noise. (i.e., temporal and species correlations are explained by
    ##             the shift)
    just_shift <- update(full, species_cor = "none", temporal_cor = "none")
    # just_shift$sd_rep

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
    just_temporal_cor <- update(shared_cor, species_cor = "none",
                                start_par = if(r == "3Ps") as.list(shared_cor$sd_rep, "Est") else NULL)
    # just_temporal_cor$sd_rep

    ## Hypothesis: population dynamics are affected by a common carrying capacity
    ##             and process errors are independent across time and species
    no_cor <- update(shared_cor, species_cor = "none", temporal_cor = "none",
                     start_par = if(r == "3Ps") as.list(shared_cor$sd_rep, "Est") else NULL)
    # no_cor$sd_rep

    null$loo <- run_loo(null)
    full$loo <- run_loo(full)
    just_shift$loo <- run_loo(just_shift)
    just_cor$loo  <- run_loo(just_cor)
    shared_cor$loo  <- run_loo(shared_cor)
    just_species_cor$loo  <- run_loo(just_species_cor)
    just_temporal_cor$loo  <- run_loo(just_temporal_cor)
    no_cor$loo  <- run_loo(no_cor)

    null$retro <- run_retro(null, folds = 20)
    full$retro <- run_retro(full, folds = 20)
    just_shift$retro <- run_retro(just_shift, folds = 20)
    just_cor$retro  <- run_retro(just_cor, folds = 20)
    shared_cor$retro  <- run_retro(shared_cor, folds = 20)
    just_species_cor$retro  <- run_retro(just_species_cor, folds = 20)
    just_temporal_cor$retro  <- run_retro(just_temporal_cor, folds = 20)
    no_cor$retro  <- run_retro(no_cor, folds = 20)

    fits <- mget(c("full", "just_shift", "just_cor", "shared_cor", "just_species_cor",
                   "just_temporal_cor", "no_cor", "null"))

    saveRDS(fits, file = paste0("analysis/NL_case_study/exports/fits_", r, ".rds"))

}

