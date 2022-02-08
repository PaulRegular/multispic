
## TODO:
## - Calculate one-step ahead residuals

## - Add data from 2J3K and 3Ps eco-regions <-- done
## - Fit to each region (single-species, 5 species, 10 species, all species)
## - Try to combine regions (K_group = ~region) <-- possible, but too many correlation parameters, many of which may not make sense to estimate (e.g. cod in 3Ps vs. redfish in 2J3K)
## - Calculate species-specific leave one out scores to assess predictive ability of each model
## - Hypothesis: multispecies >> single-species inference

## Consider:
## - make a tidy_par function and name par using factor levels
## - consider imposing a mean change in the collapse era (current fit is using K to cause the collapse)


source("analysis/001_NL_case_study_helpers.R")


## Multispecies analysis ---------------------------------------------------------------------------

for (r in c("2J3K", "3LNO", "3Ps")) {

    list2env(nl_inputs_and_priors(region = r, species = NULL), envir = globalenv())

    if (r == "3LNO") {
        survey_formula <- ~gear + season * species
    } else {
        survey_formula <- ~gear + species
    }

    ## Dev notes:
    ## - Forcing the RW structure results in unusual process errors for some species
    ## - During testing the fixed, random, and uniform_prior options rarely converged
    full <- multispic(inputs, species_cor = "all", temporal_cor = "ar1",
                      log_K_option = par_option(option = "normal_prior",
                                                mean = mean_log_K, sd = sd_log_K),
                      log_B0_option = par_option(option = "normal_prior",
                                                 mean = mean_log_B0, sd = sd_log_B0),
                      log_r_option = par_option(option = "normal_prior",
                                                mean = mean_log_r, sd = sd_log_r),
                      log_sd_B_option = par_option(option = "normal_prior",
                                                   mean = mean_log_sd_B, sd = sd_log_sd_B),
                      log_q_option = par_option(option = "normal_prior",
                                                mean = mean_log_q, sd = sd_log_q),
                      log_sd_I_option = par_option(option = "normal_prior",
                                                   mean = mean_log_sd_I, sd = sd_log_sd_I),
                      logit_rho_option = par_option(option = "normal_prior",
                                                    mean = mean_logit_rho, sd = sd_logit_rho),
                      logit_phi_option = par_option(option = "normal_prior",
                                                    mean = mean_logit_phi, sd = sd_logit_phi),
                      n_forecast = 1, K_groups = ~1, survey_groups = survey_formula,
                      pe_covariates = ~winter_nao)


    no_nao <- update(full, pe_covariates = ~0)
    just_nao <- update(full, species_cor = "none", temporal_cor = "none")
    one_species_cor <- update(no_nao, species_cor = "one")
    no_species_cor <- update(one_species_cor, species_cor = "none")
    no_temporal_cor <- update(no_species_cor, temporal_cor = "none")

    full$loo <- run_loo(full)
    no_nao$loo<- run_loo(no_nao)
    just_nao$loo <- run_loo(just_nao)
    one_species_cor$loo  <- run_loo(one_species_cor)
    no_species_cor$loo  <- run_loo(no_species_cor)
    no_temporal_cor$loo  <- run_loo(no_temporal_cor)

    fits <- mget(c("full", "no_nao", "just_nao", "one_species_cor",
                   "no_species_cor", "no_temporal_cor"))

    saveRDS(fits, file = paste0("analysis/exports/spp_fits_", r, ".rds")) # spp = multiple species

}


## Add retrospective analysis

for (r in c("2J3K", "3LNO", "3Ps")) {

    list2env(nl_inputs_and_priors(region = r, species = NULL), envir = globalenv())

    if (r == "3LNO") {
        survey_formula <- ~gear + season * species
    } else {
        survey_formula <- ~gear + species
    }

    fits <- readRDS(paste0("analysis/exports/spp_fits_", r, ".rds"))

    for (i in seq_along(fits)) {
        fits[[i]]$retro <- run_retro(fits[[i]], folds = 10)
    }

    saveRDS(fits, file = paste0("analysis/exports/spp_fits_", r, ".rds"))

}



## Species-specific analyses -----------------------------------------------------------------------


for (r in c("2J3K", "3LNO", "3Ps")) {

    if (r == "3LNO") {
        survey_formula <- ~gear + season
    } else {
        survey_formula <- ~gear
    }

    spp_fits <- readRDS(paste0("analysis/exports/spp_fits_", r, ".rds"))
    spp_region <- levels(spp_fits$full$landings$species)
    spp <- gsub(paste0("-", r), "", spp_region)
    names(spp) <- spp_region

    ## Single-species fits = NULL model
    null <- vector("list", length(spp))
    names(null) <- spp_region

    for (sr in spp_region) {

        message("\n", sr)

        list2env(nl_inputs_and_priors(region = r, species = spp[sr]), envir = globalenv())

        fit <- try(multispic(inputs, species_cor = "none", temporal_cor = "none",
                             log_K_option = par_option(option = "normal_prior",
                                                       mean = mean_log_K, sd = sd_log_K),
                             log_B0_option = par_option(option = "normal_prior",
                                                        mean = mean_log_B0, sd = sd_log_B0),
                             log_r_option = par_option(option = "normal_prior",
                                                       mean = mean_log_r, sd = sd_log_r),
                             log_sd_B_option = par_option(option = "normal_prior",
                                                          mean = mean_log_sd_B, sd = sd_log_sd_B),
                             log_q_option = par_option(option = "normal_prior",
                                                       mean = mean_log_q, sd = sd_log_q),
                             log_sd_I_option = par_option(option = "normal_prior",
                                                          mean = mean_log_sd_I, sd = sd_log_sd_I),
                             logit_rho_option = par_option(option = "normal_prior",
                                                           mean = mean_logit_rho, sd = sd_logit_rho),
                             logit_phi_option = par_option(option = "normal_prior",
                                                           mean = mean_logit_phi, sd = sd_logit_phi),
                             n_forecast = 1, K_groups = ~1, survey_groups = survey_formula,
                             pe_covariates = ~0, silent = TRUE))

        if (class(fit) == "try-error" || fit$opt$message == "false convergence (8)") {
            null[[sr]] <- "Did not converge"
        } else {
            fit$loo <- run_loo(fit)
            fit$retro <- run_retro(fit, folds = 10)
            null[[sr]] <- fit
        }

    }

    saveRDS(null, file = paste0("analysis/exports/sp_fits_", r, ".rds"))

}


