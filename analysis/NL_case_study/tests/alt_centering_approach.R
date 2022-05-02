

## Idea was to center K and log_B by species in template such that start par of 0 would be sensible.
## Turns out it isn't helping (assuming it is implemented correctly).
## Maybe the data should be scaled by species? But then what to do with K?

source("analysis/NL_case_study/001_helpers.R")

r <- "2J3K"

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

m <- multispic(inputs, species_cor = "none", temporal_cor = "none",
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
               pe_covariates = ~0, K_covariates = ~0, silent = FALSE, nlminb_loops = 2)



list2env(nl_inputs_and_priors(region = r, species = spp), envir = globalenv())

m <- multispic(inputs, species_cor = "one", temporal_cor = "ar1",
               n_forecast = 1, K_groups = ~1, survey_groups = ~species_survey,
               pe_covariates = ~0, K_covariates = ~0, silent = FALSE, nlminb_loops = 2)


m <- multispic(inputs, species_cor = "one", temporal_cor = "ar1",
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
                  pe_covariates = ~0, K_covariates = ~0, silent = FALSE, nlminb_loops = 2)



