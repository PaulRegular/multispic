
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


library(units)
library(plotly)
library(TMB)
library(multispic)
library(dplyr)
library(zoo)


#' Wrapper of multispic for NL case study
#'
#' @details Helper function for subsetting the NL case study data, defining generic priors, and fitting the multispic model.
#'
#' @param region  Subset to this region
#' @param species Subset to these species. All species with index and landings values will be
#'                kept if `NULL`
#'

nl_multispic <- function(region = "2J3K", species = NULL) {

    ## Subset data ---------------------------------------------------------------------------------

    index <- multispic::index
    landings <- multispic::landings
    covariates <- multispic::covariates
    landings <- merge(landings, covariates, by = "year", all.x = TRUE)

    index <- index[index$region == region, ]
    landings <- landings[landings$region == region, ]

    ## Limit analysis to start of survey series and to species with indices
    sp_region <- table(index$species, index$region) > 0 # present-absent
    sp_ind <- rowSums(sp_region) == max(rowSums(sp_region))
    sub_sp <- rownames(sp_region)[sp_ind]

    sub_sp <- sub_sp[sub_sp != "Silver Hake"] # Causing convergence issues...maybe because of noise early in the series
    # sub_sp <- c("American Plaice", "Atlantic Cod", "Greenland Halibut", "Redfish spp.")

    if (!is.null(species)) {
        sub_sp <- species
    }

    index <- index[index$species %in% sub_sp, ]
    landings <- landings[landings$species %in% sub_sp, ]

    ## Survey = species-region-gear-season (i.e. groups for catchability estimates)
    index$survey <- paste0(index$species, "-", index$region, "-", index$season, "-", index$gear)

    ## Species (stock) = species-region
    index$species <- paste0(index$species, "-", index$region)
    landings$species <- paste0(landings$species, "-", landings$region)

    ## Run model -----------------------------------------------------------------------------------

    ## Set-up prior settings
    ## Note: during testing the fixed, random, and uniform_prior options rarely converged

    ## Find the smallest max aggregate landings among the groups to inform lower range for K

    tot_landings <- landings %>%
        group_by(year, region) %>%
        summarise(tot_landings = sum(landings))

    max_tot_landings <- tot_landings %>%
        group_by(region) %>%
        summarise(max = max(tot_landings)) %>%
        with(min(max))

    lower_log_r <- log(0.01)
    upper_log_r <- log(1)
    mean_log_r <- (lower_log_r + upper_log_r) / 2
    sd_log_r <- (upper_log_r - lower_log_r) / 2

    lower_log_K <- log(max_tot_landings) - upper_log_r
    upper_log_K <- log(max_tot_landings * 100) - lower_log_r
    mean_log_K <- (lower_log_K + upper_log_K) / 2
    sd_log_K <- (upper_log_K - lower_log_K) / 2

    lower_log_B0 <- log(0.01 * exp(lower_log_K) / length(unique(landings$species)))
    upper_log_B0 <- upper_log_K
    mean_log_B0 <- (lower_log_B0 + upper_log_B0) / 2
    sd_log_B0 <- c(upper_log_B0 - lower_log_B0) / 2

    lower_log_sd_B <- log(0.01)
    upper_log_sd_B <- log(1)
    mean_log_sd_B <- (lower_log_sd_B + upper_log_sd_B) / 2
    sd_log_sd_B <- (upper_log_sd_B - lower_log_sd_B) / 2

    ## Use design-based estimates of cv to inform prior for observation error
    mean_log_sd_I <- mean(log(index$cv))
    sd_log_sd_I <- sd(log(index$cv))

    lower_log_q <- log(0.2)
    upper_log_q <- log(1.2)
    mean_log_q <- (lower_log_q + upper_log_q) / 2
    sd_log_q <- (upper_log_q - lower_log_q) / 2

    lower_logit_rho <- logit(-0.9, shift = TRUE)
    upper_logit_rho <- logit(0.9, shift = TRUE)
    mean_logit_rho <- (lower_logit_rho + upper_logit_rho) / 2
    sd_logit_rho <- (upper_logit_rho - lower_logit_rho) / 2

    lower_logit_phi <- logit(0.1)
    upper_logit_phi <- logit(0.9)
    mean_logit_phi <- (lower_logit_phi + upper_logit_phi) / 2
    sd_logit_phi <- (upper_logit_phi - lower_logit_phi) / 2


    ## Use full landings time series to inform prior for K but limit analysis to
    ## index time series as there is nothing to inform the process errors prior to the start
    ## of the surveys
    landings <- landings[landings$year >= min(index$year) &
                             landings$year <= max(index$year), ]

    inputs <- list(landings = landings, index = index)


    ## Multivariate AR1 process now working
    ## Forcing the RW structure results in unusual process errors for some species
    fit <- multispic(inputs, species_cor = "all", temporal_cor = "ar1",
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
                     n_forecast = 1, K_groups = NULL, pe_covariates = ~winter_nao)

    fit

}

