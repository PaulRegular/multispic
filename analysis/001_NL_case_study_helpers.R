

library(units)
library(plotly)
library(TMB)
library(multispic)
library(dplyr)
library(zoo)

options(dplyr.summarise.inform = FALSE)

## Helper function for subsetting the NL case study data and defining generic priors

nl_inputs_and_priors <- function(region = "2J3K", species = NULL) {

    ## Subset data ---------------------------------------------------------------------------------

    index <- multispic::index
    landings <- multispic::landings

    index <- index[index$region == region, ]
    landings <- landings[landings$region == region, ]

    ## Limit analysis to start of survey series and to species with indices
    sp_region <- table(index$species, index$region) > 0 # present-absent
    sp_ind <- rowSums(sp_region) == max(rowSums(sp_region))
    sub_sp <- rownames(sp_region)[sp_ind]

    # sub_sp <- c("American Plaice", "Atlantic Cod", "Greenland Halibut", "Redfish spp.")

    if (!is.null(species)) {
        sub_sp <- species
    }

    index <- index[index$species %in% sub_sp, ]
    landings <- landings[landings$species %in% sub_sp, ]

    ## Species survey = species-region-season-gear
    index$species_survey <- paste0(index$species, "-", index$region, "-", index$season, "-", index$gear)

    ## Species (stock) = species-region
    index$species <- paste0(index$species, "-", index$region)
    landings$species <- paste0(landings$species, "-", landings$region)

    ## set-up priors -------------------------------------------------------------------------------

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

    ## Relax lower range if Yankee data are included
    if (region == "2J3K") lower_log_q <- log(0.2) else lower_log_q <- log(0.1)
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

    mean_pe_betas <- 0
    sd_pe_betas <- 10

    mean_K_betas <- 0
    sd_K_betas <- 10

    ## Use full landings time series to inform prior for K but limit analysis to
    ## index time series as there is nothing to inform the process errors prior to the start
    ## of the surveys
    landings <- landings[landings$year >= min(index$year) &
                             landings$year <= max(index$year), ]

    inputs <- list(landings = landings, index = index)

    mget(c("inputs",
           "mean_log_K", "sd_log_K",
           "mean_log_B0", "sd_log_B0",
           "mean_log_r", "sd_log_r",
           "mean_log_sd_B", "sd_log_sd_B",
           "mean_log_q", "sd_log_q",
           "mean_log_sd_I", "sd_log_sd_I",
           "mean_logit_rho", "sd_logit_rho",
           "mean_logit_phi", "sd_logit_phi",
           "mean_pe_betas", "sd_pe_betas",
           "mean_K_betas", "sd_K_betas"))
}

