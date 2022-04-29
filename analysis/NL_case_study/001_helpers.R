
library(units)
library(plotly)
library(TMB)
library(multispic)
library(dplyr)
library(zoo)

options(dplyr.summarise.inform = FALSE)

## Helper function for subsetting the NL case study data and defining generic priors

nl_inputs_and_priors <- function(region = "2J3K", species = NULL, K_groups = ~region) {

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
    landings <- landings[landings$species %in% sub_sp &
                             landings$year >= min(index$year) &
                             landings$year <= max(index$year), ]

    ## Species survey = species-region-season-gear
    index$species_survey <- paste0(index$species, "-", index$region, "-", index$season, "-", index$gear)

    ## Species (stock) = species-region
    index$species <- paste0(index$species, "-", index$region)
    landings$species <- paste0(landings$species, "-", landings$region)

    ## set-up priors -------------------------------------------------------------------------------

    ## Use max aggregate landings to inform prior for K,
    ## and landings in year 1 to inform prior for the biomass

    by_yr_grp <- paste(c("year", all.vars(K_groups)), collapse = " + ")
    by_grp <- ifelse(K_groups == ~1, "0", paste(all.vars(K_groups), collapse = " + "))
    tot_landings <- aggregate(as.formula(paste("landings ~", by_yr_grp)), data = landings, FUN = sum)
    max_tot_landings <- aggregate(as.formula(paste("landings ~", by_grp)),
                                  data = tot_landings, FUN = max) %>% with(landings)

    lower_log_r <- log(0.1)
    upper_log_r <- log(1)
    mean_log_r <- (lower_log_r + upper_log_r) / 2
    sd_log_r <- (upper_log_r - lower_log_r) / 2

    lower_log_K <- log(max_tot_landings) - upper_log_r
    upper_log_K <- log(max_tot_landings * 4) - lower_log_r
    mean_log_K <- (lower_log_K + upper_log_K) / 2
    sd_log_K <- (upper_log_K - lower_log_K) / 2

    L0 <- landings[landings$year == min(index$year), ]
    L0$species <- factor(L0$species)
    L0 <- L0[order(L0$species), ]
    lower_log_B0 <- log(L0$landings) - upper_log_r
    upper_log_B0 <- log(L0$landings * 4) - lower_log_r
    mean_log_B0 <- (lower_log_B0 + upper_log_B0) / 2
    sd_log_B0 <- c(upper_log_B0 - lower_log_B0) / 2

    lower_log_sd_B <- log(0.01)
    upper_log_sd_B <- log(1)
    mean_log_sd_B <- (lower_log_sd_B + upper_log_sd_B) / 2
    sd_log_sd_B <- (upper_log_sd_B - lower_log_sd_B) / 2

    ## Use design-based estimates of cv to inform prior for observation error
    ## Note: sd of cv was widened as it is an imperfect and partial indicator of observation error.
    ##       Multiplier is greater for deep water species as the survey likely captures less of their
    ##       range (i.e. greater variation is possible as the species may move in and out of the
    ##       survey area)
    deep_spp <- "Greenland Halibut|Atlantic Halibut|Witch Flounder|Redfish spp.|White Hake|Silver Hake|Monkfish"
    log_cv_stats <- aggregate(cv ~ species_survey, data = index,
                              FUN = function(x) {
                                  c(mean = mean(log(x)), sd = sd(log(x)))
                              })
    mean_log_sd_I <-  log_cv_stats$cv[, "mean"] # mean(log(index$cv))
    sd_log_sd_I <-  log_cv_stats$cv[, "sd"] # sd(log(index$cv))
    ind <- grepl(deep_spp, log_cv_stats$species_survey)
    sd_log_sd_I <- ifelse(ind, sd_log_sd_I * 4, sd_log_sd_I * 2)
    # plot_ly(y = exp(mean_log_sd_I), x = log_cv_stats$species_survey)

    ## Use survey coverage to inform lower bound for q
    ## Values are further reduced to account for selectivity and availability.
    ## Percent reduction is greater for deep water species as the survey likely captures less of
    ## their range.
    coverage_stats <- aggregate(coverage ~ species_survey, data = index,
                                FUN = function(x) {
                                    round(unique(x), 1)
                                })
    ind <- grepl(deep_spp, coverage_stats$species_survey)
    if (r == "3Ps") {
        ## Relatively small area; mixing may therefore be more prevalent; q may therefore be lower.
        lower_log_q <- ifelse(ind, log(coverage_stats$coverage * 0.1),
                              log(coverage_stats$coverage * 0.2))
    } else {
        lower_log_q <- ifelse(ind, log(coverage_stats$coverage * 0.2),
                              log(coverage_stats$coverage * 0.6))
    }
    upper_log_q <- rep(log(1), nrow(log_cv_stats))
    mean_log_q <- (lower_log_q + upper_log_q) / 2
    sd_log_q <- (upper_log_q - lower_log_q) / 2
    # plot_ly(y = exp(mean_log_q), x = coverage_stats$species_survey)

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

