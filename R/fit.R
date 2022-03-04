

#' Helper for defining the settings of a parameter
#'
#' @param option    Should the parameter be estimated freely (`"fixed"`), coupled across groups
#'                  (`"coupled"`), revolve around an estimated mean and sd (`"random"`),
#'                  or a fixed mean and sd (`"prior"`). The `"coupled"` and `"random"` options
#'                  are not applicable if only one parameter is estimated (e.g. `log_K`).
#' @param mean      Mean value of the parameter. Treated as a starting value if
#'                  option is `"random"` or as a prior if option is option is `"prior"`.
#'                  Ignored if option is `"fixed"`, or `"coupled"`.
#' @param sd        SD value of the parameter. Treated as a starting value if
#'                  option is `"random"` or a prior if option is `"prior"`. Ignored
#'                  if option is `"fixed"`, or `"coupled"`.
#'
#' @details `mean` and `sd`inputs can be a vector of equal length as the number
#'           of parameter values if priors are being defined.
#'
#' @export
#'

par_option <- function(option = "fixed", mean = 0, sd = 1) {

    ## When used inside the multispic function, key lists for TMB are edited
    function (env, par_name, is_option = c("fixed", "coupled", "random", "prior"),
              center = FALSE) {

        par_length <- length(env$par[[par_name]])

        if (option == "random" && (length(mean) > 1 | length(sd) > 1)) {
            stop("Only one mean or sd starting value is expected using the 'random' option")
        }

        if (!option %in% is_option) {
            stop(paste0("The '", option, "' option is not applicable for parametere '", par_name, "'."))
        }

        options <- factor(option, levels = c("fixed", "coupled", "random", "prior"))
        env$dat[[paste0(par_name, "_option")]] <- as.numeric(options) - 1

        if (center) mean <- mean - env$log_center

        if (length(mean) == 1) {
            env$par[[paste0("mean_", par_name)]] <- rep(mean, par_length)
        } else {
            if (length(mean) == par_length) {
                env$par[[paste0("mean_", par_name)]] <- mean
            } else {
                stop("The number of mean values supplied to par_option for parameter '", par_name, "' does not equal the number of '", par_name, "' parameters.")
            }
        }

        if (length(sd) == 1) {
            env$par[[paste0("log_sd_", par_name)]] <- rep(log(sd), par_length)
        } else {
            if (length(sd) == par_length) {
                env$par[[paste0("log_sd_", par_name)]] <- log(sd)
            } else {
                stop("The number of sd values supplied to par_option for parameter '", par_name, "' does not equal the number of '", par_name, "' parameters.")
            }
        }

        if (option == "coupled") {
            env$map[[par_name]] <- factor(rep(1, par_length))
        }

        if (option == "random") {
            env$map[[paste0("mean_", par_name)]] <- factor(rep(1, par_length))
            env$map[[paste0("log_sd_", par_name)]] <- factor(rep(1, par_length))
            env$random <- c(env$random, par_name)
        } else {
            env$map[[paste0("mean_", par_name)]] <- factor(rep(NA, par_length))
            env$map[[paste0("log_sd_", par_name)]] <- factor(rep(NA, par_length))
        }

    }

}


#' Fit a multispecies surplus production model
#'
#' @param inputs             List that includes the following data.frames with required columns in
#'                           parentheses: `landings` (`species`, `year`, `landings`), `index` (`species`, `year`,
#'                           `index`). Process error covariates can be included in the `landings`
#'                           data.frame and specified using the `pe_covariates` arguments (optional). Additional
#'                           columns can also be included in the `index` data.frame to define `survey_groups`
#'                           (i.e. effects for catchability and observation error).
#' @param center             Center input values to aid convergence?
#' @param log_K_option       Settings for the estimation of `log_K`; define using [par_option()].
#' @param log_B0_option      Settings for the estimation of the starting biomass;
#'                           define using [par_option()].
#' @param log_r_option       Settings for the estimation of `log_r`; define using [par_option()].
#' @param log_sd_B_option    Settings for the estimation of sd for the process; define using
#'                           [par_option()].
#' @param log_q_option       Settings for the estimation of `log_q`; define using [par_option()].
#' @param log_sd_I_option    Settings for the estimation of sd for the indices; define using
#'                           [par_option()].
#' @param logit_rho_option   Setting for the estimation of the correlation across stocks; define using
#'                           [par_option()].
#' @param logit_phi_option   Setting for the estimation of temporal correlation in the process errors;
#'                           define using [par_option()].
#' @param K_betas_option    Setting for the estimation of covariate effects on K;
#'                           define using [par_option()].
#' @param pe_betas_option    Setting for the estimation of covariate effects on the process errors;
#'                           define using [par_option()].
#' @param species_cor        Correlation structure across species (`rho`). `"none"` will not estimate
#'                           correlations across species, `"one"` will estimate one shared correlation
#'                           parameter across species, and `"all"` will estimate correlation parameters
#'                           across all combinations of species.
#' @param temporal_cor       Correlation structure across time (`phi`). `"none"` assumes no temporal dependence
#'                           in the process errors, `"rw"` assumes a random walk, and `"ar1"` fits an
#'                           AR1 process, estimating an extra parameter.
#' @param survey_groups      Formula specifying the grouping variables to use to estimate catchability,
#'                           `log_q`, and observation error, `log_sd_I`. For example, a two-factor
#'                           main-effects model will be assumed by supplying `~survey + species`, and
#'                           interactive-effects will be assumed by supplying `~survey * species`. One
#'                           parameter will be estimated if set to `~1`.
#' @param K_groups           Formula specifying a grouping variable to use to estimate `K` and
#'                           aggregate biomass in the production equation. Biomass from all species
#'                           (stocks) will be aggregated and one `K` value estimated if set to `~1`.
#' @param K_covariates       Formula describing covariate effects on carrying capacity, K. Intercepts
#'                           are not estimated. No covariates are applied if set to `~0`.
#' @param pe_covariates      Formula describing relationship between surplus production (process error)
#'                           and covariates. Note that intercepts are not estimated. No covariates are
#'                           applied if set to `~0`.
#' @param n_forecast         Number of years to forecast. Assumes status quo landings and covariates
#'                           (i.e. terminal values assumed through projected years).
#' @param leave_out          Specific index values to leave out from the analysis (row number).
#'                           Useful for cross-validation. All data are kept if `NULL`.
#' @param start_par          List of starting parameter values. Start parameters are internally
#'                           defined, however, it may be useful to supply parameters from a previous
#'                           model fit to speed up convergence.
#' @param light              Skip running [TMB::sdreport()] and limit output to speed things up?
#' @param silent             Disable tracing information?
#'
#' @export
#'

multispic <- function(inputs,
                      center = TRUE,
                      log_K_option = par_option(),
                      log_B0_option = par_option(),
                      log_r_option = par_option(),
                      log_sd_B_option = par_option(),
                      log_q_option = par_option(),
                      log_sd_I_option = par_option(),
                      logit_rho_option = par_option(),
                      logit_phi_option = par_option(),
                      K_betas_option = par_option(),
                      pe_betas_option = par_option(),
                      species_cor = "none",
                      temporal_cor = "none",
                      survey_groups = ~survey,
                      K_groups = ~1,
                      K_covariates = ~0,
                      pe_covariates = ~0,
                      n_forecast = 0,
                      leave_out = NULL,
                      start_par = NULL,
                      light = FALSE,
                      silent = FALSE) {

    start_time <- Sys.time()
    call <- match.call()

    landings <- inputs$landings
    index <- inputs$index

    ## Stop if there are any NAs
    if (any(is.na(landings$landings)) | any(is.na(index$index))) {
        stop("NA landings or index values are not permitted in the current model.")
    }

    ## Values to keep (define here before any filtering or ordering below)
    index$left_out <- rep(FALSE, length(index$index))
    if (!is.null(leave_out)) {
        index$left_out[leave_out] <- TRUE
    }

    ## Drop zeros
    ind <- index$index == 0
    if (sum(ind) > 0) {
        warning("Index values equal to zero are present in the inputs. The current model lacks a way to deal with index values equal to zero. These values are being dropped from the analysis (i.e. treated as missing values).")
        index <- index[!ind, ]
    }

    ## Set-up inputs for forecasts
    if (n_forecast > 0) {
        terminal_landings <- landings[landings$year == max(landings$year), ]
        sq_landings <- lapply(seq(n_forecast), function(i) {
            terminal_landings$year <- terminal_landings$year + i
            terminal_landings
        })
        sq_landings <- do.call(rbind, sq_landings)
        landings <- rbind(landings, sq_landings)
    }

    ## Set-up factors for indexing in TMB
    if (!identical(sort(unique(landings$species)), sort(unique(index$species)))) {
        stop("One or more species were not found in both the landings and index data.frames. Please ensure there are landings and index data for all species.")
    }
    landings$species <- factor(landings$species)
    landings$y <- factor(landings$year)
    landings$sy <- factor(paste0(landings$species, "-", landings$year))
    landings <- landings[order(landings$sy), ]
    index$sy <- factor(paste0(index$species, "-", index$year), levels = levels(landings$sy))
    index$species <- factor(index$species, levels = levels(landings$species))

    ## Set-up model matrix | formula with covariates
    if (pe_covariates == ~0) {
        pe_model_mat <- matrix(rep(0, nrow(landings)), ncol = 1)
    } else {
        pe_model_mat <- model.matrix(pe_covariates, data = landings)[, -1, drop = FALSE] # drop intercept because that is defined by the process errors
    }
    if (K_covariates == ~0) {
        K_model_mat <- matrix(rep(0, nrow(landings)), ncol = 1)
    } else {
        K_model_mat <- model.matrix(K_covariates, data = landings)[, -1, drop = FALSE] # drop intercept because that is defined by K
    }
    if (K_groups == ~1) {
        K_map <- rep(0, nlevels(landings$species))
        B_group_mat <- matrix(1, nrow = nlevels(landings$species), ncol = nlevels(landings$species))
    } else {

        if (length(all.vars(K_groups)) > 1) {
            stop("Please supply only one variable to the K_groups argument; mapping by more than one variable has yet to be implemented.")
        }

        sp_group <- unique(landings[, c("species", all.vars(K_groups))])
        sp_group <- sp_group[order(sp_group$species), ]
        sp_group[, all.vars(K_groups)] <- factor(sp_group[, all.vars(K_groups)])
        if (nrow(sp_group) > nlevels(landings$species)) {
            stop("All species are not nested within the K_groups variable. Please check groupings.")
        }

        ## Set up a matrix of 0s and 1s to control the summation of biomass
        f <- as.formula(paste(Reduce(paste, deparse(K_groups)), "-1"))
        mm <- model.matrix(f, data = sp_group)
        colnames(mm) <- levels(sp_group[, all.vars(K_groups)])
        rownames(mm) <- levels(sp_group$species)
        B_group_mat <- mm[, sp_group[, all.vars(K_groups)]]

        ## And set up a map for the K parameters
        K_map <- as.numeric(factor(sp_group[[2]])) - 1

    }

    ## Arrange data by survey groups and identify unique surveys
    unique_surveys <- index[do.call(order, index[, all.vars(survey_groups), drop = FALSE]), ]
    unique_surveys <- unique(unique_surveys[, all.vars(survey_groups), drop = FALSE])
    unique_surveys$survey_id <- seq(nrow(unique_surveys)) - 1
    unique_surveys$survey <- do.call(paste, c(unique_surveys[, all.vars(survey_groups), drop = FALSE],
                                              sep = "-"))
    unique_surveys$survey <- factor(unique_surveys$survey, levels = unique_surveys$survey)
    survey_model_mat <- model.matrix(survey_groups, data = unique_surveys)
    index <- merge(index, unique_surveys, by = all.vars(survey_groups))
    index$survey <- factor(index$survey, levels = unique_surveys$survey) # ensure survey is a factor

    ## Center index and landings by the mean(log(index) ~ K_group) to aid convergence
    by_sp_yr_grp <- paste(c("species", "year", all.vars(K_groups)), collapse = " + ")
    by_yr_grp <- paste(c("year", all.vars(K_groups)), collapse = " + ")
    by_grp <- ifelse(K_groups == ~1, "0", paste(all.vars(K_groups), collapse = " + "))
    mean_index <- aggregate(as.formula(paste("index ~ ", by_sp_yr_grp)),
                            FUN = mean, data = index)
    tot_mean_index <- aggregate(as.formula(paste("index ~ ", by_yr_grp)),
                                FUN = sum, data = mean_index)
    log_center <- ifelse(center, mean(log(tot_mean_index$index)), 0)
    index$index <- exp(log(index$index) - log_center)
    landings$landings <- exp(log(landings$landings) - log_center)

    ## Compute total landings by year to inform starting value for K
    ## And mean log index to inform starting values for B
    tot_landings <- aggregate(as.formula(paste("landings ~ ", by_yr_grp)),
                              FUN = sum, data = landings)
    max_landings <- aggregate(as.formula(paste("landings ~", by_grp)),
                              FUN = max, data = tot_landings)$landings
    mean_log_index <- aggregate(index ~ species, FUN = function(x) mean(log(x)), data = index)

    ## Set-up the objects for TMB
    if (nlevels(landings$species) == 1) {
        n_rho <- 1
        if (species_cor != "none") {
            warning("Data from only one species was supplied; species_cor is therefore moot; setting to 'none'")
            species_cor <- "none"
        }
    } else {
        n_rho <- sum(lower.tri(matrix(NA, nrow = nlevels(landings$species),
                                      ncol = nlevels(landings$species))))
    }

    dat <- list(L = as.numeric(landings$landings),
                L_species = as.numeric(landings$species) - 1,
                L_year = as.numeric(landings$y) - 1,
                I = as.numeric(index$index),
                I_species = as.numeric(index$species) - 1,
                I_survey = index$survey_id,
                I_sy = as.numeric(index$sy) - 1,
                min_B = 0.00001,
                nY = nlevels(landings$y),
                nS = nlevels(landings$species),
                survey_covariates = survey_model_mat,
                K_covariates = K_model_mat,
                pe_covariates = pe_model_mat,
                K_map = K_map,
                B_groups = B_group_mat,
                keep = as.numeric(!index$left_out))

    par <- list(log_B = matrix(ceiling(mean_log_index$index),
                               ncol = dat$nS, nrow = dat$nY, byrow = TRUE),
                log_sd_B = rep(-1, nlevels(landings$species)),
                logit_rho = rep(0, n_rho),
                logit_phi = 0,
                log_K = ceiling(log(max_landings)),
                log_B0 = ceiling(mean_log_index$index),
                log_r = rep(-1, nlevels(landings$species)),
                log_m = rep(log(2), nlevels(landings$species)),
                log_q = rep(-1, nrow(unique_surveys)),
                log_q_betas = rep(0, ncol(survey_model_mat)),
                log_sd_I = rep(-1, nrow(unique_surveys)),
                log_sd_I_betas = rep(0, ncol(survey_model_mat)),
                K_betas = rep(0, ncol(K_model_mat)),
                pe_betas =  rep(0, ncol(pe_model_mat)))

    map <- list(log_m = factor(rep(NA, nlevels(landings$species))))
    if (species_cor == "one") {
        map$logit_rho <- factor(rep(1, length(par$logit_rho)))
    }
    if (species_cor == "none") {
        map$logit_rho <- factor(rep(NA, length(par$logit_rho)))
        dat$logit_rho_option <- 0 # skip prior / random effect loop
    }
    if (temporal_cor == "none") {
        map$logit_phi <- factor(NA)
        par$logit_phi <- -10 # results in a value very close to 0
    }
    if (temporal_cor == "rw") {
        map$logit_phi <- factor(NA)
        par$logit_phi <- 10 # results in a value very close to 1
    }
    if (pe_covariates == ~0) {
        map$pe_betas <- factor(NA)
    }
    if (K_covariates == ~0) {
        map$K_betas <- factor(NA)
    }

    random <- "log_B"

    ## Augment dat, par, map, and random objects using par_option closures
    log_K_option(environment(), "log_K", center = TRUE)
    log_B0_option(environment(), "log_B0", center = TRUE)
    log_r_option(environment(), "log_r")
    log_sd_B_option(environment(), "log_sd_B")
    log_q_option(environment(), "log_q", is_option = c("fixed", "random", "prior"))
    log_sd_I_option(environment(), "log_sd_I", is_option = c("fixed", "random", "prior"))
    logit_rho_option(environment(), "logit_rho")
    logit_phi_option(environment(), "logit_phi", is_option = c("fixed", "prior"))
    K_betas_option(environment(), "K_betas", is_option = c("fixed", "prior"))
    pe_betas_option(environment(), "pe_betas", is_option = c("fixed", "prior"))

    ## Extra tweaks for q and sd_I
    if (dat$log_q_option != 2) { # != "random"
        map$log_q <- factor(rep(NA, length(par$log_q)))
    } else {
        map$log_q_betas <- factor(rep(NA, length(par$log_q_betas)))
        if (length(all.vars(survey_groups)) > 1) {
            message("Fitting all unique survey_groups as random effects (i.e., a mixture of fixed main effects and random effects are currently not possible)")
        }
    }
    if (dat$log_sd_I_option != 2) {
        map$log_sd_I <- factor(rep(NA, length(par$log_sd_I)))
    } else {
        map$log_sd_I_betas <- factor(rep(NA, length(par$log_sd_I_betas)))
        if (length(all.vars(survey_groups)) > 1) {
            message("Fitting all unique survey_groups as random effects (i.e., a mixture of fixed main-effects and random effects are currently not possible)")
        }
    }

    ## Fit model
    if (!is.null(start_par)) par <- start_par
    obj <- TMB::MakeADFun(dat, par, map = map, random = random, DLL = "multispic",
                          silent = silent, checkParameterOrder = FALSE)
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 1000, iter.max = 1000))

    ## Reset scale
    landings$landings <- exp(log(landings$landings) + log_center)
    index$index <- exp(log(index$index) + log_center)

    ## Extract REPORT objects
    rep <- obj$report()

    index$log_index <- log(index$index)
    index$log_pred_index <- rep$log_pred_I + log_center
    index$std_res <- rep$log_I_std_res

    pop <- landings
    pop$pe <- exp(rep$log_pe)
    pop$res_pe <- exp(rep$log_res_pe)
    pop$log_std_res_pe <- rep$log_std_res_pe
    pop$B <- exp(log(rep$B_vec) + log_center)
    pop$F <- rep$F
    pop$K <- exp(log(rep$K_vec) + log_center)

    tot_pop <- tot_landings
    tot_pop$landings <- exp(log(tot_pop$landings) + log_center)
    tot_pop$B <- exp(log(c(rep$tot_B)) + log_center)

    se <- sd_rep <- par_lwr <- par_upr <- NULL

    if (!light) {

        ## Extract ADREPORT objects
        sd_rep <- sdreport(obj)
        par <- as.list(sd_rep, "Est")
        se <- as.list(sd_rep, "Std. Error")
        par_lwr <- lapply(seq_along(par), function(i) par[[i]] - 1.96 * se[[i]])
        par_upr <- lapply(seq_along(par), function(i) par[[i]] + 1.96 * se[[i]])
        names(par_lwr) <- names(par_upr) <- names(par)
        par$log_B0 <- par$log_B0 + log_center
        par$log_B <- par$log_B + log_center
        par$log_K <- par$log_K + log_center
        par_lwr$log_B0 <- par_lwr$log_B0 + log_center
        par_lwr$log_B <- par_lwr$log_B + log_center
        par_lwr$log_K <- par_lwr$log_K + log_center
        par_upr$log_B0 <- par_upr$log_B0 + log_center
        par_upr$log_B <- par_upr$log_B + log_center
        par_upr$log_K <- par_upr$log_K + log_center

        ## Extract and append fits
        est <- split(unname(sd_rep$value), names(sd_rep$value))
        sd <- split(sd_rep$sd, names(sd_rep$value))
        lwr <- split(unname(sd_rep$value) - 1.96 * sd_rep$sd, names(sd_rep$val))
        upr <- split(unname(sd_rep$value) + 1.96 * sd_rep$sd, names(sd_rep$val))
        index$pred <- exp(est$log_pred_I + log_center)
        index$pred_lwr <- exp(lwr$log_pred_I + log_center)
        index$pred_upr <- exp(upr$log_pred_I + log_center)

        ## Extract population estimates
        pop$B <- exp(est$log_B_vec + log_center)
        pop$B_lwr <- exp(lwr$log_B_vec + log_center)
        pop$B_upr <- exp(upr$log_B_vec + log_center)
        pop$F <- exp(est$log_F)
        pop$F_lwr <- exp(lwr$log_F)
        pop$F_upr <- exp(upr$log_F)
        pop$K <- exp(est$log_K_vec + log_center)
        pop$K_lwr <- exp(lwr$log_K_vec + log_center)
        pop$K_upr <- exp(upr$log_K_vec + log_center)

        tot_pop$B <- exp(est$log_tot_B + log_center)
        tot_pop$B_lwr <- exp(lwr$log_tot_B + log_center)
        tot_pop$B_upr <- exp(upr$log_tot_B + log_center)

        tot_pop$K <- exp(est$log_grouped_K + log_center)
        tot_pop$K_lwr <- exp(lwr$log_grouped_K + log_center)
        tot_pop$K_upr <- exp(upr$log_grouped_K + log_center)

        ## Replace log_q and log_sd_I with reported pred values
        ## (will be same values if option is "random")
        par$log_q <- est$pred_log_q
        se$log_q <- sd$pred_log_q
        par_lwr$log_q <- lwr$pred_log_q
        par_upr$log_q <- upr$pred_log_q
        par$log_sd_I <- est$pred_log_sd_I
        se$log_sd_I <- sd$pred_log_sd_I
        par_lwr$log_sd_I <- lwr$pred_log_sd_I
        par_upr$log_sd_I <- upr$pred_log_sd_I

        ## Simplify logit_rho if coupled
        if (species_cor != "all") {
            par$mean_logit_rho <- par$mean_logit_rho[1]
            par$log_sd_logit_rho <- par$log_sd_logit_rho[1]
            par$logit_rho <- par$logit_rho[1]
            se$logit_rho <- se$logit_rho[1]
            par_lwr$logit_rho <- par_lwr$logit_rho[1]
            par_upr$logit_rho <- par_upr$logit_rho[1]
        }

    }

    ## Calculate marginal AIC
    mAIC <- 2 * length(opt$par) + 2 * opt$objective

    end_time <- Sys.time()
    run_dur <- end_time - start_time

    out <- list(call = call, run_dur = run_dur, log_center = log_center, tmb_dat = dat,
                obj = obj, opt = opt, sd_rep = sd_rep, rep = rep, par = par, se = se,
                par_lwr = par_lwr, par_upr = par_upr, index = index, landings = landings,
                pop = pop, tot_pop = tot_pop, mAIC = mAIC)

}


#' Function for running leave one out cross-validation
#'
#' @param fit         Object from [fit_model()]
#' @param progress    Display progress bar? (Generated using the progressr and progress packages)
#'
#' @return Returns a list with three objects:
#'    1) `preds`  -  log observations that were left out at each step with log predictions, and
#'    2) `mse`  -  mean squared error of the predictions (leave one out cross validation score).
#'
#' @export
#'

run_loo <- function(fit, progress = TRUE) {

    ## Consider adding previous loop approach to allow a base approach as an alternate option
    pkg <- c("furrr", "progressr", "progress")
    for (p in pkg) {
        if (!requireNamespace(p, quietly = TRUE)) {
            stop(paste(p, "is needed for run_loo to work. Please install it."), call. = FALSE)
        }
    }

    n <- length(fit$index$index)

    if (!is.null(fit$sd_rep)) {
        start_par <- as.list(fit$sd_rep, "Est")
        fit$sd_rep <- NULL # drop large object b/c it is no longer needed and does not need to be sent to workers
    } else {
        stop("Object sd_rep is NA in the supplied fit object. Please re-run model with light = FALSE.")
    }

    loo <- function(i, fit, start_par, p = NULL) {
        if (!is.null(p)) p()
        f <- try(update(fit, leave_out = i, start_par = start_par,
                        light = TRUE, silent = TRUE))
        if (class(f) == "try-catch" || f$opt$message == "false convergence (8)") {
            obs <- pred <- NA
        } else {
            obs <- f$index$log_index[f$index$left_out]
            pred <- f$index$log_pred_index[f$index$left_out]
        }
        data.frame(obs = obs, pred = pred)
    }

    progressr::with_progress({
        progressr::handlers(list(
            progressr::handler_progress(
                format = "[:bar] :percent (:current / :total) in :elapsed (eta: :eta)",
                width  = 100, clear = FALSE
            )))
        if (progress) {
            p <- progressr::progressor(steps = n)
        } else {
            p <- NULL
        }
        obs_pred <- furrr::future_map(seq(n), loo, fit = fit, start_par = start_par, p = p,
                                      .options = furrr::furrr_options(packages = "multispic"))
    })

    obs_pred <- do.call(rbind, obs_pred)
    names(obs_pred) <- c("log_index", "log_pred_index")
    preds <- cbind(fit$index[, c("year", "survey", "species")], obs_pred)

    if (any(is.na(obs_pred$pred))) {
        warning(paste("When iterating across ", n, " observations, model fitting failed for ", sum(is.na(obs_pred$pred)), " cases when an observation was left out."))
    }

    list(preds = preds, mse = mean((preds$log_index - preds$log_pred_index) ^ 2, na.rm = TRUE))

}


#' Function for running a retrospective analysis
#'
#' @details This function runs a retrospective analysis whereby terminal survey index estimates
#'          are excluded from the analysis. A hindcast is also preformed in each fold where
#'          observed values are left out but predicted to measure the models forecasting skill.
#'
#' @param fit         Object from [fit_model()]
#' @param folds       Number of years to 'fold' back
#' @param progress    Display progress bar?  (Generated using the progressr and progress packages)
#'
#' @return Returns a list with three objects:
#'    1) `retro_fits`  -  a list including a series of fits from each retrospective fold
#'    2) `hindcasts` -  a data.frame with observed, but left out, and predicted indices from each fold
#'    3) `mse`  -  mean squared error of the hindcasts (measure of forecasting skill)
#'
#' @export
#'

run_retro <- function(fit, folds, progress = TRUE) {

    ## Consider adding previous loop approach to allow a base approach as an alternate option
    pkg <- c("furrr", "progressr", "progress")
    for (p in pkg) {
        if (!requireNamespace(p, quietly = TRUE)) {
            stop(paste(p, "is needed for run_loo to work. Please install it."), call. = FALSE)
        }
    }

    if (!is.null(fit$sd_rep)) {
        start_par <- as.list(fit$sd_rep, "Est")
        fit$sd_rep <- NULL # drop large object b/c it is no longer needed and does not need to be sent to workers
    } else {
        stop("Object sd_rep is NA in the supplied fit object. Please re-run model with light = FALSE.")
    }

    retro <- function(i, fit, start_par, p = NULL) {

        if (!is.null(p)) p()

        inputs <- get(fit$call$inputs)
        index <- inputs$index
        landings <- inputs$landings

        ## Subset the input data on year at a time, keeping one extra year of indices and landings
        retro_index <- index[index$year <= retro_years[i] + 1, ]
        retro_landings <- landings[landings$year <= retro_years[i] + 1, ]
        retro_inputs <- list(landings = retro_landings, index = retro_index)

        ## Identify observations to leave out, but provide predictions for these values
        ind <- which(retro_index$year == retro_years[i] + 1)

        ## Use start_par to aid convergence and speed up loops (need to adjust log_B dimensions)
        nY <- length(unique(retro_landings$year)) + 1
        retro_start_par <- start_par
        retro_start_par$log_B <- retro_start_par$log_B[seq(nY), , drop = FALSE]

        retro_fit <- try(update(fit, inputs = retro_inputs, leave_out = ind, start_par = retro_start_par,
                                light = TRUE, silent = TRUE))

        if (class(retro_fit) == "try-catch" || retro_fit$opt$message == "false convergence (8)") {

            retro_fit <- NA
            hindcast <- NULL

        } else {

            if (sum(ind) > 0) {
                hindcast <- retro_fit$index[retro_fit$index$left_out,
                                            c("year", "survey", "species", "log_index", "log_pred_index")]
                hindcast$retro_year <- retro_years[i]
            } else {
                hindcast <- NULL
                warning(paste("A hindcast is moot for", retro_years[i], "as there is no index to predict."))
            }

        }

        list(retro_fit = retro_fit, hindcast = hindcast)

    }

    progressr::with_progress({
        progressr::handlers(list(
            progressr::handler_progress(
                format = "[:bar] :percent (:current / :total) in :elapsed (eta: :eta)",
                width  = 100, clear = FALSE
            )))
        if (progress) {
            p <- progressr::progressor(steps = folds)
        } else {
            p <- NULL
        }
        terminal_year <- max(fit$index$year)
        retro_years <- terminal_year - seq(folds)
        retro_hind <- furrr::future_map(seq(folds), retro, fit = fit, p = p, start_par = start_par,
                                        .options = furrr::furrr_options(packages = "multispic"))
    })

    retro_fits <- lapply(retro_hind, `[[`, "retro_fit")
    names(retro_fits) <- retro_years
    hindcasts <- lapply(retro_hind, `[[`, "hindcast")
    hindcasts <- do.call(rbind, hindcasts)

    if (any(is.na(retro_fits))) {
        warning(paste("While folding back ", folds, " years, model fitting failed in ", sum(is.na(retro_fits)), " cases."))
    }

    mse <- mean((hindcasts$log_index - hindcasts$log_pred_index) ^ 2)
    list(retro_fits = retro_fits, hindcasts = hindcasts, mse = mse)

}


