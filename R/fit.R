

#' Helper for defining the settings of a parameter
#'
#' @param option    Should the parameter be estimated freely (`"fixed"`), coupled across groups
#'                  (`"coupled"`), revolve around an estimated mean and sd (`"random"`),
#'                  or a fixed mean and sd (`"normal_prior"`), or be constrained between
#'                  a lower and upper value (`"uniform_prior"`). The `"coupled"` and `"random"` options
#'                  are not applicable if only one parameter is estimated (e.g. `log_K`).
#' @param mean      Mean value of the parameter. Treated as a starting value if
#'                  option is `"random"` or as a prior if option is option is `"normal_prior"`.
#'                  Ignored if option is `"fixed"`, `"coupled"`, or `"uniform_prior"`.
#' @param sd        SD value of the parameter. Treated as a starting value is
#'                  option is `"random"` or a prior if option is "prior". Ignored
#'                  if option is `"fixed"`, `"coupled"`, or `"uniform_prior"`.
#' @param lower     Lower range for the parameter.
#' @param upper     Upper range for the parameter.
#'
#' @export
#'

par_option <- function(option = "fixed", mean = 0, sd = 1, lower = -10, upper = 10) {
    list(option = factor(option, levels = c("fixed", "coupled", "random",
                                            "normal_prior", "uniform_prior")),
         mean = mean, sd = sd, lower = lower, upper = upper)
}


#' Fit a multispecies surplus production model
#'
#' @param inputs             List that includes the following data.frames with required columns in
#'                           parentheses: `landings` (`species`, `year`, `landings`), `index` (`species`, `year`,
#'                           `survey`, `index`). Covariates can be included in the `landings` data.frame
#'                           and specified using the `covariates` arguments (optional).
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
#' @param species_cor        Correlation structure across species (`rho`). `"none"` will not estimate
#'                           correlations across species, `"one"` will estimate one shared correlation
#'                           parameter across species, and `"all"` will estimate correlation parameters
#'                           across all combinations of species.
#' @param temporal_cor       Correlation structure across time (`phi`). `"none"` assumes no temporal dependence
#'                           in the process errors, `"rw"` assumes a random walk, and `"ar1"` fits an
#'                           AR1 process, estimating an extra parameter.
#' @param K_groups           Formula pointing to the grouping variable to use to estimate `K` and
#'                           aggregate biomass in the production equation. Biomass from all species
#'                           (stocks) will be aggregated and one `K` value estimated if set to `NULL`.
#' @param pe_covariates      Formula describing relationship between surplus production (process error)
#'                           and covariates. No covariates are applied if set to `NULL`.
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
                      species_cor = "none",
                      temporal_cor = "none",
                      K_groups = NULL,
                      pe_covariates = NULL,
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
    index$survey <- factor(index$survey)

    ## Set-up model matrix | formula with covariates
    if (is.null(pe_covariates)) {
        pe_model_mat <- matrix(rep(0, nrow(landings)), ncol = 1)
    } else {
        pe_model_mat <- model.matrix(pe_covariates, data = landings)[, -1, drop = FALSE] # drop intercept because that is defined by the process errors
    }
    if (is.null(K_groups)) {
        K_map <- rep(0, nlevels(landings$species))
        B_group_mat <- matrix(1, nrow = nlevels(landings$species), ncol = nlevels(landings$species))
    } else {

        if (length(all.vars(K_groups)) > 1) {
            stop("Please supply only one variable to the K_groups argument; mapping by more than one variable has yet to be implemented.")
        }

        sp_group <- unique(landings[, c("species", all.vars(K_groups))])
        sp_group <- sp_group[order(sp_group$species), ]
        if (nrow(sp_group) > nlevels(landings$species)) {
            stop("All species are not nested within the K_groups variable. Please check groupings.")
        }

        ## Set up a matrix of 0s and 1s to control the summation of biomass
        f <- as.formula(paste(Reduce(paste, deparse(K_groups)), "-1"))
        mm <- model.matrix(f, data = sp_group)
        ind <- rep(seq(ncol(mm)), colSums(mm))
        B_group_mat <- unname(mm[, ind])

        ## And set up a map for the K parameters
        K_map <- as.numeric(factor(sp_group[[2]])) - 1

    }

    ## Center index and landings by the mean(log(index) ~ K_group) to aid convergence
    by_sp_yr_grp <- paste(c("species", "year", all.vars(K_groups)), collapse = " + ")
    by_yr_grp <- paste(c("year", all.vars(K_groups)), collapse = " + ")
    by_grp <- ifelse(is.null(K_groups), "0", paste(all.vars(K_groups), collapse = " + "))
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

    ## Values to keep
    keep <- rep(1L, length(index$index))
    if (!is.null(leave_out)) {
        keep[leave_out] <- 0L
    }

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
                I_survey = as.numeric(index$survey) - 1,
                I_sy = as.numeric(index$sy) - 1,
                min_B = 0.0001,
                nY = nlevels(landings$y),
                nS = nlevels(landings$species),
                log_K_option = as.integer(log_K_option$option) - 1,
                log_B0_option = as.integer(log_B0_option$option) - 1,
                log_r_option = as.integer(log_r_option$option) - 1,
                log_sd_B_option = as.integer(log_sd_B_option$option) - 1,
                log_q_option = as.integer(log_q_option$option) - 1,
                log_sd_I_option = as.integer(log_sd_I_option$option) - 1,
                logit_rho_option = as.integer(logit_rho_option$option) - 1,
                logit_phi_option = as.integer(logit_phi_option$option) - 1,
                lower_log_K = log_K_option$lower - log_center, # lots of repetition - todo: find better solution
                upper_log_K = log_K_option$upper - log_center,
                lower_log_B0 = log_B0_option$lower - log_center,
                upper_log_B0 = log_B0_option$upper - log_center,
                lower_log_r = log_r_option$lower,
                upper_log_r = log_r_option$upper,
                lower_log_sd_B = log_sd_B_option$lower,
                upper_log_sd_B = log_sd_B_option$upper,
                lower_log_q = log_q_option$lower,
                upper_log_q = log_q_option$upper,
                lower_log_sd_I = log_sd_I_option$lower,
                upper_log_sd_I = log_sd_I_option$upper,
                lower_logit_rho = logit_rho_option$lower,
                upper_logit_rho = logit_rho_option$upper,
                mean_logit_phi = logit_phi_option$mean,
                sd_logit_phi = logit_phi_option$sd,
                lower_logit_phi = logit_phi_option$lower,
                upper_logit_phi = logit_phi_option$upper,
                dmuniform_sd = 0.1, # controls how sharp the approximate uniform distribution is
                pe_covariates = pe_model_mat,
                K_map = K_map,
                B_groups = B_group_mat,
                keep = keep)

    par <- list(log_B = matrix(ceiling(mean_log_index$index),
                               ncol = dat$nS, nrow = dat$nY, byrow = TRUE),
                mean_log_sd_B = log_sd_B_option$mean,
                log_sd_log_sd_B = log(log_sd_B_option$sd),
                log_sd_B = rep(-1, nlevels(landings$species)),
                mean_logit_rho = logit_rho_option$mean,
                log_sd_logit_rho = log(logit_rho_option$sd),
                logit_rho = rep(0, n_rho),
                logit_phi = 0,
                log_K = ceiling(log(max_landings)),
                mean_log_K = log_K_option$mean - log_center,
                log_sd_log_K = log(log_K_option$sd),
                mean_log_B0 = log_B0_option$mean - log_center,
                log_sd_log_B0 = log(log_B0_option$sd),
                log_B0 = ceiling(mean_log_index$index),
                mean_log_r = log_r_option$mean,
                log_sd_log_r = log(log_r_option$sd),
                log_r = rep(-1, nlevels(landings$species)),
                log_m = rep(log(2), nlevels(landings$species)),
                mean_log_q = log_q_option$mean,
                log_sd_log_q = log(log_q_option$sd),
                log_q = rep(-1, nlevels(index$survey)),
                mean_log_sd_I = log_sd_I_option$mean,
                log_sd_log_sd_I = log(log_sd_I_option$sd),
                log_sd_I = rep(-1, nlevels(index$survey)),
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
        par$logit_phi <- -15 # results in a value very close to 0
    }
    if (temporal_cor == "rw") {
        map$logit_phi <- factor(NA)
        par$logit_phi <- 15 # results in a value very close to 1
    }

    if (logit_phi_option$option %in% c("random", "coupled")) {
        stop("The 'random' or 'coupled' options are not applicable for parametere 'logit_phi'.")
    }

    if (log_K_option$option != "random") {
        map$mean_log_K <- map$log_sd_log_K <- factor(NA)
        if (log_K_option$option == "coupled") {
            stop("The 'coupled' option is not applicable for parametere 'log_K'. Please control parameter grouping using the K_groups argument.")
        }
    }

    if (log_sd_B_option$option != "random") {
        map$mean_log_sd_B <- map$log_sd_log_sd_B <- factor(NA)
        if (log_sd_B_option$option == "coupled") {
            map$log_sd_B <- factor(rep(1, length(par$log_sd_B)))
        }
    }


    if (log_B0_option$option != "random") {
        map$mean_log_B0 <- map$log_sd_log_B0 <- factor(NA)
        if (log_B0_option$option == "coupled") {
            map$log_B0 <- factor(rep(1, length(par$log_B0)))
        }
    }

    if (log_r_option$option != "random") {
        map$mean_log_r <- map$log_sd_log_r <- factor(NA)
        if (log_r_option$option == "coupled") {
            map$log_r <- factor(rep(1, length(par$log_r)))
        }
    }

    if (log_q_option$option != "random") {
        map$mean_log_q <- map$log_sd_log_q <- factor(NA)
        if (log_q_option$option == "coupled") {
            map$log_q <- factor(rep(1, length(par$log_q)))
        }
    }

    if (log_sd_I_option$option != "random") {
        map$mean_log_sd_I <- map$log_sd_log_sd_I <- factor(NA)
        if (log_sd_I_option$option == "coupled") {
            map$log_sd_I <- factor(rep(1, length(par$log_sd_I)))
        }
    }

    if (logit_rho_option$option != "random") {
        map$mean_logit_rho <- map$log_sd_logit_rho <- factor(NA)
        if (logit_rho_option$option == "coupled") {
            map$logit_rho <- factor(rep(1, length(par$logit_rho)))
        }
    }

    if (is.null(pe_covariates)) {
        map$pe_betas <- factor(NA)
    }

    random <- "log_B"
    if (log_K_option$option == "random") {
        random <- c(random, "log_K")
    }
    if (log_sd_B_option$option == "random") {
        random <- c(random, "log_sd_B")
    }
    if (log_B0_option$option == "random") {
        random <- c(random, "log_B0")
    }
    if (log_r_option$option == "random") {
        random <- c(random, "log_r")
    }
    if (log_q_option$option == "random") {
        random <- c(random, "log_q")
    }
    if (log_sd_I_option$option == "random") {
        random <- c(random, "log_sd_I")
    }
    if (logit_rho_option$option == "random") {
        random <- c(random, "logit_rho")
    }

    ## Fit model
    if (!is.null(start_par)) par <- start_par
    obj <- MakeADFun(dat, par, map = map, random = random, DLL = "multispic",
                     silent = silent)
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
    index$left_out <- !as.logical(keep)

    pop <- landings
    pop$pe <- rep$log_B_std_res
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

        if (is.null(K_groups)) {
            tot_pop$K <- exp(par$log_K)
            tot_pop$K_lwr <- exp(par_lwr$log_K)
            tot_pop$K_upr <- exp(par_upr$log_K)
        } else {
            ind <- as.numeric(factor(tot_pop[, all.vars(K_groups)]))
            tot_pop$K <- exp(par$log_K[ind])
            tot_pop$K_lwr <- exp(par_lwr$log_K[ind])
            tot_pop$K_upr <- exp(par_upr$log_K[ind])
        }

    }

    ## Calculate marginal AIC
    mAIC <- 2 * length(opt$par) + 2 * opt$objective

    end_time <- Sys.time()
    run_dur <- end_time - start_time

    out <- list(call = call, run_dur = run_dur, log_center = log_center, obj = obj, opt = opt,
                sd_rep = sd_rep, rep = rep, par = par, se = se, par_lwr = par_lwr,
                par_upr = par_upr, index = index, landings = landings,
                pop = pop, tot_pop = tot_pop, mAIC = mAIC)

}


#' Function for running leave one out cross-validation
#'
#' @param fit         Object from [fit_model()]
#' @param progress    Display progress bar? (Generated using [progress::progress_bar])
#'
#' @return Returns a list with three objects:
#'    1) `obs`  -  log observations that were left out at each step
#'    2) `pred` -  log predictions at each step, and
#'    3) `mse`  -  mean squared error of the predictions (leave one out cross validation score).
#'
#' @export
#'

run_loo <- function(fit, progress = TRUE) {

    n <- length(fit$index$index)
    obs <- pred <- numeric(n)

    if (!is.null(fit$sd_rep)) {
        start_par <- as.list(fit$sd_rep, "Est")
    } else {
        stop("Object sd_rep is NA in the supplied fit object. Please re-run model with light = FALSE.")
    }

    if (progress) {
        if (!requireNamespace("progress", quietly = TRUE)) {
            stop("The progress package is needed to display a progress bar. Please install it.", call. = FALSE)
        }
        pb <- progress::progress_bar$new(
            format = "[:bar] :percent (:current / :total) in :elapsed (eta: :eta)",
            total = n, clear = FALSE, show_after = 0, width = 100)
    }

    for (i in seq(n)) {
        f <- try(update(fit, leave_out = i, start_par = start_par,
                        light = TRUE, silent = TRUE))
        if (class(f) == "try-catch" | fit$opt$message == "false convergence (8)") {
            obs[i] <- pred[i] <- NA
        } else {
            obs[i] <- f$index$log_index[f$index$left_out]
            pred[i] <- f$index$log_pred_index[f$index$left_out]
        }
        if (progress) pb$tick()
    }

    if (any(is.na(pred))) {
        warning(paste("When iterating across ", n, " observations, model fitting failed for ", sum(is.na(pred)), " cases when an observation was left out."))
    }

    list(obs = obs, pred = pred, mse = mean((obs - pred) ^ 2, na.rm = TRUE))

}


#' Function for running a retrospective analysis
#'
#' @details This function runs a retrospective analysis whereby terminal survey index estimates
#'          are excluded from the analysis. A hindcast is also preformed in each fold where
#'          observed values are left out but predicted to measure the models forecasting skill.
#'
#' @param fit         Object from [fit_model()]
#' @param folds       Number of years to 'fold' back
#' @param progress    Display progress bar? (Generated using [progress::progress_bar])
#'
#' @return Returns a list with three objects:
#'    1) `retro_fits`  -  a list including a series of fits from each retrospective fold
#'    2) `hindcasts` -  a data.frame with observed, but left out, and predicted indices from each fold
#'    3) `mse`  -  mean squared error of the hindcasts (measure of forecasting skill)
#'
#' @export
#'

run_retro <- function(fit, folds = 10, progress = TRUE) {

    if (!is.null(fit$sd_rep)) {
        start_par <- as.list(fit$sd_rep, "Est")
    } else {
        stop("Object sd_rep is NA in the supplied fit object. Please re-run model with light = FALSE.")
    }

    if (progress) {
        if (!requireNamespace("progress", quietly = TRUE)) {
            stop("The progress package is needed to display a progress bar. Please install it.", call. = FALSE)
        }
        pb <- progress::progress_bar$new(
            format = "[:bar] :percent (:current / :total) in :elapsed (eta: :eta)",
            total = folds, clear = FALSE, show_after = 0, width = 100)
    }

    terminal_year <- max(fit$index$year)
    retro_years <- terminal_year - seq(folds)
    retro_fits <- hindcasts <- vector("list", folds)
    names(retro_fits) <- names(hindcasts) <- as.character(retro_years)

    for (i in seq(folds)) {

        index <- fit$index
        landings <- fit$landings

        ## Subset the input data on year at a time, keeping one extra year of indices and landings
        retro_index <- index[index$year <= retro_years[i] + 1, ]
        retro_landings <- landings[landings$year <= retro_years[i] + 1, ]
        retro_inputs <- list(landings = retro_landings, index = retro_index)

        ## Identify observations to leave out, but provide predictions for these values
        ind <- retro_index$year == retro_years[i] + 1

        ## TODO:
        ## - allow start_par to be used  (requires careful change to par dimensions)
        fit <- try(update(fit, inputs = retro_inputs, leave_out = ind,
                          light = TRUE, silent = TRUE))

        if (class(fit) == "try-catch" | fit$opt$message == "false convergence (8)") {

            retro_fits[[i]] <- NA
            hindcasts[[i]] <- NULL

        } else {
            retro_fits[[i]] <- fit

            if (sum(ind) > 0) {
                hindcasts[[i]] <- fit$index[fit$index$left_out,
                                            c("year", "survey", "species", "log_index", "log_pred_index")]
                hindcasts[[i]]$retro_year <- retro_years[i]
            } else {
                hindcasts[[i]] <- NULL
                warning(paste("A hindcast is moot for", retro_years[i], "as there is no index to predict."))
            }

        }

        if (progress) pb$tick()

    }

    if (any(is.na(retro_fits))) {
        warning(paste("While folding back ", folds, " years, model fitting failed in ", sum(is.na(retro_fits)), " cases."))
    }

    hindcasts <- do.call(rbind, hindcasts)

    mse <- mean((hindcasts$log_index - hindcasts$log_pred_index) ^ 2)
    list(retro_fits = retro_fits, hindcasts = hindcasts, mse = mse)

}


