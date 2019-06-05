


#' Helper for defining the settings of a parameter
#'
#' @param option    Should the parameter be estimated freely ("fixed") or revolve
#'                  around an estimated mean and sd ("random"), a fixed mean
#'                  ("prior_mean") or a fixed mean and sd ("prior").
#' @param mean      Mean value of the parameter. Treated as a starting value if
#'                  option is "random" or as a prior if option is option is "prior".
#'                  Ignored if option is "fixed".
#' @param sd        SD value of the parameter. Treated as a starting value is
#'                  option is "random" or a prior if option is "prior". Ignored
#'                  if option is "fixed".
#'
#' @export
#'

par_option <- function(option = "fixed", mean = 0, sd = 1) {
    list(option = factor(option, levels = c("fixed", "random", "prior")),
         mean = mean, sd = sd)
}


#' Fit a multispecies surplus production model
#'
#' @param inputs             List that includes landings, index and covariate (optional) data.
#' @param survey_group       Name of column in the index data to group the survey parameter estimates by
#' @param scaler             Number to scale values by to aid convergence.
#' @param log_B0_option      Settings for the estimation of the starting biomass;
#'                           define using \code{\link{par_option}}.
#' @param log_r_option       Settings for the estimation of log_r; define using \code{\link{par_option}}.
#' @param log_sd_B_option    Settings for the estimation of sd for the process; define using
#'                           \code{\link{par_option}}.
#' @param log_q_option       Settings for the estimation of log_q; define using \code{\link{par_option}}.
#' @param log_sd_I_option    Settings for the estimation of sd for the indices; define using
#'                           \code{\link{par_option}}.
#' @param logit_cor_option   Setting for the estimation of the correlation across stocks; define using
#'                           \code{\link{par_option}}.
#' @param cor_str            Correlation structure across species. "none" will not estimate
#'                           correlations across species, "one" will estimate one shared correlation
#'                           parameter across species, and "all" will estimate correlation parameters
#'                           across all combinations of species.
#' @param formula            Formula describing relationship between surplus production and covariates.
#'                           Not used if set to NULL.
#' @param leave_out          Specific index values to leave out from the analysis (row number).
#'                           Useful for cross-validation. All data are kept if NULL.
#'
#' @export
#'

fit_model <- function(inputs,
                      scaler = sd(inputs$index$index),
                      survey_group = "gear_season",
                      log_B0_option = par_option(),
                      log_r_option = par_option(),
                      log_sd_B_option = par_option(),
                      log_q_option = par_option(),
                      log_sd_I_option = par_option(),
                      logit_cor_option = par_option(),
                      cor_str = "one",
                      formula = NULL,
                      leave_out = NULL) {

    call <- match.call()

    landings <- inputs$landings
    index <- inputs$index
    covariates <- inputs$covariates

    ## Set-up model matrix | formula with covariates
    if (is.null(formula)) {
        model_mat <- matrix(rep(0, nrow(landings)), ncol = 1)
    } else {
        f <- as.formula(paste(Reduce(paste, deparse(formula)), "-1")) # drop intercept
        model_mat <- model.matrix(f, data = covariates)
    }

    ## Scale index and landings to aid convergence
    index$index <- index$index / scaler
    landings$landings <- landings$landings / scaler

    ## Values to keep
    keep <- rep(1L, length(index$index))
    if (!is.null(leave_out)) {
        keep[leave_out] <- 0L
    }

    ## Set-up the objects for TMB
    n_cor <- sum(lower.tri(matrix(NA, nrow = nlevels(landings$species),
                                  ncol = nlevels(landings$species))))
    dat <- list(L = as.numeric(landings$landings),
                L_species = as.numeric(landings$species) - 1,
                L_year = as.numeric(landings$y) - 1,
                I = as.numeric(index$index),
                I_species = as.numeric(index$species) - 1,
                I_survey = as.numeric(index[, survey_group]) - 1,
                I_sy = as.numeric(index$sy) - 1,
                min_B = 0.001,
                nY = max(as.numeric(landings$y)),
                nS = max(as.numeric(landings$species)),
                log_B0_option = as.integer(log_B0_option$option) - 1,
                log_r_option = as.integer(log_r_option$option) - 1,
                log_sd_B_option = as.integer(log_sd_B_option$option) - 1,
                log_q_option = as.integer(log_q_option$option) - 1,
                log_sd_I_option = as.integer(log_sd_I_option$option) - 1,
                logit_cor_option = as.integer(logit_cor_option$option) - 1,
                covariates = model_mat,
                keep = keep)
    par <- list(log_B = matrix(0, nrow = dat$nY, ncol = dat$nS),
                mean_log_sd_B = log_sd_B_option$mean,
                log_sd_log_sd_B = log(log_sd_B_option$sd),
                log_sd_B = rep(-1, nlevels(landings$species)),
                mean_logit_cor = logit_cor_option$mean,
                log_sd_logit_cor = log(logit_cor_option$sd),
                logit_cor = rep(0, n_cor),
                logit_phi = 0,
                mean_log_B0 = log_B0_option$mean,
                log_sd_log_B0 = log(log_B0_option$sd),
                log_B0 = rep(0, nlevels(landings$species)),
                log_K = 2,
                mean_log_r = log_r_option$mean,
                log_sd_log_r = log(log_r_option$sd),
                log_r = rep(-2, nlevels(landings$species)),
                log_m = rep(log(2), nlevels(landings$species)),
                mean_log_q = log_q_option$mean,
                log_sd_log_q = log(log_q_option$sd),
                log_q = rep(-0.5, nlevels(index[, survey_group])),
                mean_log_sd_I = log_sd_I_option$mean,
                log_sd_log_sd_I = log(log_sd_I_option$sd),
                log_sd_I = rep(-1, nlevels(index[, survey_group])),
                betas =  rep(0, ncol(model_mat)))

    map <- list(log_m = factor(rep(NA, nlevels(landings$species))))
    if (cor_str == "one") {
        map$logit_cor <- factor(rep(1, length(par$logit_cor)))
    }
    if (cor_str == "none") {
        map$logit_cor <- factor(rep(NA, length(par$logit_cor)))
        dat$logit_cor_option <- 0 # skip prior / random effect loop
    }

    if (log_sd_B_option$option %in% c("fixed", "prior")) {
        map$mean_log_sd_B <- map$log_sd_log_sd_B <- factor(NA)
    }

    if (log_B0_option$option %in% c("fixed", "prior")) {
        map$mean_log_B0 <- map$log_sd_log_B0 <- factor(NA)
    }

    if (log_r_option$option %in% c("fixed", "prior")) {
        map$mean_log_r <- map$log_sd_log_r <- factor(NA)
    }

    if (log_q_option$option %in% c("fixed", "prior")) {
        map$mean_log_q <- map$log_sd_log_q <- factor(NA)
    }

    if (log_sd_I_option$option %in% c("fixed", "prior")) {
        map$mean_log_sd_I <- map$log_sd_log_sd_I <- factor(NA)
    }

    if (logit_cor_option$option %in% c("fixed", "prior")) {
        map$mean_logit_cor <- map$log_sd_logit_cor <- factor(NA)
    }

    if (is.null(formula)) {
        map$betas <- factor(NA)
    }

    random <- "log_B"
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
    if (logit_cor_option$option == "random") {
        random <- c(random, "logit_cor")
    }

    ## Fit model
    obj <- MakeADFun(dat, par, map = map, random = random, DLL = "multispic")
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 1000, iter.max = 1000))
    sd_rep <- sdreport(obj)
    rep <- obj$report()

    ## Extract par and re-scale
    par <- as.list(sd_rep, "Est")
    par$log_B0 <- log(exp(par$log_B0) * scaler)
    par$log_B <- log(exp(par$log_B) * scaler)
    par$log_K <- log(exp(par$log_K) * scaler)
    se <- as.list(sd_rep, "Std. Error")

    ## Reset scale
    landings$landings <- landings$landings * scaler
    index$index <- index$index * scaler

    ## Extract and append fits
    est <- split(unname(sd_rep$value), names(sd_rep$value))
    sd <- split(sd_rep$sd, names(sd_rep$value))
    lwr <- split(unname(sd_rep$value) - 1.96 * sd_rep$sd, names(sd_rep$val))
    upr <- split(unname(sd_rep$value) + 1.96 * sd_rep$sd, names(sd_rep$val))
    index$pred <- exp(est$log_pred_I) * scaler
    index$pred_lwr <- exp(lwr$log_pred_I) * scaler
    index$pred_upr <- exp(upr$log_pred_I) * scaler
    index$std_res <- rep$log_I_std_res
    index$left_out <- !as.logical(keep)

    ## Extract population estimates
    pop <- data.frame(year = landings$year,
                      species = landings$species,
                      stock = landings$stock,
                      pe = rep$log_B_std_res,
                      B = exp(est$log_B_vec) * scaler,
                      B_lwr = exp(lwr$log_B_vec) * scaler,
                      B_upr = exp(upr$log_B_vec) * scaler,
                      F = exp(est$log_F),
                      F_lwr = exp(lwr$log_F),
                      F_upr = exp(upr$log_F))

    ## Calculate marginal AIC
    mAIC <- 2 * length(opt$par) + 2 * opt$objective

    list(call = call, scaler = scaler, obj = obj, opt = opt, sd_rep = sd_rep,
         rep = rep, par = par, se = se, index = index, landings = landings,
         pop = pop, mAIC = mAIC)

}



