


#' Fit a multispecies surplus production model
#'
#' @param inputs        List that includes landings and index data
#' @param q_groups      Name of column in the index data to group the q parameter estimates by
#' @param q_option      Vector of three: a mean value of log q, an sd value of log q, and a value
#'                      of 0, 1 or 2. The last value toggles whether the first two values are
#'                      ignored (0), used as starting values (1), or priors (2).
#' @param cor_str       Correlation structure across species. "none" will not estimate
#'                      correlations across species, "one" will estimate one shared correlation
#'                      parameter across species, and "all" will estimate correlation parameters
#'                      across all combinations of species.
#'
#' @return
#' @export
#'

fit_model <- function(inputs,
                      q_groups = "gear_season",
                      q_option = c(mu_log_q = 0, sd_log_q = 1, option = 0),
                      cor_str = "one") {

    landings <- inputs$landings
    index <- inputs$index

    ## Scale index and landings to aid convergence
    scaler <- sd(index$index)
    index$index <- index$index / scaler
    landings$landings <- landings$landings / scaler

    n_cor <- sum(lower.tri(matrix(NA, nrow = nlevels(landings$species),
                                  ncol = nlevels(landings$species))))
    dat <- list(L = as.numeric(landings$landings),
                L_species = as.numeric(landings$species) - 1,
                L_year = as.numeric(landings$y) - 1,
                I = as.numeric(index$index),
                I_species = as.numeric(index$species) - 1,
                I_survey = as.numeric(index[, q_groups]) - 1,
                I_sy = as.numeric(index$sy) - 1,
                min_B = 0.001,
                nY = max(as.numeric(landings$y)),
                nS = max(as.numeric(landings$species)),
                q_option = as.integer(q_option[3]))
    par <- list(log_B = matrix(0, nrow = dat$nY, ncol = dat$nS),
                log_sd_B = rep(-1, nlevels(landings$species)),
                logit_cor = rep(0, n_cor),
                log_K = 2,
                log_r = rep(-1, nlevels(landings$species)),
                log_m = rep(log(2), nlevels(landings$species)),
                mu_log_q = q_option[1],
                log_sd_log_q = log(q_option[2]),
                log_q = rep(-1, nlevels(index[, q_groups])),
                log_sd_I = rep(-1, nlevels(index[, q_groups])))

    map <- list(log_m = factor(rep(NA, nlevels(landings$species))))
    if (cor_str == "one") {
        map$logit_cor <- factor(rep(1, length(par$logit_cor)))
    }
    if (cor_str == "none") {
        map$logit_cor <- factor(rep(NA, length(par$logit_cor)))
    }
    if (q_option[3] %in% c(0, 2)) {
        map$mu_log_q <- map$log_sd_log_q <- factor(NA)
    }

    random <- "log_B"
    if (q_option[3] == 1) {
        random <- c(random, "log_q")
    }

    ## Fit model
    obj <- MakeADFun(dat, par, map = map, random = random, DLL = "multispic")
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 1000, iter.max = 1000))
    sd_rep <- sdreport(obj)
    rep <- obj$report()

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

    ## Extract process error
    pe <- data.frame(year = landings$year,
                     species = landings$species,
                     stock = landings$stock,
                     pe = rep$log_B_std_res)

    ## Extract biomass
    biomass <- data.frame(year = landings$year,
                          species = landings$species,
                          stock = landings$stock,
                          B = exp(est$log_B_vec) * scaler,
                          B_lwr = exp(lwr$log_B_vec) * scaler,
                          B_upr = exp(upr$log_B_vec) * scaler)

    list(scaler = scaler, obj = obj, opt = opt, sd_rep = sd_rep, rep = rep,
         index = index, landings = landings, pe = pe, biomass = biomass)

}



