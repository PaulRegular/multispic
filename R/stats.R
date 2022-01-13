

#' Logit and invers logit functions
#'
#' @param x     Probability or logit to transform
#' @param shift Shift transformation such that the domain is between -1 and 1? Defaults to domain of 0 to 1.
#'
#' @export
#'
#' @rdname logit
#'
#' @examples
#'
#' logit(0.9)
#' inv_logit(2.197225)
#'

logit <- function(x, shift = FALSE) {
    if (!shift) { log(x / (1 - x)) } else { log((1 + x) / (1 - x))}
}

#' @rdname logit
inv_logit <- function(x, shift = FALSE) {
    if (!shift) { 1 / (1 + exp(-x)) } else { 2 / (1 + exp(-x)) - 1 }
}


#' Mollified uniform density using a Gaussian distribution
#'
#' @description Density distribution function for a mollified uniform using the
#'              Gaussian distribution.
#'
#' @param x     Vector of quantiles.
#' @param lower Lower range for the distribution.
#' @param upper Upper range for the distribution.
#' @param sd    Standard deviation to control the shape of the distribution. A value of 0 will
#'              be equivalent to a uniform distribution and values approach the normal as sd
#'              approach infinity.
#'
#' @export
#'
#' @example
#' x <- seq(-1, 1, length = 100)
#' d <- dmuniform(x, lower = -0.5, upper = 0.5, sd = 0.05)
#' plot(x, d, type = "l", ylab = "Density")
#'

dmuniform <- function(x, lower = 0, upper = 1, sd = 0.1, log = FALSE) {

    ## R translation of C++ function in scr/utils.hpp file
    mean <- (lower + upper) / 2
    z <- (x - mean) / sd
    delta_u <- (upper - mean) / sd
    delta_l <- -((lower - mean) / sd)
    d <- ifelse(z > 0,
                exp(log(pnorm(-(z - delta_l)))) - exp(log(pnorm(-(z + delta_u)))),
                exp(log(pnorm(z + delta_u))) - exp(log(pnorm(z - delta_l))))
    if (log) d <- log(d)
    d

}


