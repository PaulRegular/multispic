

#' Logit and inverse logit functions
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
#' @export
inv_logit <- function(x, shift = FALSE) {
    if (!shift) { 1 / (1 + exp(-x)) } else { 2 / (1 + exp(-x)) - 1 }
}

