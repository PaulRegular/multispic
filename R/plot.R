

#' Add observed and fit values to plotly visualization
#'
#' @export
#'

add_fit <- function(p, x = NULL, yline = NULL, ymarker = NULL,
                    ymin = NULL, ymax = NULL,
                    ..., data = NULL, inherit = TRUE) {
    if (!is.null(ymin)) {
        p <- p %>% add_ribbons(p, x = x, ymin = ymin, ymax = ymax, ...,
                               data = data, inherit = TRUE,
                               line = list(width = 0), opacity = 0.4,
                               showlegend = FALSE)
    }
    if (!is.null(yline)) {
        p <- p %>% add_lines(x = x, y = yline, z = NULL, ...,
                             data = data, inherit = TRUE)
    }
    if (!is.null(ymarker)) {
        p <- p %>% add_markers(p, x = x, y = ymarker, z = NULL, ...,
                               data = data, inherit = TRUE,
                               showlegend = FALSE)
    }
    p
}



#' Plot parameter distributions
#'
#' @description Plot posterior and prior (optional) distributions.
#' Currently limited to the normal distribution.
#'
#' @param prior_mean   Mean of the prior.
#' @param prior_sd     SD of the prior.
#' @param show_prior   Display prior distribution?
#' @param post_mean    Mean of the posterior (can be one or more values).
#' @param post_sd      SD of the posterior (can be one or more values).
#' @param post_names   Names for the posterior.
#' @param length_out   Number of points to plot (more points will result in smoother curves).
#' @param xlab         Label for the x axis.
#' @param ylab         Label for the y axis.
#' @param title        Plot title.
#' @param min_den      Drop points with densities below this threshold.
#' @param trans_fun    Function for transforming values (e.g. `[exp]`, `[inv_logit]`). Ignored if `NULL`.
#'
#' @export
#'
#' @examples
#'
#'

plot_post <- function(prior_mean = 0, prior_sd = 1, show_prior = TRUE,
                            post_mean = 0, post_sd = 1, post_names = "Posterior",
                            length_out = 1000, xlab = "x", ylab = "Density", title = NULL,
                            min_den = 0.00001, trans_fun = NULL) {

    if (length(post_names) != length(post_mean)) {
        stop("length(post_names) != length(post_mean) - Please specify a name for each posterior distribution")
    }

    if (show_prior) {
        den_min <- prior_mean - (10 * prior_sd)
        den_max <- prior_mean + (10 * prior_sd)
    } else {
        den_min <- min(post_mean) - (10 * max(post_sd))
        den_max <- max(post_mean) + (10 * max(post_sd))
    }
    x_min <- min(post_mean) - (4 * max(post_sd))
    x_max <- max(post_mean) + (4 * max(post_sd))
    x <- seq(den_min, den_max, length.out = length_out)
    prior_y <- dnorm(x, mean = prior_mean, sd = prior_sd)

    names(post_mean) <- names(post_sd) <- post_names

    post_lab <- rep(post_names, each = length_out)
    post_y <- lapply(post_names, function(nm) {
        dnorm(x, mean = post_mean[nm], sd = post_sd[nm])
    })

    post <- data.frame(x = rep(x, length(post_mean)),
                       y = unlist(post_y),
                       lab = as.character(post_lab))
    prior <- data.frame(x = x, y = prior_y, lab = "prior")

    prior <- prior[prior$y > min_den, ]
    post <- post[post$y > min_den, ]
    if (!is.null(trans_fun)) {
        prior$x <- trans_fun(prior$x)
        post$x <- trans_fun(post$x)
    }

    plot_ly(x = ~x, y = ~y, color = ~lab) %>%
        add_lines(data = prior, name = "Prior", line = list(width = 3),
                  color = I("black"), colors = viridis::viridis(length(post_mean)),
                  visible = show_prior) %>%
        add_lines(data = post, line = list(width = 1)) %>%
        layout(xaxis = list(title = xlab, range = list(x_min, x_max)),
               yaxis = list(title = ylab),
               title = title)

}



