

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
#' @param post_mean    Mean of the posterior.
#' @param post_sd      SD of the posterior.
#' @param names        Names of the posteriors (and priors) when multiple values are provided.
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
                      post_mean = 0, post_sd = 1, names = NULL,
                      length_out = 1000, xlab = "x", ylab = "Density", title = NULL,
                      min_den = 0.00001, trans_fun = NULL) {

    if (!is.null(names)) {
        if (length(names) != length(prior_mean) |
            length(names) != length(post_mean)) {
            stop("The number of names provided do not equal the number of mean values supplied.")
        }
    } else {
        names <- "noname"
    }

    names(prior_mean) <- names(prior_sd) <- names
    prior_lab <- rep(names, each = length_out)

    prior_den_min <- min(prior_mean) - (10 * max(prior_sd))
    prior_den_max <- max(prior_mean) + (10 * max(prior_sd))

    prior_x <- seq(prior_den_min, prior_den_max, length.out = length_out)
    prior_y <- lapply(names, function(nm) {
        dnorm(prior_x, mean = prior_mean[nm], sd = prior_sd[nm])
    })

    prior <- data.frame(x = rep(prior_x, length(prior_mean)),
                        y = unlist(prior_y),
                        name = as.character(prior_lab))

    names(post_mean) <- names(post_sd) <- names
    post_lab <- rep(names, each = length_out)

    post_den_min <- min(post_mean) - (10 * max(post_sd))
    post_den_max <- max(post_mean) + (10 * max(post_sd))

    post_x <- seq(post_den_min, post_den_max, length.out = length_out)
    post_y <- lapply(names, function(nm) {
        dnorm(post_x, mean = post_mean[nm], sd = post_sd[nm])
    })

    post <- data.frame(x = rep(post_x, length(post_mean)),
                       y = unlist(post_y),
                       name = as.character(post_lab))

    prior <- prior[prior$y > min_den, ]
    post <- post[post$y > min_den, ]
    if (!is.null(trans_fun)) {
        prior$x <- trans_fun(prior$x)
        post$x <- trans_fun(post$x)
    }

    x_min <- min(post_mean) - (4 * max(post_sd))
    x_max <- max(post_mean) + (4 * max(post_sd))

    cols <- viridis::viridis(2)
    p <- plot_ly(x = ~x, y = ~y, frame = ~name) %>%
        add_lines(data = prior, line = list(width = 3), visible = show_prior,
                  color = I(cols[2]), name = "Prior") %>%
        add_lines(data = post, line = list(width = 1),
                  color = I(cols[1]), name = "Posterior", showlegend = show_prior) %>%
        layout(xaxis = list(title = xlab, range = list(x_min, x_max)),
               yaxis = list(title = ylab),
               title = title)
    if (length(post_mean) > 1) {
        p <- p %>% animation_button(visible = FALSE)
    }
    p


}



