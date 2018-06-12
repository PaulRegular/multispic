
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

