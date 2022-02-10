#' Make a flexdashboard for visualizing the model fits
#'
#' @param fit            Object produced by [multispic()].
#' @param output_file    Name of file to export using [rmarkdown::render()].
#'                       If `NULL`, flexdashboard will be rendered using [rmarkdown::run()]
#' @param ...            Additional arguments to send to [rmarkdown::run()] or
#'                       [rmarkdown::render()].
#'
#' @export
#'

vis_multispic <- function(fit, output_file = NULL, ...) {

  pkg <- c("rmarkdown", "shiny", "flexdashboard", "plotly", "viridis")
  for (p in pkg) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(paste(p, "is needed for vis_multispic to work. Please install it."), call. = FALSE)
    }
  }

  rmd_file <- system.file("rmd", "vis_multispic.Rmd", package = "multispic")

  rmd_env <- new.env()
  rmd_env$fit <- fit

  if (is.null(output_file)) {
    rmarkdown::run(file = rmd_file,
                   render_args = list(envir = rmd_env), ...)
  } else {
    output_dir <- normalizePath(dirname(output_file))
    output_file <- basename(output_file)
    rmarkdown::render(input = rmd_file,
                      output_file = output_file,
                      output_dir = output_dir, ...)
  }

}

