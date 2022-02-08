#' Make a flexdashboard for visualizing the model fits
#'
#' @param fit            Object produced by [multispic()].
#' @param ...            Additional arguments to send to [rmarkdown::run()]
#'
#' @export
#'

vis_multispic <- function(fit, ...) {

  pkg <- c("rmarkdown", "shiny", "flexdashboard", "plotly", "viridis")
  for (p in pkg) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(paste(p, "is needed for vis_multispic to work. Please install it."), call. = FALSE)
    }
  }

  rmd_env <- new.env()
  rmd_env$fit <- fit
  rmarkdown::run(file = system.file("rmd", "vis_multispic.Rmd", package = "multispic"),
                 render_args = list(envir = rmd_env), ...)

}

