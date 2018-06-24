
#' @import TMB
#' @useDynLib multispic

.onUnload <- function(lib) {
    library.dynam.unload("multispic", lib)
}
