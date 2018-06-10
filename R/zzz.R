
#' @import TMB
#' @useDynLib MSP

.onUnload <- function(lib) {
    library.dynam.unload("MSP", lib)
}
