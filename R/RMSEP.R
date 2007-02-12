RMSEP <- function(object, ...) UseMethod("RMSEP")

RMSEP.default <- function(object, ...)
  {
    stop("No default method for \"RMSEP\"")
  }

RMSEP.bootstrap <- function(object, type = c("birks1990", "standard"), ...) {
  if(!inherits(object, "bootstrap"))
    stop("'object' is not of class \"bootstrap\".")
  if(missing(type))
    type <- "birks1990"
  if(type == "birks1990")
    rmsep <- object$bootstrap$rmsep[k(object)]
  else
    rmsep <- sqrt(mean(object$bootstrap$residuals[, k(swap.boot)]))
  return(rmsep)
}
