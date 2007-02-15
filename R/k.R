## functions for extracting and setting the number
## of analogues to be used
k <- function(object, ...) UseMethod("k")

k.default <- function(object, ...) {
  stop("No default method for 'k'")
}

k.mat <- function(object, weighted=FALSE, ...){
  ## check that this is a mat object
  if(class(object) != "mat")
    stop("'object' must be of class 'mat'")
  if(weighted){
    retval <- object$weighted$k
    attr(retval, "auto") <- object$weighted$auto
    attr(retval, "weighted") <- TRUE
  } else {
    retval <- object$standard$k
    attr(retval, "auto") <- object$standard$auto
    attr(retval, "weighted") <- FALSE
  }
  return(retval)
}

k.bootstrap <- function(object, ...) {
  if (!inherits(object, "bootstrap")) 
    stop("Use only with \"bootstrap\" objects")
  retval <- object$bootstrap$k
  attr(retval, "auto") <- object$auto
  attr(retval, "weighted") <- object$weighted
  return(retval)
}

"k<-" <- function(object, weighted=FALSE, value) UseMethod("k<-")

"k<-.default" <- function(object, weighted=FALSE, value) {
  stop("no default replacement method for 'k'")
}

"k<-.mat" <- function(object, weighted=FALSE, value) {
  ## check that this is a mat object
  if(class(object) != "mat")
    stop("'object' must be of class 'mat'")
  ## check that value is not NULL
  if(is.null(value))
    stop("attempt to set NULL number of analogues, 'k'")
  ## check that value is numeric
  ## need to correct this, is.numeric is not integer
  if(!is.numeric(value))
    stop("attempt to set non-integer number of analogues, 'k'")
  if(weighted) {
    object$weighted$k <- value
    object$weighted$auto <- FALSE
  } else {
    object$standard$k <- value
    object$standard$auto <- FALSE
  }
  object
}
