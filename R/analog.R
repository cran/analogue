###########################################################################
##                                                                       ##
## analog - function to perform analogue matching between samples        ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## x, y          - two data frames with same columns. 'x' is training    ##
##                 data and 'y', the test data                           ##
## method        - character string naming the dissimilarity methods to  ##
##                 be used.                                              ##
## keep.train    - logical; keep the training set dissimilarity matrix?  ##
##                                                                       ##
###########################################################################
analog <- function(x, ...) UseMethod("analog")

analog.default <- function(x, y, method = c("euclidean", "SQeuclidean",
                                   "chord", "SQchord", "bray",
                                   "chi.square", "SQchi.square",
                                   "information", "chi.distance",
                                   "manhattan", "kendall", "gower",
                                   "alt.gower", "mixed"),
                           keep.train = TRUE, ...)
  {
    if(!is.matrix(x))
      x <- as.matrix(x)
    if(!is.matrix(y))
      y <- as.matrix(y)
    dissim <- distance(x = x, y = y, method = method)
    retval <- list(analogs = dissim)
    if(keep.train)
        retval$train <- distance(x = x, method = method)
    retval$call <- match.call()
    retval$method <- method
    class(retval) <- "analog"
    return(retval)
  }

print.analog <- function(x, probs = c(0.01, 0.02, 0.05, 0.1, 0.2),
                         digits = min(3, getOption("digits") - 4), ...)
  {
    method <- x$method
    .call <- deparse(x$call)
    cat("\n")
    writeLines(strwrap("Analogue matching for fossil samples", prefix = "\t"))
    cat(paste("\nCall:", .call, "\n"))
    cat(paste("Dissimilarity:", method, "\n"))
    minDij <- apply(x$analogs, 2, min)
    which.minDij <- apply(x$analogs, 2, which.min)
    if(!is.null(x$train))
       {
         cat("\nQuantiles of the dissimilarities for the training set:\n\n")
         print(quantile(as.vector(as.dist(x$train)), probs = probs),
               digits = digits)
       }
    cat("\nClosest modern analogue from training set:\n\n")
    print(data.frame(ID = rownames(x$analogs)[which.minDij],
                     minDC = format(minDij, digits = digits)),
          print.gap = 2)
    cat("\n")
    invisible(x)
  }
