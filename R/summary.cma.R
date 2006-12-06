###########################################################################
##                                                                       ##
## cma           - extracts and formats close modern analogues           ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object for method dispatch. Only class 'analog'.      ##
##                                                                       ##
###########################################################################
summary.cma <- function(object, ...)
  {
    class(object) <- "summary.cma"
    object
  }

print.summary.cma <- function(x,
                              digits = min(3, getOption("digits") - 4), ...)
  {
    class(x) <- "cma"
    print(x)
    cat("\nDistances:\n\n")
    print(x$distances, digits = digits, na.print = "")
    cat("\nSamples:\n\n")
    print(x$samples, quote = FALSE, right = TRUE, na.print = "")
    invisible(x)
  }
