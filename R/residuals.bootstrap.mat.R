###########################################################################
##                                                                       ##
## residuals.bootstrap.mat() - 'residuals' method for MAT models         ##
##                                                                       ##
## Created       : 13-Jun-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 13-Jun-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object on which method dispatch applied (Only 'mat')  ##
## k             - number of analogues to use. If missing 'k' is chosen  ##
##                 automatically as the 'k' that achieves lowest RMSE.   ##
## weighted      - Logical. Should the analysis use weighted mean of env ##
##                 data of analogues as fitted/estimated values?         ##
##                                                                       ##
###########################################################################
residuals.bootstrap.mat <- function(object, which = c("apparent", "bootstrap"),
                                    ...)
  {
    which <- match.arg(which, several.ok = TRUE)
    res <- vector("list")
    if("apparent" %in% which)
      res$apparent <- object$apparent$residuals
    if("bootstrap" %in% which)
      res$bootstrap <- object$bootstrap$residuals
    res$k <- object$k
    if(!is.null(object$bootstrap))
      res$n.boot <- object$n.boot
    res$auto <- object$auto
    res$weighted <- object$weighted
    class(res) <- "residuals.bootstrap.mat"
    return(res)
  }

print.residuals.bootstrap.mat <- function(x,
                                          digits = min(3, getOption("digits") - 3),
                                          ...)
  {
    cat("\n")
    writeLines(strwrap("Modern Analogue Technique bootstrap residuals",
                       prefix = "\t"))
    cat(paste("\nResiduals based on a",
              ifelse(x$weighted, " weighted", ""),
              " model with ", x$k,
              "-closest analogues\n", sep = ""))
    if(x$auto)
      cat("(k chosen from model with lowest RMSE)\n\n")
    else
      cat("\n")
    if(!is.null(x$apparent)) {
      cat("Apparent residuals:\n")
      print(x$apparent, digits = digits)
      cat("\n")
    }
    if(!is.null(x$bootstrap)) {
      cat("Bootstrap residuals:\n")
      print(x$bootstrap, digits = digits)
      cat("\n")
    }
    invisible(x)
  }
