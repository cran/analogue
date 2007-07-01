###########################################################################
##                                                                       ##
## residuals.bootstrap.mat() - 'residuals' method for MAT models         ##
##                                                                       ##
## Created       : 13-Jun-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.2                                                   ##
## Last modified : 27-Jun-2007                                           ##
##                                                                       ##
###########################################################################
residuals.bootstrap.mat <- function(object, which = c("model", "bootstrap"),
#residuals.bootstrap <- function(object, which = c("model", "bootstrap"),
                                    ...)
  {
    which <- match.arg(which, several.ok = TRUE)
    res <- vector("list")
    if("model" %in% which)
      res$model <- object$model$residuals
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
#print.residuals.bootstrap <- function(x,
                                          digits = min(3, getOption("digits") - 3),
                                          ...)
  {
    cat("\n")
    writeLines(strwrap("Bootstrap residuals", prefix = "\t"))
    cat(paste("Model type:", x$model, "\n"))
    cat(paste("\nResiduals based on a",
              ifelse(x$weighted, " weighted", ""),
              " model with ", x$k,
              "-closest analogues\n", sep = ""))
    if(x$auto)
      cat("(k chosen from model with lowest RMSEP)\n\n")
    else
      cat("\n")
    if(!is.null(x$model)) {
      cat("Model residuals:\n")
      print(x$model, digits = digits)
      cat("\n")
    }
    if(!is.null(x$bootstrap)) {
      cat("Bootstrap residuals:\n")
      print(x$bootstrap, digits = digits)
      cat("\n")
    }
    invisible(x)
  }
