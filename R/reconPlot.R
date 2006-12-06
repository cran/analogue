###########################################################################
##                                                                       ##
## reconPlot - function to draw reconstructed environmental variables    ##
##             from transfer function models                             ##
##                                                                       ##
## Created       : 05-Nov-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.0-1                                                 ##
## Last modified : 05-Nov-2006                                           ##
##                                                                       ##
###########################################################################
##
## Arguments
##
## S3 method
reconPlot <- function(x, ...) UseMethod("reconPlot")

reconPlot.default <- function(x, ...)
  {
    stop("No default method for \"reconPlot\"")
  }

reconPlot.predict.mat <- function(x, depths, use.labels = FALSE,
                                  predictions = c("apparent",
                                    "bootstrap"),
                                  error.bars = FALSE,
                                  sample.specific = TRUE,
                                  rev.x = TRUE,
                                  type = "l",
                                  xlim, ylim,
                                  xlab = "", ylab = "", main = "",
                                  ...)
  {
    if(missing(predictions))
      predictions <- "apparent"
    predictions <- match.arg(predictions)
    if(missing(depths))
      {
        if(use.labels) {
          if(predictions == "apparent") {
            depths <- as.numeric(colnames(x$predictions$apparent$predicted))
            n.analogues <- x$predictions$apparent$k
            preds <- x$predictions$apparent$predicted[n.analogues, ]
            errors <- x$apparent$rmse[n.analogues]
          } else {
            depths <- as.numeric(rownames(x$predictions$bootstrap$predicted))
            n.analogues <- x$predictions$bootstrap$k
            preds <- x$predictions$bootstrap$predicted[,n.analogues]
            if(sample.specific)
              errors <- x$predictions$sample.errors$rmsep[, n.analogues]
            else
              errors <- x$bootstrap$rmsep[n.analogues]
          }
        } else {
          stop("If \"use.labels = FALSE\", then \"depths\" must be provided.")
        }
      }
    if(missing(xlim))
      xlim <- range(depths)
    if(rev.x)
      xlim <- rev(xlim)
    upper <- preds + errors
    lower <- preds - errors
    if(missing(ylim)) {
      if(error.bars)
        ylim <- range(preds, upper, lower)
      else
        ylim <- range(preds)
    }
    plot(depths, preds, ylim = ylim, xlim = xlim, type = "n",
         ylab = ylab, xlab = xlab, main = main, ...)
    if(error.bars)
      arrows(depths, upper, depths, lower, length = 0.02, angle = 90,
             code = 3, col = "grey")
    lines(depths, preds, type = type, ...)
    invisible()
  }
    

#reconPlot.mat <- function(x, fossil, k, depths, use.labels = FALSE, ...)
#  {
#    if(missing(k))
#      env <- predict(x, newdata = fossil)
#    if(missing(depths)) {
#      if(use.labels)
#        depths <- x$predictions$apparent$predicted
#    }
#}
#reconPlot(rlgh.mat, use.labels = TRUE, xlab = "Depth", ylab = "pH",
#          error.bars = TRUE, predictions = "bootstrap")

