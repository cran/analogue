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
## cutoff        - numeric. Critical value determining level above which ##
##                 samples are defined as close modern analogues         ##
##                                                                       ##
###########################################################################
cma <- function(object, ...) UseMethod("cma")

cma.default <- function(object, ...)
  {
    stop("No default method for \"cma\"")
  }

cma.analog <- function(object, cutoff, prob = c(0.01, 0.025, 0.05), ...)
  {
    if (!inherits(object, "analog")) 
      stop("Use only with \"analog\" objects")
    if(missing(cutoff)) {
      if(is.null(object$train))
        stop("If 'cutoff' is not provided, 'object' must contain\ncomponent \"train\"")
      else
        cutoff <- quantile(dissim(object), probs = 0.025)
    } else {
      if (!is.numeric(cutoff))
        stop("Argument \"cutoff\" must be numeric")
    }
    if(!any(apply(object$analogs, 2, function(x) any(x <= cutoff))))
      stop(paste("No analogues as close or closer than \"cutoff = ",
                 cutoff, "\":\n\tChoose a more suitable value", sep = ""))
    n.samp <- ncol(object$analogs)
    close <- vector("list", length = n.samp)
    nams <- colnames(object$analogs)
    for(i in 1:n.samp)
      {
        x <- sort(object$analogs[,i])
        x <- x[x <= cutoff]
        if(length(x) == 0) {
          x <- NA
          names(x) <- "None"
        }
        close[[i]] <- x
      }
    each.analogs <- sapply(close, length)
    max.analogs <- max(each.analogs)
    samples <- distances <- matrix(NA, nrow = max.analogs,
                               ncol = n.samp)
    for (i in seq(along = close))
      {
        len <- each.analogs[i]
        distances[1:len,i] <- close[[i]]
        samples[1:len,i] <- names(close[[i]])
      }
    rownames(distances) <- rownames(samples) <- 1:max.analogs
    names(each.analogs) <- colnames(distances) <- colnames(samples) <- nams
    each.analogs[each.analogs == 1] <- 0
    structure(list(distances = distances, samples = samples,
                   call = match.call(), cutoff = cutoff,
                   quant = quantile(dissim(object), probs = prob),
                   prob = prob,
                   method = object$method,
                   n.analogs = each.analogs),
              class = "cma")
  }

print.cma <- function(x,
                      digits = min(3, getOption("digits") - 4), ...)
  {
    method <- x$method
    .call <- deparse(x$call)
    cat("\n")
    writeLines(strwrap("Close modern analogues of fossil samples",
                       prefix = "\t"))
    cat(paste("\nCall:", .call, "\n"))
    cat(paste("\nDissimilarity:", method, "\n"))
    cat(paste("Cutoff:", round(x$cutoff, digits), "\n\n"))
    writeLines(strwrap("Number of analogues per fossil sample:",
                       prefix = "\t"))
    cat("\n")
    print(x$n.analogs, digits = digits)
    cat("\n")
    invisible(x)
  }
