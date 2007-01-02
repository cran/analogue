###########################################################################
##                                                                       ##
## mat - function to perform the modern analogue technique for           ##
##       environmental reconstruction                                    ##
##                                                                       ##
## Created       : 17-Apr-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.0-1                                                 ##
## Last modified : 17-Apr-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
##                                                                       ##
###########################################################################
## x = training data, y = env var of interest
mat <- function(x, ...) UseMethod("mat")

mat.default <- function(x, y,
                        method = c("euclidean", "SQeuclidean", "chord",
                          "SQchord", "bray", "chi.square", "SQchi.square",
                          "information", "chi.distance", "manhattan",
                          "kendall", "gower", "alt.gower", "mixed"),
                        ...)
  {
    dims <- dim(x) # the numbers of samples / species
    site.nams <- rownames(x) # store sample names for later
    x <- as.matrix(x) # convert to matrix for speed (?)
    dimnames(x) <- NULL # clear the dimnames for speed (?)
    .call <- match.call()
    method <- match.arg(method)
    dis <- distance(x, method = method) # calculate the distances
    Wmeans <- apply(dis, 2, cumWmean, y) # Estimated values
    means <- apply(dis, 2, cummean, y)
    minDC <- apply(dis, 2, minDij) # minimum Dij per sample
    error <- Werror <- matrix(ncol = length(y), nrow = nrow(Wmeans))
    for(i in seq(along = y))
      {
        Werror[, i] <- y[i] - Wmeans[, i] # residuals for Wmeans
        error[, i] <- y[i] - means[, i] # residuals for mean
      }
    WRMSE <- apply(Werror^2, 1, function(x) sqrt(mean(x))) # RMSE
    k.w <- which.min(WRMSE)
    RMSE <- apply(error^2, 1, function(x) sqrt(mean(x)))
    k <- which.min(RMSE)
    Wbias <- apply(Werror, 1, mean)  # average bias
    bias <- apply(error, 1, mean)
    Wmax.bias <- apply(Werror, 1, maxBias, y) # maximum bias
    max.bias <- apply(error, 1, maxBias, y)
    r2.mean <- apply(means, 1, function(x, y) {cor(x, y)^2}, y) # r.squared
    r2.Wmean <- apply(Wmeans, 1, function(x, y) {cor(x, y)^2}, y)
    ## re-apply samples names and n. closest
    colnames(Wmeans) <- colnames(means) <- site.nams
    colnames(Werror) <- colnames(error) <- site.nams
    rownames(Wmeans) <- rownames(means) <- 1:(dims[1] -1)
    rownames(Werror) <- rownames(error) <- 1:(dims[1] -1)
    ## return results
    structure(list(standard = list(est = means, resid = error,
                     rmse = RMSE, avg.bias = bias, max.bias = max.bias,
                     r.squared = r2.mean, k = k, auto = TRUE),
                   weighted = list(est = Wmeans, resid = Werror,
                     rmse = WRMSE, avg.bias = Wbias, max.bias = Wmax.bias,
                     r.squared = r2.Wmean, k = k.w, auto = TRUE),
                   Dij = dis,
                   orig.x = x,
                   orig.y = y,
                   call = .call,
                   method = method),
              class = "mat")
  }

print.mat <- function(x, k = 10,
                      digits = min(3, getOption("digits") - 4),
                      ...)
  {
    ##if(is.null(k))
    ##  k <- k(x)
    cat("\n")
    writeLines(strwrap("Modern Analogue Technique", prefix = "\t"))
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    tbl <- cbind(x$standard$rmse[1:k], x$standard$r.squared[1:k],
                 x$standard$avg.bias[1:k], x$standard$max.bias[1:k])
    tbl.w <- cbind(x$weighted$rmse[1:k], x$weighted$r.squared[1:k],
                   x$weighted$avg.bias[1:k], x$weighted$max.bias[1:k])
    tbl <- as.matrix(format(tbl, digits = digits))
    tbl.w <- as.matrix(format(tbl.w, digits = digits))
    tbl <- cbind(as.integer(1:k), tbl)
    tbl.w <- cbind(as.integer(1:k), tbl.w)
    rownames(tbl) <- rownames(tbl.w) <- rep("", nrow(tbl))
    colnames(tbl) <- colnames(tbl.w) <- c("k",
                                          "RMSE","R2","Avg Bias","Max Bias")
    cat("\nQuantiles of the dissimilarities for the training set:\n\n")
    print(quantile(as.dist(x$Dij), probs = c(0.01, 0.02, 0.05, 0.1, 0.2)),
          digits = digits)
    cat("\nInferences based on the mean of k-closest analogues:\n\n")
    print(tbl, quote = FALSE, right = TRUE)
    cat("\nInferences based on the weighted mean of k-closest analogues:\n\n")
    print(tbl.w, quote = FALSE, right = TRUE)
    cat("\n")
    invisible(x)
  }
