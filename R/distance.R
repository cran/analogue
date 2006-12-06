###########################################################################
##                                                                       ##
## distance - function to compute distances between samples              ##
##                                                                       ##
## Created       : 17-Apr-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.0-1                                                 ##
## Last modified : 05-Nov-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
##                                                                       ##
###########################################################################
## x = training data, y = fossil data
distance <- function(x, y,
                     method = c("euclidean", "SQeuclidean", "chord",
                       "SQchord", "bray", "chi.square", "SQchi.square",
                       "information", "chi.distance", "manhattan",
                       "kendall", "gower"))
  {
    euclidean <- function(x, y)
      {
        sqrt(sum((x - y)^2))
      }
    SQeuclidean <- function(x, y)
      {
        sum((x - y)^2)
      }
    chord <- function(x, y)
      {
        x <- sqrt(x); y <- sqrt(y)
        euclidean(x, y)
      }
    SQchord <- function(x, y)
      {
        x <- sqrt(x); y <- sqrt(y)
        SQeuclidean(x, y)
      }
    bray <- function(x, y)
      {
        sum(abs(x - y)) / sum(x + y)
      }
    chi.square <- function(x, y)
      {
        inds <- !(x == 0 & y == 0)
        sqrt(sum(((x[inds] - y[inds])^2) / (x[inds] + y[inds])))
      }
    SQchi.square <- function(x, y)
      {
        inds <- !(x == 0 & y == 0)
        sum(((x[inds] - y[inds])^2) / (x[inds] + y[inds]))
      }
    information <- function(x, y)
      {
        XY <- x + y
        A <- x * log2((2 * x) / XY)
        B <- y * log2((2 * y) / XY)
        sum(A, na.rm = TRUE) + sum(B, na.rm = TRUE)
      }
    chi.distance <- function(x, y, colsum)
      {
        sqrt(sum(((x - y)^2) / (colsum / sum(colsum))))
      }
    manhattan <- function(x, y)
      {
        sum(abs(x - y))
      }
    kendall <- function(x, y, maxi)
      {
        sum(maxi - min(x, y))
      }
    gower <- function(x, y, maxi, mini)
      {
        sqrt( 2 * sum(abs(x - y) / (maxi - mini)))
      }
    Dist <- function(y, x, method, ...)#colsum = NULL)
      {
        switch(method,
               euclidean = apply(x, 1, euclidean, y),
               SQeuclidean = apply(x, 1, SQeuclidean, y),
               chord = apply(x, 1, chord, y),
               SQchord = apply(x, 1, SQchord, y),
               bray = apply(x, 1, bray, y),
               chi.square = apply(x, 1, chi.square, y),
               SQchi.square = apply(x, 1, SQchi.square, y),
               information = apply(x, 1, information, y),
               chi.distance = apply(x, 1, chi.distance, y, colsum = colsum),
               manhattan = apply(x, 1, manhattan, y),
               kendall = apply(x, 1, kendall, y, maxi = maxi),
               gower = apply(x, 1, gower, y, maxi = maxi, mini = mini)
               )
      }
    x <- as.matrix(x)
    x.names <- rownames(x)
    if(missing(y)) {
      colsumx <- colSums(x)
      if(any(colsumx <= 0)) {
        x <- x[, colsumx > 0, drop = FALSE]
        warning("some species contain no data and were removed from data matrix x\n")
      }
      y <- x
      y.names <- x.names
    } else {
      y <- as.matrix(y)
      y.names <- rownames(y)
    }
    if(missing(method))
      method <- "euclidean"
    method <- match.arg(method)
    if(method == "chi.distance")
      colsum <- colSums(join(as.data.frame(x),as.data.frame(y)))
    if(method %in% c("kendall", "gower")){
      maxi <- apply(join(as.data.frame(x),as.data.frame(y)), 2, max)
      mini <- apply(join(as.data.frame(x),as.data.frame(y)), 2, min)
    }
    dimnames(x) <- dimnames(y) <- NULL
    if(method == "chi.distance") {
      y <- y / rowSums(y)
      x <- x / rowSums(x)
      res <- apply(y, 1, Dist, x, method, colsum = colsum)
    } else if (method == "kendall") {
      res <- apply(y, 1, Dist, x, method, maxi = maxi)
    } else if (method == "gower") {
      res <- apply(y, 1, Dist, x, method, maxi = maxi, mini = mini)
    } else{
      res <- apply(y, 1, Dist, x, method)
    }
    colnames(res) <- y.names
    rownames(res) <- x.names
    return(res)
  }
