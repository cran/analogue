###########################################################################
##                                                                       ##
## Internal functions for package analogue - not meant to be used by     ##
## users.                                                                ##
##                                                                       ##
###########################################################################

###########################################################################
##                                                                       ##
## cumWmean - calculates the cumulative weighted mean of y               ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## weights           - the weights to use                                ##
## y                 - the vector of values to calculate weighted mean   ##
##                     of                                                ##
## drop              - drop spurious zero distance                       ##
##                                                                       ##
###########################################################################
cumWmean <- function(weights, y, drop = TRUE)
  {
    #if (length(weights) != length(y)) 
    #  stop("'y' and 'weights' must have the same length")
    ord <- order(weights)
    if(drop) {
      weights <- 1 / weights[ord][-1]
      env <- y[ord][-1]
    } else {
      weights <- 1 / weights[ord]
      env <- y[ord]
    }
    cumsum(weights * env) / cumsum(weights)
  }
###########################################################################
##                                                                       ##
## cummean - calculates the cumulative mean of y                         ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## dis               - the distances to sort by                          ##
## y                 - the vector of values to calculate mean of         ##
## drop              - drop spurious zero distance                       ##
##                                                                       ##
###########################################################################
cummean <- function(dis, y, drop = TRUE)
  {
    ord <- order(dis)
    if(drop) {
      dis <- dis[ord][-1]
      y <- y[ord][-1]
    } else {
      dis <- dis[ord]
      y <- y[ord]
    }
    cumsum(y) / 1:length(dis)
  }
###########################################################################
##                                                                       ##
## minDij - returns the non-zero minimum distance                        ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## x                 - the vector of distances for which the non-zero    ##
##                     minimum is required                               ##
##                                                                       ##
###########################################################################
minDij <- function(x)
  {
    ord <- order(x)
    x[ord][2] # we don't want the first zero distance
  }
###########################################################################
##                                                                       ##
## maxBias - returns the maximum bias statistic of mat residuals         ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## error             - model residuals                                   ##
## y                 - the vector original observed env data             ##
## n                 - number of sections to break env gradient into     ##
##                                                                       ##
###########################################################################
maxBias <- function(error, y, n = 10)
  {
    groups <- cut(y, breaks = n, labels = 1:n)
    bias <- aggregate(error, list(group = groups), mean)$x
    bias[which.max(abs(bias))]
  }
