\name{Screeplot}
\alias{Screeplot}
\alias{Screeplot.mat}
\alias{Screeplot.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Screeplots of model results}
\description{
  Draws screeplots of performance statistics for models of varying complexity.
}
\usage{
Screeplot(x, \ldots)

\method{Screeplot}{default}(x, \ldots)

\method{Screeplot}{mat}(x,
          k,
          restrict = 20,
          display = c("rmse", "rmsep", "avg.bias",
                      "max.bias", "r.squared"),
          weighted = TRUE,
          xlab = NULL,
          ylab = NULL,
          main = NULL,
          sub = NULL,
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object on which method dispatch applied; currently only for
    class \code{\link{mat}}.}
  \item{k}{number of analogues to use. If missing 'k' is chosen
    automatically as the 'k' that achieves lowest RMSE.} 
  \item{restrict}{logical; restrict comparison of k-closest model to k
    \eqn{<=}{<=} \code{restrict}.}
  \item{display}{which aspect of \code{\link{mat}} object to
    plot? Partial match.}
  \item{weighted}{logical; should the analysis use weighted mean of env
    data of analogues as fitted/estimated values?}
  \item{xlab, ylab}{x- and y-axis labels respectively.}
  \item{main, sub}{main and subtitle for the plot.}
  \item{\dots}{arguments passed to other graphics functions.}
}
\details{
  Screeplots are often used to graphically show the results of
  cross-validation or other estimate of model performance across a range
  of model complexity.

  Five measures of model performance are currently available: i) root
  mean square error (RMSE); ii) root mean square error of prediction
  (RMSEP); iii) average bias --- the mean of the model residuals; iv)
  maximum bias --- the maximum average bias calculated for each of
  \emph{n} sections of the gradient of the environmental variable; and
  v) model \eqn{R^2}.
}
\author{Gavin Simpson}
\note{
  Currently only models of class \code{\link{mat}} supported.
}
\seealso{\code{\link{screeplot}} }
\examples{
## load the example data
data(swapdiat)
data(swappH)
data(rlgh)
## process so common set of columns for training and test
## number of training set samples
n.train <- nrow(swapdiat)
## merge training and test set on columns
dat <- join(swapdiat, rlgh, verbose = TRUE)
## convert to proportions
dat <- dat / 100
## subset data back into training and test sets
swapdiat <- dat[1:n.train, ]
rlgh <- dat[(n.train+1):nrow(dat), ]
## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")
swap.mat

##
Screeplot(swap.mat)
}
\keyword{hplot}
\keyword{methods}
