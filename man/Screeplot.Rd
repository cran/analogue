\name{Screeplot}
\alias{Screeplot}
\alias{Screeplot.mat}
\alias{Screeplot.default}
%\alias{Screeplot.bootstrap.mat}
\alias{Screeplot.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Screeplots of model results}
\description{
  Draws screeplots of performance statistics for models of varying complexity.
}
\usage{
Screeplot(x, \ldots)

\method{Screeplot}{default}(x, \ldots)

\method{Screeplot}{mat}(x, k, restrict = 20,
          display = c("rmse", "rmsep", "avg.bias",
                      "max.bias", "r.squared"),
          weighted = TRUE,  xlab = NULL, ylab = NULL,
          main = NULL, sub = NULL, \ldots)

\method{Screeplot}{bootstrap}(x, k, restrict = 20,
          display = c("rmse","rmsep","avg.bias","max.bias",
                      "r.squared"),
          legend = TRUE, loc.legend = "topright",
          xlab = NULL, ylab = NULL, main = NULL, sub = NULL,
          ..., lty = c("solid","dashed"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object on which method dispatch applied; currently only for
    class \code{\link{mat}} and \code{bootstrap}.}
  \item{k}{number of analogues to use. If missing 'k' is chosen
    automatically as the 'k' that achieves lowest RMSE.} 
  \item{restrict}{logical; restrict comparison of k-closest model to k
    \eqn{<=}{<=} \code{restrict}.}
  \item{display}{which aspect of \code{x} to plot? Partial match.}
  \item{weighted}{logical; should the analysis use weighted mean of env
    data of analogues as fitted/estimated values?}
  \item{xlab, ylab}{x- and y-axis labels respectively.}
  \item{main, sub}{main and subtitle for the plot.}
  \item{legend}{logical; should a legend be displayed on the figure?}
  \item{loc.legend}{character; a keyword for the location of the
    legend. See \code{\link{legend}} for details of allowed keywords.}
  \item{lty}{vector detailing the line type to use in drawing the
    screeplot of the apparent and bootstrap statistics,
    respectively. Code currently assumes that \code{length(lty)} is 2.}
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

  For the \code{\link{bootstrap}} method, apparent and bootstrap
  versions of these statistics are available and plotted.
}
\author{Gavin Simpson}
\note{
  Currently only models of class \code{\link{mat}} and
  \code{\link{bootstrap}} are supported.
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
