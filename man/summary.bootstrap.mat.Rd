%\name{summary.bootstrap.mat}
%\alias{summary.bootstrap.mat}
%\alias{print.summary.bootstrap.mat}
\name{summary.bootstrap}
\alias{summary.bootstrap}
\alias{print.summary.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarising bootstrap resampling for MAT models}
\description{
  \code{\link{summary}} method for class \code{"bootstrap"}.
}
\usage{
\method{summary}{bootstrap}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"bootstrap"}, usually the
    result of a call to \code{\link{bootstrap}}.}
  \item{\dots}{arguments passed to or from other methods.}
}
\value{
  A data frame with the following components:
  
  \item{Obs}{observed response value for training set samples.}
  \item{Est}{estimated (fitted) values from the \code{\link{mat}}
    model.}
  \item{Resid}{residuals from the \code{\link{mat}} model.}
  \item{Boot.Est}{bootstrap estimated (fitted) values for the
    \code{\link{mat}} model.}
  \item{Boot.Resid}{residuals of the bootstrap estimates of a
    \code{\link{mat}} model.}
  \item{s1}{bootstrap-derived s1 error component for each sample in the
    training set.}
  \item{s2}{bootstrap-derived s2 error component for each sample in the
    training set.}
  \item{RMSEP}{bootstrap-derived RMSEP for each sample in the
    training set.}
}
\author{Gavin L. Simpson}
\seealso{\code{\link{bootstrap}}, \code{\link{mat}},
  \code{\link{summary}}.}
\examples{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## bootstrap training set
swap.boot <- bootstrap(swap.mat, k = 10, n.boot = 1000)
swap.boot
summary(swap.boot)

## bootstrap with predictions:
rlgh.boot <- bootstrap(swap.mat, rlgh, k = 10, n.boot = 1000)
summary(rlgh.boot)
}
\keyword{methods}% at least one, from doc/KEYWORDS
