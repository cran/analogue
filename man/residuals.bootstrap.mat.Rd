\name{residuals.bootstrap.mat}
%\name{residuals.bootstrap}
\alias{residuals.bootstrap.mat}
%\alias{residuals.bootstrap}
\alias{resid.bootstrap.mat}
%\alias{resid.bootstrap}
\alias{print.residuals.bootstrap.mat}
%\alias{print.residuals.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract Bootstrapped Model Residuals}
\description{
  \code{\link{residuals}} is a generic function which extracts
  model residuals from objects returned by modeling
  functions. \code{residuals.bootstrap.mat} is a \code{\link{residuals}} method
  for \code{\link{mat}} models.

  The abbreviated form of this function is \code{resid}.
}
\usage{
\method{residuals}{bootstrap.mat}(object, which = c("model", "bootstrap"), \dots) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{bootstrap.mat}}.}
  \item{which}{character; which set of residuals to return, the
    model residuals or the residuals of the bootstrap-derived
    estimates?}
  \item{\dots}{arguments pass to or from other methods.}
}
\details{
  Typical usages are
  \preformatted{
    resid(object, which = c("model", "bootstrap"), \dots)}
  }
\value{
  A list containg the requested residuals and metadata, with the
  following components:
  \item{model}{Leave one out residuals for the MAT-estimated model.}
  \item{bootstrap}{residuals for the bootstrapped MAT model.}
  \item{k}{numeric; indicating the size of model used in estimates and
    predictions.}
  \item{n.boot}{numeric; the number of bootstrap samples taken.}
  \item{auto}{logical; whether \code{"k"} was choosen automatically or
    user-selected.}
  \item{weighted}{logical; whether the weighted mean was used instead of
    the mean of the environment for \emph{k}-closest analogues.}
}
\author{Gavin L. Simpson}
\seealso{\code{\link{residuals}}, \code{\link{mat}},
  \code{\link{bootstrap.mat}}.}
\examples{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## bootstrap training set
swap.boot <- bootstrap(swap.mat, k = 10, n.boot = 100)

## residuals
resid(swap.boot) # uses abbreviated form
}
\keyword{methods}% at least one, from doc/KEYWORDS
