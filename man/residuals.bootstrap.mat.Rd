\name{residuals.bootstrap.mat}
\alias{residuals.bootstrap.mat}
\alias{resid.bootstrap.mat}
\alias{print.residuals.bootstrap.mat}
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
\method{residuals}{bootstrap.mat}(object, which = c("apparent", "bootstrap"), ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{mat}}.}
  \item{which}{character; which set of residuals to return, the
    apparanet residuals or the residuals of the bootstrap-derived
    estimates?}
  \item{\dots}{arguments pass to or from other methods.}
}
\value{
  A list containg the requested residuals and metadata, with the
  following components:
  \item{apparent}{residuals for the MAT-estimated model.}
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
\seealso{\code{\link{residuals}}, \code{\link{mat}}.}
\examples{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## bootstrap training set
swap.boot <- bootstrap(swap.mat, k = 10, n.boot = 1000)

## residuals
resid(swap.boot) # uses abbreviated form
}
\keyword{methods}% at least one, from doc/KEYWORDS
