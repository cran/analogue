\name{residuals.mat}
\alias{residuals.mat}
\alias{resid.mat}
\alias{print.residuals.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract Model Residuals}
\description{
  \code{\link{residuals}} is a is a generic function which extracts
  model residuals from objects returned by modeling
  functions. \code{residuals.mat} is a \code{\link{residuals}} method
  for \code{\link{mat}} models.

  The abbreviated form of this function is \code{resid}.
}

\usage{
\method{residuals}{mat}(object, k, weighted = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{mat}}.}
  \item{k}{number of analogues to use. If missing, \code{k} is chosen
    automatically as the \code{k} that achieves lowest RMSE.}
  \item{weighted}{logical; should the weighted mean of the environment
    for the \code{"k"} modern analogues be used instead of the mean?}
  \item{\dots}{arguments pass to or from other methods.}
}
\details{
  Typical usages are
  \preformatted{
    resid(object, k, weighted = FALSE, \dots)}
  }
\value{
  A list with the following components:
  \item{residuals}{numeric; a vector of model residuals.}
  \item{k}{numeric; either the user supplied number of analogues used,
    or the \emph{k}-closest analogue model with lowest apparent RMSE.}
  \item{weighted}{logical; are the residuals from a model where fitted
    values are the weighted averages of the environment for the
    \emph{k}-closest analogues. If \code{FALSE}, the residuals are from
    a model using the average of the environment for the
    \emph{k}-closest analogues.} 
}
\author{Gavin L. Simpson}
\seealso{\code{\link{residuals}}, \code{\link{mat}}.}
\examples{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## model residuals
resid(swap.mat) # uses abbreviated form

## model residuals using 10 analogues and
## weighted averages
resid(swap.mat, k = 10, weighted = TRUE)
}
\keyword{methods}% at least one, from doc/KEYWORDS
