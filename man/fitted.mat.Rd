\name{fitted.mat}
\alias{fitted.mat}
\alias{print.fitted.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitted values of a MAT model}
\description{
  Extracts the fitted values for the training set of a MAT model using
  the, possibly weighted, mean of the environment for the
  \emph{k}-closest analogues.
}
\usage{
\method{fitted}{mat}(object, k, weighted = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{mat}}.}
  \item{k}{numeric; choose which of the \emph{k}-closest analogue
    models' fitted values is printed. Overides the default stored in
    the object returned by \code{fitted.mat}.}
  \item{weighted}{logical; should weighted averages be used instead of
    simple averages?}
  \item{\dots}{arguments to be passed to other methods.}
}
\details{
  The fitted values are for the training set and are taken as the,
  possibly weighted, mean of the environmental variable in question
  across the \emph{k}-closest analogues. The fitted value for each
  sample does \bold{not} include a contribution from itself --- it is
  the closest analogue, having zero dissimilarity. This spurious
  distance is ignored and analogues are ordered in terms of the non-zero
  distances to other samples in the training set, with the
  \emph{k}-closest contributing to the fitted value.
}
\value{
  A list with the following components:
  \item{estimated}{numeric; a vector of fitted values.}
  \item{k}{numeric; this is the \emph{k}-closest analogue model with
    lowest apparent RMSE.}
  \item{weighted}{logical; are the fitted values the weighted averages
    of the environment for the \emph{k}-closest analogues. If
    \code{FALSE}, the fitted values are the average of the environment
    for the \emph{k}-closest analogues.} 
}
\author{Gavin L. Simpson}
\seealso{\code{\link{mat}} for the model fitting function,
  \code{\link{residuals.mat}} for extraction of residuals. Related
  functions are \code{\link{predict.mat}} for predicting for unknown
  samples and \code{\link{bootstrap}} for bootstrap estimates for the
  training set.}
\examples{
## continue the RLGH and SWAP example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")
swap.mat

## the fitted values:
fitted(swap.mat)

## the fitted values for a MAT model
## with 10 analogues
fitted(swap.mat, k = 10)
}
\keyword{methods}% at least one, from doc/KEYWORDS
\keyword{multivariate}
