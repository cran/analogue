\name{summary.predict.mat}
\alias{summary.predict.mat}
\alias{print.summary.predict.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarising MAT model predictions}
\description{
  \code{\link{summary}} method for objects of class \code{"predict.mat"}.
}
\usage{
\method{summary}{predict.mat}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"predict.mat"}, usually the
    result of a call to \code{\link{predict.mat}}.}
  \item{\dots}{arguments passed to or from other methods.}
}
\details{
  Currently a simple wrapper for \code{\link{summary.bootstrap.mat}}.
}
\value{
  An object of class \code{"summary.bootstrap.mat"}, see
  \code{\link{summary.bootstrap.mat}} for more details.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{predict.mat}}, \code{\link{mat}},
  \code{\link{bootstrap}} and \code{\link{summary}}.}
\examples{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## predict for RLGH data
swap.pred <- predict(swap.mat, rlgh, bootstrap = FALSE)
summary(swap.pred)

## again but with bootstrap estimates and sample specific errors
swap.pred <- predict(swap.mat, rlgh)
summary(swap.pred)
}
\keyword{methods}% at least one, from doc/KEYWORDS
