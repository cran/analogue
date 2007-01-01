\name{k}
\alias{k}
\alias{k.default}
\alias{k.mat}
\alias{k<-}
\alias{k<-.default}
\alias{k<-.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract and set the number of analogues}
\description{
  An extractor function to access the number of analogues used in
  particular models. Can also be used to set or update the stored value
  of \eqn{k}.
}
\usage{
k(object, ...)

\method{k}{mat}(object, weighted=FALSE, ...)

k(object, weighted=FALSE) <- value

\method{k}{mat}(object, weighted=FALSE) <- value

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object; currently only for objects of class
    \code{\link{mat}}.}
  \item{weighted}{logical; extract/set number of analogues for a
    weighted or un-weighted model?}
  \item{\dots}{further arguments to other methods.}
  \item{value}{integer; replacement value for \eqn{k}.}
}
\details{
  \code{k} is a generic accessor function, and \code{k<-} is a generic
  replacement function.
}
\value{
  For \code{k}, an integer value that is the number of analogues stored
  for use. The returned object has attributes \dQuote{auto} and
  \dQuote{weighted}. \dQuote{auto} refers to whether the extracted value
  of \eqn{k} was set automatically (\code{TRUE}) or by the user
  (\code{FALSE}). \dQuote{weighted} states if the returned value is for
  a \code{weighted} analysis or an un-\code{weighted} analysis (\code{FALSE}).

  For \code{k<-}, the updated object.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{mat}}}
\examples{
## continue the example from join
example(join)

## fit a MAT model
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## How many analogues gives lowest RMSE?
k(swap.mat)
## note that this value was chosen automatically

## Now set k to be 10
k(swap.mat) <- 10

## check
k(swap.mat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{manip}% __ONLY ONE__ keyword per line
