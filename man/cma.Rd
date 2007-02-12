\name{cma}
\alias{cma}
\alias{cma.default}
\alias{cma.analog}
\alias{print.cma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Close modern analogues}
\description{
  Extracts and formats close modern analogue samples from a modern
  reference set that are closer than a defined cut off threshold.
}
\usage{
cma(object, \dots)

\method{cma}{default}(object, \dots)

\method{cma}{analog}(object, cutoff, prob = c(0.01, 0.025, 0.05), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object for which close modern analogues are to be
    returned. Currently only for objects of class \code{\link{analog}}.}
  \item{cutoff}{numeric; critical value determining level below which
    samples from the modern reference set are defined as close modern
    analogues. May be missing, in which case the 2.5\% quantile of the
    training set dissimilarities is used unless \code{object$train} is
    \code{NULL}, in which case \code{"cutoff"} must be supplied.}
  \item{prob}{numeric vector of probabilities with values in [0,1], for
    which quantiles of the distribution of training set dissimilarities
    will be calculated. See \code{\link{quantile}}.}
  \item{\dots}{arguments to be passed to other \code{\link{cma}} methods.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of class \code{"cma"} with the following components:
  \item{distances}{a matrix of distances/dissimilarities. Individual
    columns contain the ordered close modern analogues for individual
    fossil samples. Rows of this matrix refer to the
    \eqn{k^{\mathrm{th}}}{k-th} closest analogue for each fossil
      sample. See notes below.}
  \item{samples }{a matrix of sample names from the reference set that
    are close modern analogues for a fossil sample. Individual
    columns contain the ordered close modern analogues for individual 
    fossil samples. Rows of this matrix refer to the
    \eqn{k^{\mathrm{th}}}{k-th} closest analogue for each fossil sample.
    See notes below.}
  \item{cutoff}{the cutoff threshold used to define close modern
    analogues.}
  \item{quant}{numeric vector of the requested quantiles.}
  \item{probs}{the probabilities of the requested quantiles.}
  \item{method}{character; the dissimilarity coefficient used}
  \item{n.analogs}{numeric vector of the number of analogues per fossil
    sample.}
}
\references{
  
  Flower, R.J., Juggins, S. and Battarbee, R.W. (1997) Matching diatom
  assemblages in lake sediment cores and modern surface sediment
  samples: the implications for lake conservation and restoration with
  special reference to acidified systems. \emph{Hydrobiologia}
  \strong{344}; 27--40.

  Simpson, G.L., Shilland, E.M., Winterbottom, J. M. and Keay, J. (2005)
  Defining reference conditions for acidified waters using a modern
  analogue approach. \emph{Environmental Pollution} \strong{137};
  119--133.
}
\author{Gavin L. Simpson}
\note{
  Currently, only objects of class \code{\link{analog}} are supported.
  The number of rows in the returned matrices is equal to the maximum
  number of close modern analogues identified for an individual fossil
  sample. If no close modern analogues exist for an individual fossil
  sample, then the relevant column in \code{"distances"}will contain all
  missing values and in \code{"samples"} the string \code{"none"}. Rows
  of individual columns will be padded with missing values if the number
  of close modern analogues for that sample is less than the maximum
  number of close modern analogues identified for a single sample.
}
\section{Warning}{
  The matrix representation of the returned objects \code{"distances"}
  and \code{"samples"} is used for convenience \bold{only}. Do not
  attempt to compare values across rows between columns.
}
\seealso{\code{\link{analog}}}
\examples{
## continue the RLGH example from ?join
example(join)

## analog matching between SWAP and RLGH core
swap.analog <- analog(swapdiat, rlgh, method = "chord")
swap.analog
summary(swap.analog)

## close modern analogues
swap.cma <- cma(swap.analog)
swap.cma
summary(swap.cma)

}
\keyword{methods}% at least one, from doc/KEYWORDS
\keyword{manip}
