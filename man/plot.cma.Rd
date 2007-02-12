\name{plot.cma}
\alias{plot.cma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draws stripcharts of close modern analogue results}
\description{
  Draws stripcharts of close modern analogue
  results. \code{\link{stripchart}} produces one dimensional scatter
  plots (or dot plots) of the given data.  These plots are a good
  alternative to \code{\link{boxplot}}'s when sample sizes are small.
}
\usage{
\method{plot}{cma}(x, method = c("overplot", "jitter", "stack"),
   jitter = 0.1, vertical = FALSE,
   draw.quant = TRUE, xlab = NULL, ylab = "",
   main = "", cex.axis = NULL, ...,
   col.quant = "red", lty.quant= "dashed")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"cma"}.}
  \item{method}{the method to be used to separate coincident points. The
    default method\cr\code{"overplot"} causes such points to be
    overplotted, but it is also possible to specify \code{"jitter"} to 
    jitter the points, or \code{"stack"} have coincident points
    stacked.  The last method only makes sense for very granular data.}
  \item{jitter}{when \code{method="jitter"} is used, \code{jitter} gives
    the amount of jittering applied.}
  \item{vertical}{when vertical is \code{TRUE} the plots are drawn
    vertically rather than the default horizontal.}
  \item{draw.quant}{logical; should the quantiles be drawn on the stripchart?}
  \item{xlab,ylab,main}{Graphical parameters}
  \item{cex.axis}{The magnification to be used for axis annotation
    relative to the current setting of \code{cex}. See
    \code{\link{par}}.}
  \item{\dots}{additional arguments passed to \code{\link{stripchart}}}
  \item{col.quant,lty.quant}{colour and line type in which to drawn the
    quantile lines.}
}
\details{
  Simple wrapper to \code{\link{stripchart}}.
}
\value{
  A plot is produced on the currently active device. Also invisibly
  returns a list with the following components:
  \item{distances}{a vector of stacked distances extracted from
    \code{object}.}
  \item{groups}{a factor listing the fossil sample for which the
    distances are the distances to the close modern analogues for the
    training set.}
}
\author{Gavin L. Simpson }
\seealso{\code{\link{stripchart}}, or \code{\link{boxplot}} for an
  alternative representation.}
\examples{
## continue the RLGH and SWAP example from ?join
example(join)

## analog matching between SWAP and RLGH core
swap.analog <- analog(swapdiat, rlgh, method = "chord")
swap.analog
summary(swap.analog)

## close modern analogues
swap.cma <- cma(swap.analog)

## plot the results
## many samples so reduce the label and point size
plot(swap.cma, cex.axis = 0.6, cex = 0.6)
}
\keyword{hplot}% at least one, from doc/KEYWORDS
