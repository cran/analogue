\name{plot.bayesF}
\alias{plot.bayesF}
\concept{likelihood ratios}
\concept{Bayes factors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayes factor plots}
\description{
  Plotting function for \code{\link{bayesF}}, to draw the posterior
  probability of analogue.
}
\usage{
\method{plot}{bayesF}(x, xlab = NULL, ylab = "Pr (A+ | d)",
        col = "red", abline.col = "lightgrey", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \dQuote{\code{bayesF}}.}
  \item{xlab,ylab}{the x- and y-axis labels for the plot.}
  \item{col}{colour of the line used to draw the posterior probability.}
  \item{abline.col}{colour of the vertical line drawn to indicate the
    optimal dissimilarity determined from the ROC curve.}
  \item{\dots}{other plot arguments passed to plotting
    functions. Currently ignored.}
}
\details{
  This is a simple plotting function that currently only draws the
  posterior probability of analogue based on the Bayes facotr or
  likelihood ratio of a positive event (analogue).
}
\value{
  A plot on the current device.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{bayesF}} and \code{\link{roc}}}
\examples{
## continue the example from ?roc
example(roc)

## calculate the Bayes factors of analogue and no-analogue
## (uses observed probabilities of analogue/no-analogue
swap.bayes <- bayesF(swap.roc)
swap.bayes

## plot the probability of analogue
plot(swap.bayes)

## calculate the Bayes factors of analogue and no-analogue
## with prior probabilities c(0.5, 0.05)
swap.bayes2 <- bayesF(swap.roc, prior = c(0.5, 0.05))
swap.bayes

## plot the probability of analogue
plot(swap.bayes2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{methods}
