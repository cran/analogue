\name{reconPlot}
\alias{reconPlot}
\alias{reconPlot.default}
\alias{reconPlot.predict.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stratigraphic plots of palaeoenvironmental reconstructions}
\description{
  Draws a palaeoenvironmental reconstruction of predicted environmental
  values for sub-fossil assemblages.
}
\usage{

reconPlot(x, ...)

\method{reconPlot}{default}(x, ...)

\method{reconPlot}{predict.mat}(x, depths, use.labels = FALSE,
          predictions = c("apparent", "bootstrap"),
          error.bars = FALSE, sample.specific = TRUE,
          rev.x = TRUE, type = "l", xlim, ylim,
          xlab = "", ylab = "", main = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An R object. Currently only objects of classes
    \code{"predict.mat"} and \code{"mat"} are supported.}
  \item{depths}{numeric; a vector of depths for which predicted values
    exist or will be generated. Can be missing, in which case,
    \strong{if} \code{use.labels = TRUE}, the function will attempt to
    derive suitable values for you. See Details below.}
  \item{use.labels}{logical; should \code{reconPlot} attempt to derive
    values for argument \code{depths} from the names of the predicted
    values? Only use if \code{depths} is missing. See Details below.}
  \item{predictions}{character; one of \code{"apparent"} or
    \code{"bootstrap"}. Which type of predicted values should be
    plotted? The actual predictions (\code{"apparent"}) or the
    bootstrap-derived predictions (\code{"bootstrap"}).}
  \item{error.bars}{logical; should error bars be drawn on the plot?}
  \item{sample.specific}{logical; should sample specific errors be used?
  Only for \code{predictions = "bootstrap"}.}
  \item{rev.x}{logical; should the depth/age axis be reversed (drawn
    from high to low)?}
  \item{type}{type of line drawn. See \code{\link{par}} and argument
    \code{"type"}.}
  \item{xlab, ylab}{character; the x- and y-axis labels respectively.}
  \item{main}{character; main title for the plot.}
  \item{xlim, ylim}{numeric, length 2; the x- and y-limits for the
    plotted axes. If not provided, the function will calculate
    appropriate values to cover the range of plotted values and any
    error bars (of requested via \code{"error.bars"}.}
  \item{\dots}{arguments passed to other graphics functions.}
}
\details{
  Conventionally, these plots are drawn on a depth or an age
  scale. Argument \code{depths} is used to provide the depth or age
  axis, against which the predicted values are plotted.

  If \code{depths} is not provided, then the function will try to
  derive the appropriate values from the labels of the predictions if
  \code{use.labels = TRUE}. You must provide \code{depths} or set
  \code{use.labels = TRUE} otherwise an error will result. The derived
  labels will be coerced to numerics. If your labels are coercible, then
  you'll either get nonsense on the plot or an error from R. If so,
  provide suitable values for \code{depths}.
}
\value{
  A plot on the currently active device.
}
\author{Gavin L. Simpson }
\seealso{\code{\link{mat}}, and \code{\link{predict.mat}} }
\examples{
## continue the RLGH example from ?mat
example(mat)

## draw the reconstruction
reconPlot(rlgh.mat, use.labels = TRUE, error.bars = TRUE,
          xlab = "Depth", ylab = "pH")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
