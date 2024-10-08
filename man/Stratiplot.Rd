\name{Stratiplot}
\alias{Stratiplot}
\alias{Stratiplot.default}
\alias{Stratiplot.formula}
\alias{Stratiplot.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Palaeoecological stratigraphic diagrams}
\description{
  Draws palaeoecological stratigraphic diagrams of one or more variables
  as a function of depth/age, with the time dimension flowing from the
  bottom to the top of the y-axis, using the \pkg{Lattice} graphics
  package.
}
\usage{
Stratiplot(x, \dots)

\method{Stratiplot}{default}(x, y, type = "l", ylab = NULL, xlab = "",
           pages = 1, rev = TRUE, ylim, sort = c("none", "wa", "var"),
           svar = NULL, rev.sort = FALSE, strip = FALSE, topPad =6,
           varTypes = "relative", absoluteSize = 0.5,
           zoneNames = NULL, drawLegend = TRUE, na.action = "na.omit",
           labelValues = NULL, labelAt = NULL, labelRot = 60, yticks,
           \dots)

\method{Stratiplot}{formula}(formula, data, subset,
           na.action = "na.pass", type = "l",
           ylab = NULL, xlab = "", pages = 1, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix-like object; the variables to be plotted.}
  \item{y}{numeric vector of depths/ages corresponding to rows in
    \code{x}. Length of \code{y} must be the same as \code{nrow(x)} or
    exactly equal to \code{nrow(x) / ncol(x)}. See Details.}
  \item{formula}{an object of class \code{"formula"} (or one that can be
    coerced to that class): a symbolic description of the model to be
    fitted. The details of plot specification are given under
    \sQuote{Details}.}
  \item{type}{character; The type of plotting. Can be a vector. Note
    that not all Lattice \sQuote{type}s are supported and some new types
    are allowed. See \code{\link{panel.Stratiplot}} for further
    details.}
  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables to plot. If not found in \code{data}, the variables
    are taken from \code{environment(formula)}, typically the
    environment from which \code{Stratiplot} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be  used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s. The default is \code{"na.omit"} for the
    default method, which strips \code{NA}s from the stacked data,
    whereas the default for the formula method is \code{"na.pass"} which
    results in \code{NA} being passed on to the plotting function. See
    Details for further information.}
  \item{ylab, xlab}{the x- and y-axis labels.}
  \item{pages}{numeric; the number of pages to draw the plot over. May
    be useful for data sets with many species.}
  \item{rev}{logical; should the y-axis limits be reversed}
  \item{ylim}{user supplied limits for the y-axis (time/depth). If not
    supplied, suitable limits will be determined from the data. As such,
    in general use \code{ylim} need not be supplied. If you choose to
    supply your own \code{ylim} values, note the default for argument
    \code{rev}; the default will reverse the values you supply to
    \code{ylim}.}
  \item{sort}{character; how should the variables (columns) of \code{x}
    be sorted on the plot. \code{"wa"} sorts by weighted averages
    of variable \code{svar} if not \code{NULL} or of \code{y}
    otherwise. The default when \code{"wa"} is specified is to
    order by wiehgted average of the depth/time axis -- \code{y}. If
    \code{"var"}, then ordering is done as per the \strong{order} of
    \code{svar}.}
  \item{svar}{vector; optional variable to sort columns of \code{x} by.}
  \item{rev.sort}{logical; should the sorting order be reversed.}
  \item{strip}{logical; Should panels have strips containing variable
    labels drawn on them? Default is \code{FALSE}, which labels each
    panel with a label resulting in a more familiar plot style.}
  \item{topPad}{numeric; additional padding for the top axis to
    accomodate long variable names. This is a temporary fudge until the
    actual space required can be automagically calculated from the
    variable names themselves. The currently gets most of the way there,
    but \code{topPad} is used to add some extra space if required.}
  \item{varTypes}{a character vector of length 1 or equal in length to
    the number of variables plotted. If length 1, the vector is expanded
    to the required length. Two values are allowed; i.
    \code{"relative"}, and ii. \code{"absolute"}. \code{"relative"}
    treats the indicated variable as a percentage type variable and the
    panel width is scaled relative to the percentage type variables
    plotted. \code{"absolute"} treats the indicated variable as an
    absolute variable whose panel width should be independent of the
    other panels. Use \code{"absolute"} for variables that are not
    species compositional data or other percentage type data.}
  \item{absoluteSize}{numeric, length 1. This controls the width of
    panels for variables marked as \code{"absolute"}, and is the
    proportion of the largest non-\code{"absolute"} panel.}
  \item{zoneNames}{character vector of labels, one per zone, with which
    to label the zone legend, if drawn (see argument
    \code{drawLegend}). See Details.} 
  \item{drawLegend}{logical; should a legend for the zones}
  \item{labelValues}{a vector of labels for the variables
    plotted. Should be equal in length to the number or variables in the
    resulting plot. The main use for \code{labelValues} is to provide
    non-standard labels for the variables, including a vector of
    expressions. See Examples for an illustration of this.}
  \item{labelAt, labelRot}{these control the placement and rotation,
    respectively, of the variable labels. \code{labelAt} is the
    coordinate at which the label is drawn; currently only one value is
    used so you can't place labels in different locations depending on
    which panel is drawn. This will be fixed in a future version. The
    default location for the label is the panel mid-point.
    \code{labelAt} controls the rotation of the label; it is a numeric
    value in degree.}
  \item{yticks}{This is passed to the \code{scales} argument of
    \code{\link[lattice]{xyplot}} as component \code{at}. This should be a
    numeric vector of tick locations for the y (depth/age) axis. Setting
    this to \code{NULL} or \code{FALSE} suppresses ticks on the y
    axis. The default uses \code{TRUE}, which uses the default choices
    for \code{at} used by \code{\link[lattice]{xyplot}}.}
  \item{\dots}{additional arguments passed to
    \code{\link{panel.Stratiplot}} and the underlying
    \code{\link[lattice]{xyplot}} function.}
}
\details{
  The function now includes preliminary code to handle both relative
  (proportional or percentage data) and absolute data types, and
  mixtures thereof. Mixtures can be specified by supplying a vector of
  types to \code{varTypes}, in the same order as the variables are drawn
  on the plot.

  Plots can be specified symbolically using a formula. A typical model
  has the form \code{Y ~ variables}, where \code{Y} is either the core
  depths or sample ages/dates (to be plotted on the y-axis) and
  \code{variables} is a series of terms which specifies the variables to
  plot against \code{Y}. Terms should be specified with the form
  \code{var1 + var2 + var3} to plot only those variables. Other,
  standard, notation for formulae apply, such as model formulae used in
  \code{\link{lm}}.

  For the formula method the default for argument \code{na.action} is
  \code{"na.pass"}, which results in any \code{NA} values being passed
  on to the plotting code. This allows for plotting of proxies that been
  measured on different levels of the stratigraphy. Should you wish to
  have \code{NA} removed from the data before plotting, use
  \code{na.action = "na.omit"}, though do note this will remove all rows
  where any column/variable takes the value \code{NA}. The default
  \code{Stratiplot} method, which is used by the formula method for
  plotting, will strip any \code{NA} values from the data provided to
  it. This allows the function to correctly handle the situation where
  proxies are measured on different levels of the core \emph{and} you
  are displaying the data using lines of polygons. If the \code{NA} were
  not dropped by \code{Stratiplot.default}, any \code{NA} values would
  show up as breaks in the line or polygon drawn for each panel.

  In summary, the two methods have different defaults for
  \code{na.action} to allow them to handle proxies measured on different
  levels of the same core. This does mean that you can not use the
  formula interface \strong{and} strip NA's at the
  \code{Stratiplot.default} level. If you need that level of control use
  \code{Stratiplot.default} directly by not providing a formula as
  argument \code{x} and by supplying data for the y-axis via argument
  \code{y}. See Examples for an illustration of these features.

  Note that \code{formula} is \strong{not} passed on to
  \code{\link[lattice]{xyplot}}. Instead, the formula is parsed and
  evaluated within \code{Stratiplot} and an appropriate data structure
  formed to facilitate plotting via \code{\link[lattice]{xyplot}}. As
  such, the special features of \pkg{Lattice} formulae cannot be used.

  If zones are drawn on the stratigraphic plot, the \code{zoneNames}
  argument can be used to supply a set of names with which to label the
  zones using a legend. This legend is drawn on the right-hand side of
  the the straigraphic diagram if \code{drawLegend = TRUE} is
  supplied. The \code{zoneNames} must be supplied in stratigraphic
  order, as that is the order in which they are drawn in the
  legend. Whether this ordering is reversed or not will depend on the
  value of argument \code{rev}. It is up to the user to provide the
  labels in the correct order. Zones are specified by the zone
  boundaries (excluding the core sequence top and bottom), and as a
  result 1 more label is required than the number of zone boundaries
  supplied. If no \code{zoneNames} is not supplied, but a legend is
  requested, suitable names will be produced. If you do not wish to have
  any labelling at all, then set \code{zoneNames = ""} as this will get
  recycled to the correct length. See the Example section for an
  illustration of how this drawing zones works. 
}
\value{
  A side effect of calling \code{Stratiplot} is that a plot is drawn on
  the currently active device. A Lattice plot object of class
  \code{"trellis"} is returned invisibly. This is a change from pre
  0.17-0 version of the package.
}
\author{Gavin L. Simpson.}
\note{
  The function currently doesn't know about ages/dates and will
  interpret these as \sQuote{depths} instead. This will be fixed in a
  future version.
}
\seealso{\code{\link[lattice]{xyplot}},
\code{\link{panel.Stratiplot}}, \code{\link{panel.Loess}}.}
\examples{
data(V12.122)
Depths <- as.numeric(rownames(V12.122))

(plt <- Stratiplot(Depths ~ O.univ + G.ruber + G.tenel + G.pacR,
                   data = V12.122,  type = c("h","l","g","smooth")))

## Order taxa by WA in depth --- ephasises change over time
(plt <- Stratiplot(Depths ~ O.univ + G.ruber + G.tenel + G.pacR,
                   data = V12.122, type = c("h"), sort = "wa"))

## Using the default interface
spp.want <- c("O.univ","G.ruber","G.tenel","G.pacR")
(plt <- Stratiplot(V12.122[, spp.want], y = Depths,
                   type = c("poly", "g")))

## Adding zones to a Stratigraphic plot
## Default labelling and draw zone legend
## Here we choose 4 arbitrary Depths as the zone boundaries
set.seed(123)
Zones <-sample(Depths, 4)
Stratiplot(Depths ~ O.univ + G.ruber + G.tenel + G.pacR,
           data = V12.122, type = c("poly","g"),
           zones = Zones)

## As before, but supplying your own zone labels
zone.labs <- c("A","B","C","D","E")
Stratiplot(Depths ~ O.univ + G.ruber + G.tenel + G.pacR,
           data = V12.122, type = c("poly","g"),
           zones = Zones, zoneNames = zone.labs)

## Suppress the drawing of the zone legend
Stratiplot(Depths ~ O.univ + G.ruber + G.tenel + G.pacR,
           data = V12.122, type = c("poly","g"),
           zones = Zones, drawLegend = FALSE)

## Add zones and draw a legend, but do not label the zones
Stratiplot(Depths ~ O.univ + G.ruber + G.tenel + G.pacR,
           data = V12.122, type = c("poly","g"),
           zones = Zones, zoneNames = "")

## Show illustration of NA handling
set.seed(42)
dat <- data.frame(Depth = 1:20, LOI = runif(20), TC = NA)
dat <- within(dat, TC[sample(20, 10)] <- runif(10))
## default is 'na.action = "na.pass"'
Stratiplot(Depth ~ LOI + TC, data = dat, type = c("l","p"))
## to remove rows with NA, use 'na.action = "na.omit"'
Stratiplot(Depth ~ LOI + TC, data = dat, type = c("l","p"),
           na.action = "na.omit")

## Example of two proxies measured on different levels of core
## (Here measurements on alternate levels)
set.seed(5)
dat2a <- data.frame(Depth = seq(1, by = 2, length = 20), LOI = runif(20))
dat2b <- data.frame(Depth = seq(0, by = 2, length = 20), TC = runif(20))
dat2 <- join(dat2a, dat2b, na.replace = FALSE, split = FALSE)
dat2 <- dat2[order(dat2$Depth), ]
head(dat2)

## Default is to allow NA through formula, but drop them when plotting
Stratiplot(Depth ~ LOI + TC, data = dat2, type = c("l","p"))

## compare with this if we didn't suppress NA in default Stratiplot
## method (can't use formula interface for this yet
Stratiplot(dat2[,-1], dat2[,1], type = c("l","p"),
           na.action = "na.pass")
## Notice no lines are draw as there a no "sections" ithout missing
## levels. If you want/desire this behaviour then you can't use formula
## interface yet as there is no way to specify the na.action separately

## Works with matrices
M <- as.matrix(V12.122)
Stratiplot(M, Depths, type = c("h"))

## Custom variable labels using expressions
df <- data.frame(Age = 1:10, Var1 = rnorm(10), Var2 = rnorm(10),
                 Var3 = rnorm(10))
## Use a vector of expressions to label variables on plot
## See ?plotmath for syntax of expressions
exprs <- expression(delta^{15}*N,       # label for Var1
                    delta^{18}*O,       # label for Var2
                    delta^{13}*C)       # label for Var3
Stratiplot(Age ~ ., data = df, labelValues = exprs, varTypes = "absolute")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
