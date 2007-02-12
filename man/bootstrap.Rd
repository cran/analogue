\name{bootstrap}
\alias{bootstrap}
\alias{bootstrap.default}
\alias{bootstrap.mat}}
\alias{print.bootstrap}
%\alias{print.bootstrap.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap estimation and errors }
\description{
  Function to calculate bootstrap statistics for transfer function
  models such as bootstrap estimates, model RMSEP, sample specific
  errors for predictions and summary statistics such as bias and
  \eqn{R^2} between oberved and estimated environment.
}
\usage{

bootstrap(object, ...)

\method{bootstrap}{default}(object, ...)

\method{bootstrap}{mat}(object, newdata, newenv, k,
          weighted = FALSE, n.boot = 1000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object for which bootstrap statistics are to be
    generated. Only objects of class \code{"mat"} currently supported.}
  \item{newdata}{a data frame containing samples for which bootstrap
    predictions and sample specific errors are to be generated. May be
    missing --- See Details. \code{"newdata"} must have the same number
    of columns as the training set data.}
  \item{newenv}{a vector containing environmental data for samples
    in \code{"newdata"}. Used to calculate full suite of errors for new
    data such as a test set with known environmental values. May be
    missing --- See Details. \code{"newenv"} must have the same number
    of rows as \code{"newdata"}.}
  \item{k}{numeric; how many modern analogues to use to generate the
    bootstrap statistics and, if requested, the predictions.}
  \item{weighted}{logical; should the weighted mean of the environment
    for the \code{"k"} modern analogues be used instead of the mean?}
  \item{n.boot}{Number of bootstrap samples to take.}
  \item{\dots}{arguments passed to other methods.}
}
\details{
  \code{bootstrap} is a fairly flexible function, and can be called with
  or without arguments \code{newdata} and \code{newenv}.

  If called with only \code{object} specified, then bootstrap estimates
  for the training set data are returned. In this case, the returned
  object will not include component \code{predictions}.

  If called with both \code{object} and \code{newdata}, then in addition
  to the above, bootstrap estimates for the new samples are also
  calculated and returned. In this case, component \code{predictions}
  will contain the apparent and bootstrap derived predictions and
  sample-specific errors for the new samples.

  If called with \code{object}, \code{newdata} and \code{newenv}, then
  the full \code{bootstrap} object is returned (as described in the
  Value section below). With environmental data now available for the
  new samples, residuals, RMSE(P) and \eqn{R^2} and bias statistics can
  be calculated.

  The individual components of \code{predictions} are the same as those
  described in the components relating to the training set data. For
  example, \code{returned.object$predictions$bootstrap} contains the
  components as \code{returned.object$bootstrap}.

  It is not usual for environmental data to be available for the new
  samples for which predictions are required. In normal
  palaeolimnological studies, it is more likely that \code{newenv} will
  not be available as we are dealing with sediment core samples from the
  past for which environmental data are not available. However, if
  sufficient training set samples are available to justify producing a
  training and a test set, then \code{newenv} will be available, and
  \code{bootstrap} can accomodate this extra information and calculate
  apparent and bootstrap estimates for the test set, allowing an
  independent assessment of the RMSEP of the model to be performed.
}
\value{
  A large object is returned with some or all of the following depending
  on whether \code{newdata} and \code{newenv} are supplied or not.

  \item{observed}{vector of observed environmental values.}
  \item{apparent}{a list containing the apparent or non-bootstrapped
    estimates for the training set. With the following components:
    \item{\code{estimated}}{estimated values for \code{"y"}, the
      environment.}
    \item{\code{residuals}}{model residuals.}
    \item{\code{r.squared}}{Apparent \eqn{R^2} between observed and
      estimated values  of \code{"y"}.}
    \item{\code{avg.bias}}{Average bias of the model residuals.}
    \item{\code{max.bias}}{Maximum bias of the model residuals.}
    \item{\code{rmse}}{Apparent error (RMSE) for the model.}
    \item{\code{k}}{numeric; indicating the size of model used in
      estimates and predictions.}
  }
  \item{bootstrap}{a list containing the bootstrap estimates for the
    training set. With the following components:
    \item{\code{estimated}}{Bootstrap estimates for \code{"y"}.}
    \item{\code{residuals}}{Bootstrap residuals for \code{"y"}.}
    \item{\code{r.squared}}{Bootstrap derived \eqn{R^2} between observed
      and estimated values of \code{"y"}.}
    \item{\code{avg.bias}}{Average bias of the bootstrap derived model
      residuals.}
    \item{\code{max.bias}}{Maximum bias of the bootstrap derived model
      residuals.}
    \item{\code{rmsep}}{Bootstrap derived RMSEP for the model.}
    \item{\code{s1}}{Bootstrap derived S1 error component for the
      model.}
    \item{\code{s2}}{Bootstrap derived S2 error component for the
      model.} 
    \item{\code{k}}{numeric; indicating the size of model used in
      estimates and predictions.}
  }
  \item{sample.errors}{a list containing the bootstrap-derived sample
    specific errors for the training set. With the following components:
    \item{\code{rmsep}}{Bootstrap derived RMSEP for the training set
      samples.} 
    \item{\code{s1}}{Bootstrap derived S1 error component for training
      set samples.}
    \item{\code{s2}}{Bootstrap derived S2 error component for training
      set samples.}
  }
  \item{weighted}{logical; whether the weighted mean was used instead of
    the mean of the environment for \emph{k}-closest analogues.}
  \item{auto}{logical; whether \code{"k"} was choosen automatically or
    user-selected.}
  \item{n.boot}{numeric; the number of bootstrap samples taken.}
  \item{predictions}{a list containing the apparent and
    bootstrap-derived estimates for the new data, with the following
    components:
    \item{\code{observed}}{the observed values for the new samples ---
      only if \code{newenv} is provided.}
    \item{\code{apparent}}{a list containing the apparent or
      non-bootstrapped estimates for the new samples. A list with the
      same components as \code{apparent}, above.
    }
    \item{\code{bootstrap}}{a list containing the bootstrap estimates
      for the new samples, with some or all of the same components as
      \code{bootstrap}, above.}
    \item{\code{sample.errors}}{a list containing the bootstrap-derived
      sample specific errors for the new samples, with some or all of
      the same components as \code{sample.errors}, above.}
  }
}
\references{

  Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C. and ter Braak,
  C.J.F. (1990). Diatoms and pH reconstruction. \emph{Philosophical
    Transactions of the Royal Society of London; Series B}, \strong{327};
  263--278.

}
\author{Gavin L. Simpson }
\seealso{\code{mat}, \code{plot.mat} \code{summary.bootstrap}}
\examples{
## continue the RLGH and SWAP example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## bootstrap training set
swap.boot <- bootstrap(swap.mat, n.boot = 100)
swap.boot
summary(swap.boot)

## bootstrap with predictions:
rlgh.boot <- bootstrap(swap.mat, rlgh, n.boot = 100)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
