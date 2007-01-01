\name{mat}
\alias{mat}
\alias{mat.default}
\alias{print.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Modern Analogue Technique transfer function models}
\description{
  Modern Analogue Technique (MAT) transfer function models for
  palaeoecology. The fitted values are the, possibly weighted, averages
  of the environment for the \emph{k}-closest modern analogues. MAT is a
  \emph{k}-NN method.
}
\usage{
mat(x, \dots)

\method{mat}{default}(x, y,
    method = c("euclidean", "SQeuclidean", "chord", "SQchord",
               "bray", "chi.square", "SQchi.square",
               "information", "chi.distance", "manhattan",
               "kendall", "gower", "alt.gower", "mixed"),
    \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data frame containing the training set data, usually
    species data.}
  \item{y}{a vector containing the response variable, usually
    environmental data to be predicted from \code{x}.}
  \item{method}{a character string indicating the dissimilarity
    (distance) coefficient to be used to define modern analogues. See
    Details, below.}
  \item{\dots}{arguments to or from other methods.}
}
\details{
  The Modern Analogue Technique (MAT) is perhaps the simplest of the
  transfer function models used in palaeoecology. An estimate of the
  environment, \eqn{x}, for the response for a fossil sample, \eqn{y},
  is the, possibly weighted, mean of that variable across the
  \emph{k}-closest modern analogues selected from a modern  training set
  of samples. If used, weights are the reciprocal of the dissimilarity
  between the fossil sample and each modern analogue.

  Pairwise sample dissimilarity is defined by dissimilarity or
  distance coefficients. A variety of coefficients are supported --- see
  \code{\link{distance}} for details of the supported coefficients.

  \emph{k} is chosen by the user. The simplest choice for \emph{k} is to
  evaluate the RMSE of the difference between the predicted and observed
  values of the environmental variable of interest for the training set
  samples for a sequence of models with increasing \emph{k}. The number
  of analogues chosen is the value of \emph{k} that has lowest
  RMSE. However, it should be noted that this value is biased as the
  data used to build the model are also used to test the predictive
  power.

  An alternative approach is to employ an optimisation data set on which
  to evaluate the size of \eqn{k} that provides the lowest RMSEP. This
  may be impractical with smaller sample sizes.

  A third option is to bootstrap re-sample the training set many times. At
  each bootstrap sample, predictions for samples in the bootstrap test
  set can be made for \eqn{k = 1, ..., n}, where \eqn{n} is the
  number of samples in the training set. \eqn{k} can be chosen from the
  model with the lowest RMSEP. See function \code{\link{bootstrap}} for
  further details on choosing \eqn{k}.

  The output from \code{\link{summary.mat}} can be used to choose
  \eqn{k} in the first case above. For predictions on an optimsation or
  test set see \code{\link{predict.mat}}. For bootstrap resampling of
  \code{\link{mat}} models, see \code{\link{bootstrap}.}
}
\value{
  Returns an object of class \code{mat} with the following components:
  \item{standard}{list; the model statistics based on simple
    averages of \emph{k}-closest analogues. See below.}
  \item{weighted}{list; the model statistics based on weighted  of
    \emph{k}-closest analogues. See below.}
  \item{Dij}{matrix of pairwise sample dissimilarities for the training
    set \code{x}.}
  \item{orig.x}{the original training set data.}
  \item{orig.y}{the original environmental data or response, \code{y}.}
  \item{call}{the matched function call.}
  \item{method}{the dissimilarity coefficient used.}
}
\note{
  Lists \code{"standard"} and \code{"weighted"} both contain the
  following elements:
  \itemize{
    \item{\code{est}}{a matrix of estimated values for the training set
      samples for models using \eqn{k} analogues, where \eqn{k = 1,
	..., n}. \eqn{n} is the number of smaples in the training
      set. Rows contain the values for each model of size \eqn{k}, with
      colums containing the estimates for each training set sample.}
    \item{\code{resid}}{matrix; as for \code{"est"}, but containing the
      model residuals.}
    \item{\code{rmse}}{vector; containing the apparent root mean square
      error.}
    \item{\code{avg.bias}}{vector; contains the average bias (mean of
      residuals) for models using \emph{k} analogues, where \eqn{k = 1,
	..., n}. \eqn{n} is the number of smaples in the training set.}
    \item{\code{max.bias}}{vector; as for \code{"avg.bias"}, but
      containing the maximum bias statistics.}
    \item{\code{r.squared}}{vector; as for \code{"avg.bias"}, but
      containing the \eqn{R^2} statistics.}
    }
}
\references{
  
  Gavin, D.G., Oswald, W.W., Wahl, E.R. and Williams, J.W. (2003) A
  statistical approach to evaluating distance metrics and analog
  assignments for pollen records. \emph{Quaternary Research}
  \strong{60}, 356--367.
  
  Overpeck, J.T., Webb III, T. and Prentice I.C. (1985) Quantitative
  interpretation of fossil pollen spectra: dissimilarity coefficients and
  the method of modern analogues. \emph{Quaternary Research} \strong{23},
  87--108.

  Prell, W.L. (1985) The stability of low-latitude sea-surface
  temperatures: an evaluation of the CLIMAP reconstruction with emphasis
  on the positive SST anomalies, Report TR 025. U.S. Department of
  Energy, Washington, D.C.
  
  Sawada, M., Viau, A.E., Vettoretti, G., Peltier, W.R. and Gajewski,
  K. (2004) Comparison of North-American pollen-based temperature and
  global lake-status with CCCma AGCM2 output at 6 ka. \emph{Quaternary
    Science Reviews} \strong{23}, 87--108.

}
\author{Gavin L. Simpson}
\seealso{\code{\link{summary.mat}}, \code{\link{bootstrap}} for boostrap
  resampling of MAT models, \code{\link{predict.mat}} for making
  predictions from MAT models, \code{\link{fitted.mat}} and
  \code{\link{resid.mat}} for extraction of fitted values and residuals
  from MAT models respectively. \code{\link{plot.mat}} provides a
  \code{\link{plot.lm}}-like plotting tool for MAT models.}
\examples{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")
swap.mat

## model summary
summary(swap.mat)

## fitted values
fitted(swap.mat)

## model residuals
resid(swap.mat)

## draw summary plots of the model
par(mfrow = c(2,2))
plot(swap.mat)
par(mfrow = c(1,1))

## reconstruct for the RLGH core data
rlgh.mat <- predict(swap.mat, rlgh, k = 10)
rlgh.mat
summary(rlgh.mat)
rlgh.Wmat <- predict(swap.mat, rlgh, k = 10, weighted  = TRUE)
rlgh.Wmat
summary(rlgh.Wmat)

## plot of pH change in the RLGH
depths <- as.numeric(colnames(rlgh.mat$predictions$apparent$predicted))
n.analogues <- rlgh.mat$predictions$apparent$k
plot(rlgh.mat$predictions$apparent$predicted[n.analogues, ], depths,
     ylim = rev(range(depths)),
     xlab = "pH",
     ylab = "Depth (cm)",
     main = "Estimated pH",
     type = "l")
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
