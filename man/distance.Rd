\name{distance}
\alias{distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flexibly calculate dissimilarity or distance measures }
\description{
  Flexibly calculates distance or dissimilarity measures between a
  training set \code{"x"} and a fossil or test set \code{"y"}. If
  \code{"y"} is not supplied then the pairwise dissimilarities between
  samples in the training set, \code{"x"}, are calculated.
}
\usage{
distance(x, y, method = c("euclidean", "SQeuclidean",
         "chord", "SQchord", "bray", "chi.square",
         "SQchi.square", "information", "chi.distance",
         "manhattan", "kendall", "gower"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data frame or matrix containing the training set samples.}
  \item{y}{data frame or matrix containing the fossil or test set
    samples.}
  \item{method}{character; which choice of dissimilarity coefficient to
    use. One of the listed options. See Details below}
}
\details{
  A range of dissimilarity coefficients can be used to calculate
  dissimilarity between samples. The following are currently available:

  \tabular{ll}{
    \code{euclidean}
    \tab \eqn{d_{jk} = \sqrt{\sum_i (x_{ij}-x_{ik})^2}}{d[jk] = sqrt(sum
      (x[ij]-x[ik])^2)}
    \cr
    \code{SQeuclidean}
    \tab \eqn{d_{jk} = \sum_i (x_{ij}-x_{ik})^2}{d[jk] = sum
      (x[ij]-x[ik])^2}
    \cr
    \code{chord}
    \tab \eqn{d_{jk} = \sqrt{\sum_i
	(\sqrt{x_{ij}}-\sqrt{x_{ik}})^2}}{d[jk] = sqrt(sum((sqrt(x) -
      sqrt(y))^2))}
    \cr
    \code{SQchord}
    \tab \eqn{d_{jk} = \sum_i (\sqrt{x_{ij}}-\sqrt{x_{ik}})^2}{d[jk] =
      sum((sqrt(x) - sqrt(y))^2)}
    \cr
    \code{bray}
    \tab \eqn{d_{jk} = \frac{\sum_i |x_{ij} - x_{ik}|}{\sum_i (x_{ij} +
	x_{ik})}}{d[jk] = sum(abs(x - y)) / sum(x + y)}
    \cr
    \code{chi.square}
    \tab \eqn{d_{jk} = \sqrt{\sum_i \frac{(x_{ij} - x_{ik})^2}{x_{ij} +
	x_{ik}}}}{d[jk] = sqrt(sum(((x - y)^2) / (x + y)))}
    \cr
    \code{SQchi.square}
    \tab \eqn{d_{jk} = \sum_i \frac{(x_{ij} - x_{ik})^2}{x_{ij} +
	x_{ik}}}{d[jk] = sum(((x - y)^2) / (x + y))}
    \cr
    \code{information}
    \tab \eqn{d_{jk} = \sum_i (p_{ij}log(\frac{2p_{ij}}{p_{ij} + p_{ik}})
      + p_{ik}log(\frac{2p_{ik}}{p_{ij} + p_{ik}}))}{d[jk] = sum((x[ij] *
      log((2 * x[ij]) / (x[ij] + x[ik]))) + (x[ik] * log((2 * x[ik]) /
      (x[ij] + x[ik]))))}
    \cr
    \code{chi.distance}
    \tab \eqn{d_{jk} = \sqrt{\sum_i (x_{ij}-x_{ik})^2 / (x_{i+} /
	x_{++})}}{d[jk] = sqrt(sum((x[ij] - x[ik])^2 / (x[i+] / x[++])))}
    \cr
    \code{manhattan}
    \tab \eqn{d_{jk} = \sum_i (|x_{ij}-x_{ik}|)}{d[jk] = sum
      (|x[ij]-x[ik]|)}
    \cr
    \code{kendall}
    \tab \eqn{d_{jk} = \sum_i MAX_i - minimum(x_{ij}, x_{ik})}{d[jk] = sum
      (MAX[i] - min(x[ij]-x[ik]))}
    \cr
    \code{gower}
    \tab \eqn{d_{jk} = \sqrt{2\sum_i\frac{|p_{ij} -
	  p_{ik}|}{R_i}}}{sqrt(2 * sum(abs(x[ij] - x[ik]) / R[i]))}
    \cr
    \tab where \eqn{R_i}{R[i]} is the range of proportions for pollen type
    \emph{i}
  }
}
\value{
  A matrix of dissimilarities where columns are the samples in
  \code{"y"} and the rows the samples in \code{"x"}. If \code{"y"} is
  not provided then a square, symmetric matrix of pairwise sample
  dissimilarities for the training set \code{"x"} is returned.
}
\note{
  The dissimilarities are calculated in native R code. As such, other
  implementations (see See Also below) will be quicker. This is done for
  one main reason - it is hoped to allow a user defined function to be
  supplied as argument \code{"method"} to allow for user-extension of
  the available coefficients.

  The other advantage of \code{distance} over other implementations, is
  the simplicity of calculating only the required pairwise sample
  dissimilarities between each fossil sample (\code{"y"}) and each
  training set sample (\code{"x"}). To do this in other implementations,
  you would need to merge the two sets of samples, calculate the full
  dissimilarity matrix and then subset it to achieve similar results.
}
\references{

  Faith, D.P., Minchin, P.R. and Belbin, L. (1987) Compositional
  dissimilarity as a robust measure of ecological
  distance. \emph{Vegetatio} \strong{69}, 57--68.
  
  Gavin, D.G., Oswald, W.W., Wahl, E.R. and Williams, J.W. (2003) A
  statistical approach to evaluating distance metrics and analog
  assignments for pollen records. \emph{Quaternary Research}
  \strong{60}, 356--367.

  Kendall, D.G. (1970) A mathematical approach to
  seriation. \emph{Philosophical Transactions of the Royal Society of
    London - Series B} \strong{269}, 125--135.
  
  Overpeck, J.T., Webb III, T. and Prentice I.C. (1985) Quantitative
  interpretation of fossil pollen spectra: dissimilarity coefficients and
  the method of modern analogues. \emph{Quaternary Research} \strong{23},
  87--108. 
  
  Prentice, I.C. (1980) Multidimensional scaling as a research tool in
  Quaternary palynology: a review of theory and methods. \emph{Review of
    Palaeobiology and Palynology} \strong{31}, 71--104.
 
}
\author{Gavin L. Simpson }
\seealso{\code{\link[vegan]{vegdist}} in package \code{vegan},
  \code{\link[cluster]{daisy}} in package \code{cluster}, and
  \code{\link[stats]{dist}} provide comparable functionality for the
  case of missing \code{"y"} and are implemented in compiled code, so
  will be faster.} 
\examples{
## simple example using dummy data
train <- data.frame(matrix(runif(200), ncol = 10))
rownames(train) <- LETTERS[1:20]
colnames(train) <- as.character(1:10)
fossil <- data.frame(matrix(runif(100), ncol = 10))
colnames(fossil) <- as.character(1:10)
rownames(fossil) <- letters[1:10]

## calculate distances/dissimilarities between train and fossil
## samples
test <- distance(train, fossil)

## calculate pairwise distances/dissimilarities for training
## set samples
test2 <- distance(train)

## using a different coefficient, chi-square distance
test <- distance(train, fossil, method = "chi.distance")
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
