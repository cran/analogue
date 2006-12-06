\name{bayesF}
\alias{bayesF}
\alias{print.bayesF}
\concept{likelihood ratios}
\concept{Bayes factors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayes factors}
\description{
  Calculates Bayes factors or likelihood ratios of analogue and
  no-analogue results.
}
\usage{
bayesF(x, prior = NULL)
%bayesF(x, which = "both", prior = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \dQuote{\code{roc}}.}
  %\item{which}{character; one of \code{"both"}, \code{"pos"}, or
  %  \code{"neg"} describing the Bayes factors required. \code{"pos"}
  %  calculates the Bayes factor or likelihood ratio of a positive event
  %  (analogue) and \code{"neg"} of a negative event (no-analogue).}
  \item{prior}{numeric; the prior probabilities of analogue and
    no-analogue, provided as a vector of length 2 whose elements sum to
    1. If not provided, the function will use the relative occurences of
    analogue and no analogue situations used to evaluate the ROC curve.}
}
\details{
  LR(+), is the likelihood ratio of a positive test result, that
  the value of \emph{d} assigns the sample to the group it belongs
  to. LR(-) is the likelihood ratio of a negative test result, that the
  value of \emph{d} assigns the sample to the wrong group.

  LR(+) is defined as \eqn{LR(+) = TPF / FPF} (or sensitivity / (1 -
  specificity)), and LR(-) is defined as \eqn{LR(-) = FPF / TNF} (or (1
  - sensitivity) / specificity), in Henderson (1993).

  The posterior probability of analogue given a dissimilarity is the 
  LR(+) likelihood ratio values multiplied by the prior odds of
  analogue, for given values of the dissimilarity, and is then converted
  to a  probability.
}
\value{
  A list with the followin components, some of which may be \code{NULL}
  depending on argument \code{which}:
  \item{pos}{Bayes factor or likelihood ratio of a positive event
    (analogue).}
  \item{neg}{Bayes factor or likelihood ratio of anegative event
    (analogue).}
  \item{posterior}{list with components \code{pos} and \code{neg}
    containing the posterior probabilities of positive and negative
    events, respectively.}
  \item{prior}{list with components \code{pos} and \code{neg}
    containing the prior probabilities of positive and negative
    events, respectively.}
  \item{roc.points}{vector of points at which the ROC curve was
    evaluated and for which Bayes factors and prior and posterior
    probabilities are available.}
  \item{optimal}{numeric; the optimal dissimilarity, as assessed by the
    ROC curve.}
  %\item{method}{one of \code{"both"}, \code{"pos"}, or \code{"neg"}.}
  \item{object}{name of the object passed as argument \code{x}.}
}
\references{
  
  Brown, C.D., and Davis, H.T. (2006) Receiver operating characteristics
  curves and related decision measures: A tutorial. \emph{Chemometrics
    and Intelligent Laboratory Systems} \bold{80}, 24--38.
  
  Gavin, D.G., Oswald, W.W., Wahl, E.R. and Williams, J.W. (2003) A
  statistical approach to evaluating distance metrics and analog
  assignments for pollen records. \emph{Quaternary Research}
  \strong{60}, 356--367.

  Henderson, A.R. (1993) Assessing test accuracy and its clinical
  consequences: a primer for receiver operating characteristic curve
  analysis. \emph{Annals of Clinical Biochemistry} \strong{30},
  834--846.

}
\author{Gavin L. Simpson}
\seealso{\code{\link{roc}} and \code{\link{plot.bayesF}}.}
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
\keyword{univar}
\keyword{methods}
