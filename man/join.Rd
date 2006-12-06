\name{join}
\alias{join}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Merge species data sets on common columns (species)}
\description{
  Merges any number of species matrices on their common columns to
  create a new data set with number of columns equal to the number of
  unqiue columns across all data frames. Needed for analysis of fossil
  data sets with respect to training set samples. 
}
\usage{
join(\dots, verbose = FALSE, na.replace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{data frames containing the data sets to be merged}
  \item{verbose}{logical; if \code{TRUE}, the function prints out the
    dimensions of the data frames in \code{"\dots"}, as well as those of
    the returned, merged data frame.}
  \item{na.replace}{logical; samples where a column in one data frame
    that have no matching column in the other will contain missing
    values (\code{NA}). If \code{na.replace} is \code{TRUE}, these
    missing values are replaced with zeros. This is standard practice in
    ecology and palaeoecology. If you want to replace with another
    value, then set \code{na.replace} to \code{FALSE} and do the
    replacement later.}
}
\value{
  A data frame containing the merged data sets.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{merge}}}
\examples{
## load the example data
data(swapdiat)
data(swappH)
data(rlgh)
## process so common set of columns for training and test
## number of training set samples
n.train <- nrow(swapdiat)

## merge training and test set on columns
dat <- join(swapdiat, rlgh, verbose = TRUE)

## convert to proportions
dat <- dat / 100
## subset data back into training and test sets
swapdiat <- dat[1:n.train, ]
rlgh <- dat[(n.train+1):nrow(dat), ]
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
