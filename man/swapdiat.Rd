\name{swapdiat}
\alias{swapdiat}
\docType{data}
\title{SWAP sub-fossil diatom and pH training set}
\description{
  The Surface Waters Acidifcation Project (SWAP) Palaeolimnology
  Programme diatom surface sample training set.
}
\usage{data(swapdiat)}
\format{
  A data frame of observations on 111 diatom taxa for a 74 lake subset
  of the SWAP diatom training set.
}
\details{
  This data set contains a 74-lake subset of the original SWAP diatom
  training set.

  The variable names (\code{\link{colnames}}) are DIATCODE codes for
  individual taxa.
}
\references{

  Stevenson, A.C., Juggins, S., Birks, H.J.B., Anderson, D.S., Anderson,
  N.J., Battarbee, R.W., Berge, F., Davis, R.B., Flower, R.J., Haworth,
  E.Y., Jones, V.J., Kingston, J.C., Kreiser, A.M., Line, J.M., Munro,
  M.A.R., and Renberg, I. (1995). \emph{The Surface Waters Acidification
  Project Palaeolimnology programme: modern diatom/lake-water chemistry
  data-set}. ENSIS Publishing, London.

}
\seealso{\code{\link{swappH}}}
\examples{
data(swapdiat)
str(swapdiat)
}
\keyword{datasets}
