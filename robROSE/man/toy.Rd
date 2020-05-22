\name{toy}
\docType{data}
\alias{toy}
\title{Toy data}
\description{
The small data set that is used in the paper to illustrate robROSE.
}
\usage{data(toy)}
\format{
A data frame containing 997 observations and 3 variables.
Variables \code{X1} and \code{X2} are numeric, and factor variable \code{class} is the binary response target taking values 0 and 1.
}
\keyword{datasets}
\examples{
data(toy)
str(toy)
head(toy)
summary(toy)
}
