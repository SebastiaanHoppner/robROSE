\name{robROSE}

\alias{robROSE}

\title{Robust over-sampling of minority cases}

\description{Over-sampling minority cases while ignoring anomalies and incorporating the covariance structure of the minority cases.}

\usage{robROSE(formula, data, r, dup_size, alpha = 0.5, const = 1, seed = NULL)}

\arguments{
\item{formula}{an object of class "\code{\link{formula}}". The left-hand-side (response) should be a vector specifying the class labels.
The right-hand-side should be a series of vectors with the predictors. Transformations of variables are not allowed}

\item{data}{a data frame containing the variables.}

\item{r}{the desired fraction of minority cases. Either parameter \code{r} or \code{dup_size} must be specfied, not both.}

\item{dup_size}{this integer answers the question how many times \code{robROSE} should loop through the existing, real minority instances. So \code{dup_size} is the number representing the desired times of synthetic minority instances over the original number of minority instances. Either parameter \code{r} or \code{dup_size} must be specfied, not both.}

\item{alpha}{numeric parameter used by the \code{covMcd} function for controlling the size of the subsets over which the determinant is minimized (default is \code{0.5}).}

\item{const}{tuning constant that changes the volume of the elipsoids (default is 1).}

\item{seed}{optional positive integer. A single value, interpreted as an integer, recommended to specify seeds and keep trace of the generated sample.}
}

\details{robROSE over-samples only "clean" minority cases and does not include "outlying" (i.e. outcast) minority cases.}

\value{
\code{robROSE} returns an object of class "\code{robROSE}" which is a list containing the following components:
\item{call}{the matched call.}
\item{hmult}{internal values used.}
\item{ind.mino.out}{indices of orignal minority cases that are considered as outliers by robROSE.}
\item{tab.id.mino.new}{table that contains the number of synthetic cases created for each "clean" minority case.}
\item{mahdists}{robust Mahalanobis distances based on Minimum Covariance Discriminant (MCD) estimator (with default values \code{alpha = 0.5} and \code{nsamp = 500}) used to find outlying minority cases.}
\item{data}{data frame with synthetic minority cases.}
}

\references{Baesens, B., Hoppner, S., Ortner, I., and Verdonck T. (2020). robROSE: A robust approach for dealing with imbalanced data in fraud detection. \emph{arXiv:2003.11915}}

\author{Bart Baesens, Sebastiaan Hoppner, Irene Ortner, and Tim Verdonck}

\examples{
library(robROSE)
data(toy)

plot(toy$X1, toy$X2,
     pch = ifelse(toy$class == 0, 1, 18),
     cex = ifelse(toy$class == 0, 1, 1.8),
     col = ifelse(toy$class == 0, "dodgerblue", "black"),
     xlim = c(-3, 2), ylim = c(-4, 3))

result <- robROSE(class ~ ., data = toy, r = 0.5, const = 0.2, seed = 2020)

points(result$data$X1, result$data$X2, pch = 16, col = "red")
}
