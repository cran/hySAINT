% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Crossover.R
\name{Crossover}
\alias{Crossover}
\title{Performing crossover\cr}
\usage{
Crossover(X, myParent, EVAoutput, heredity = "Strong", r1, r2, numElite = 40)
}
\arguments{
\item{X}{Input data. An optional data frame, or numeric matrix of dimension
\code{n} observations by \code{p} main effects.}

\item{myParent}{A numeric matrix with dimension \code{numElite} by \code{r1 + r2}.}

\item{EVAoutput}{The output from function \code{EVA}.}

\item{heredity}{Whether to enforce Strong, Weak, or No heredity. Default is "Strong".}

\item{r1}{At most how many main effects do you want to include in your model?.
For high-dimensional data, \code{r1} cannot be larger than the number of screened main effects.}

\item{r2}{At most how many interaction effects do you want to include in your model?}

\item{numElite}{Number of elite parents. Default is 40.}
}
\value{
Offspring. If crossover occurred, it returns a numeric matrix with dimensions
\code{choose(numElite,2)} by \code{r1+r2}. Otherwise, \code{numElite} by \code{r1 + r2}.
}
\description{
This function gives offspring from parents. It performs crossover at a fixed
probability of 0.6.
}
\examples{
set.seed(0)
interaction.ind <- t(combn(10,2))
X <- matrix(rnorm(100*10,1,0.1), 100, 10)
epl <- rnorm(100,0,0.01)
y <- 1+X[,1]+X[,2]+X[,3]+X[,1]*X[,2]+X[,1]*X[,3]+epl
EVAoutput <- EVA(X, y, r1 = 5, sigma = 0.01,
  interaction.ind = interaction.ind)
myParent <- Initial(X = X, y = y, EVAoutput, r1 = 5, r2 = 2)
Offsprings <- Crossover(X, myParent, EVAoutput, r1 = 5, r2 = 2)
}
\seealso{
\code{\link{EVA}}, \code{\link{Initial}}.
}
