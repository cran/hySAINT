#' Creating initial parents\cr
#'
#' This function gives initial parents.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects.
#' @param y Response variable. A \code{n}-dimensional vector.
#' @param EVAoutput The output from function \code{EVA}
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param r1 At most how many main effects do you want to include in your model?.
#'        For high-dimensional data, \code{r1} cannot be larger than the number of
#'        screened main effects.
#' @param r2 At most how many interaction effects do you want to include in your model?
#' @param numElite Number of elite parents. Default is 40.
#'
#' @return Initial parents. A numeric matrix with dimensions \code{numElite} by \code{r1+r2}.
#' @export
#'
#' @seealso \code{\link{EVA}}
#'
#' @examples
#' set.seed(0)
#' interaction.ind <- t(combn(10,2))
#' X <- matrix(rnorm(100*10,1,0.1), 100, 10)
#' epl <- rnorm(100,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,3]+X[,1]*X[,2]+X[,1]*X[,3]+epl
#' EVAoutput <- EVA(X, y, r1 = 5, sigma = 0.01,
#'   interaction.ind = interaction.ind)
#' myParent <- Initial(X = X, y = y, EVAoutput, r1 = 5, r2 = 2)

Initial <- function(X, y, EVAoutput, heredity = "Strong", r1, r2, numElite = 40){

  # DEFINE INITIAL PARENT MATRIX
  initial_parents <- matrix(0, nrow = numElite, ncol = (r1+r2))

  # BEGIN FILL R1 + R2 IDX. THE HEREDITY CONDITION IS AUTOMATICALLY SATISFIED HERE
  for (i in 1:mydim(initial_parents)[1]) {
    filler.initial <- sampleGpool(X, EVAoutput, heredity, r1, r2)
    initial_parents[i,1:(length(filler.initial))] <- c(filler.initial)
  }
  return(initial_parents)
}
