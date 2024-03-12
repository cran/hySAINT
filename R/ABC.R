#' ABC Evaluation\cr
#'
#' Gives ABC score for each fitted model. For a model I, the ABC is defined as
#' \deqn{ABC(I)=\sum\limits_{i=1}^n\bigg(Y_i-\hat{Y}_i^{I}\bigg)^2+2r_I\sigma^2+\lambda\sigma^2C_I.}
#' When comparing ABC of fitted models to the same dataset, the smaller the ABC, the better fit.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects.
#' @param y Response variable. A \code{n}-dimensional vector.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param sigma The standard deviation of the noise term. In practice, sigma is usually
#'        unknown. Users can estimate sigma from function \code{selectiveInference::estimateSigma},
#'        then use the output as the sigma value. See examples for details.
#' @param varind A numeric vector that specifies the indices of variables to be extracted from \code{X}.
#'        Default is "No".
#' @param interaction.ind A two-column numeric matrix. Each row represents a unique
#'        interaction pair, with the columns indicating the index numbers of the variables
#'        involved in each interaction. Note that interaction.ind must be generated
#'        outside of this function using \code{t(utils::combn(p,2))}. See Example section for
#'        details.
#' @param lambda A numeric value defined by users. The number needs to satisfy the condition:
#'        \eqn{\lambda\geq 5.1/log(2)}. Default is 10.
#'
#' @return A numeric value is returned. It represents the ABC score of the fitted model.
#' @export
#'
#' @references
#' Ye, C. and Yang, Y., 2019. \emph{High-dimensional adaptive minimax sparse estimation with interactions.}
#'
#' @examples # When sigma is known
#' set.seed(0)
#' interaction.ind <- t(combn(4,2))
#' X <- matrix(rnorm(50*4,1,0.1), 50, 4)
#' epl <- rnorm(50,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,1]*X[,2] + epl
#'  ABC(X, y, sigma = 0.01, varind = c(1,2,5), interaction.ind = interaction.ind)
#'
#' @examples # When sigma is not known
#' full <- Extract(X, varind = c(1:(dim(X)[2]+dim(interaction.ind)[1])), interaction.ind)
#' sigma <- selectiveInference::estimateSigma(full, y)$sigmahat # Estimate sigma

ABC <- function(X, y, heredity = "Strong", sigma, varind = NULL, interaction.ind = NULL, lambda = 10){
  if (is.null(interaction.ind)){
    stop("Interaction.ind is missing. Use t(utils::combn()) to generate interaction matrix.")
  }
  if (as.logical(any(duplicated(varind[which(varind!=0)])))){
    stop("There cannot be duplicated values in varind.")
  }

  p <- dim(X)[2]
  pi0 <- 1e-10
  pi1 <- pi2 <- pi3 <- (1-pi0)/3
  n <- dim(X)[1]
  sigma <- sigma
  data_extract <- Extract(X, varind, interaction.ind)
  data <- data_extract
  r.I <- Matrix::rankMatrix(data)[1]

  # else{
  #   r.I <- Matrix::rankMatrix(X)[1]
  #   data <- X
  #   # colnames(data) <- paste0("X", c(1:p))
  # }
  k2 <- sum(as.numeric(gsub(".*?([0-9]+).*", "\\1",  colnames(data))) > p)
  k1 <- sum(as.numeric(gsub(".*?([0-9]+).*", "\\1",  colnames(data))) <= p)
  if (r.I < ncol(data)){
    yhat <- pracma::orth(data)%*%solve(crossprod(pracma::orth(data)),
                                       t(pracma::orth(data))%*%y, tol = 1e-50)
  }else{
    yhat <- (data)%*%solve(crossprod(data), t(data)%*%y, tol = 1e-50)
  }
  SSE <- sum((y - yhat)^2)
  if (heredity == "Strong"){
    C.I.strong <- -log(pi1)+log(min(p,n))+log(min(mychoose(k1),n)) + log(choose(p,k1)) + log(choose(choose(k1,2),k2))
    ABC <- SSE+2*r.I*sigma^2+lambda*sigma^2*C.I.strong
  }
  else if (heredity == "Weak"){
    K <- k1*p-choose(k1,2)-k1
    C.I.weak <- -log(pi2)+log(min(p,n))+log(min(K,n)) + log(choose(p,k1))+ log(choose(K,k2))
    ABC <- SSE+2*r.I*sigma^2+lambda*sigma^2*C.I.weak
  }
  else if (heredity == "No"){
    C.I.no <- -log(pi3)+log(min(p,n))+log(min(choose(p,2),n)) + log(choose(p,k1)) + log(choose(choose(p,2),k2))
    ABC <- SSE+2*r.I*sigma^2+lambda*sigma^2*C.I.no
  }
  return(ABC)
}
