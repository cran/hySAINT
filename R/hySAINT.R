#' Hybrid Genetic and Simulated Annealing Algorithm\cr
#'
#' This is the main function of package hySAINT. It implements both genetic
#' algorithm and simulated annealing. The simulated annealing technique is
#' used within mutation operator.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects.
#' @param y Response variable. A \code{n}-dimensional vector.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param r1 At most how many main effects do you want to include in your model?.
#'        For high-dimensional data, \code{r1} cannot be larger than the number of screened main effects.
#' @param r2 At most how many interaction effects do you want to include in your model?
#' @param sigma The standard deviation of the noise term. In practice, sigma is usually
#'        unknown. Users can estimate sigma from function \code{selectiveInference::estimateSigma},
#'        then use the output as the sigma value.
#' @param interaction.ind A two-column numeric matrix. Each row represents a unique
#'        interaction pair, with the columns indicating the index numbers of the variables
#'        involved in each interaction. Note that interaction.ind must be generated
#'        outside of this function using \code{t(utils::combn(p,2))}. See Example section for
#'        details.
#' @param varind A numeric vector that specifies the indices of variables to be extracted from \code{X}.
#' @param numElite Number of elite parents. Default is 40.
#' @param max.iter Maximum number of iterations. Default is 500.
#' @param initial.temp Initial temperature. Default is 1000.
#' @param cooling.rate A numeric value represents the speed at which the
#'        temperature decreases. Default is 0.95.
#' @param lambda A numeric value defined by users. The number needs to satisfy the condition:
#'        \eqn{\lambda\geq 5.1/log(2)}. Default is 10.
#'
#' @return An object with S3 class \code{"hySAINT"}.
#' \item{Final.variable.names}{Name of the selected effects.}
#' \item{Final.variable.idx}{Index of the selected effects.}
#' \item{Final.model.score}{Final Model ABC.}
#' \item{All.iter.score}{Best ABC scores from initial parents and all iterations. }
#' @export
#'
#' @seealso \code{\link{ABC}}, \code{\link{EVA}}, \code{\link{Initial}},
#'          \code{\link{Crossover}}, \code{\link{Mutation}}
#'
#' @importFrom utils combn
#' @importFrom selectiveInference estimateSigma
#' @importFrom Matrix rankMatrix
#' @importFrom pracma orth
#' @importFrom VariableScreening screenIID
#' @importFrom SIS SIS
#' @importFrom stats rbinom
#' @importFrom stats runif
#'
#' @examples
#' set.seed(0)
#' interaction.ind <- t(combn(10,2))
#' X <- matrix(rnorm(100*10,1,0.1), 100, 10)
#' epl <- rnorm(100,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,3]+X[,1]*X[,2]+X[,1]*X[,3]+epl
#' hySAINT(X, y, r1 = 5, r2 = 2, sigma = 0.01, interaction.ind = interaction.ind, max.iter = 5)

hySAINT <- function(X, y, heredity = "Strong", r1, r2, sigma, interaction.ind = NULL, varind = NULL, numElite = 40, max.iter = 500, initial.temp = 1000, cooling.rate = 0.95, lambda = 10){
  if (is.null(interaction.ind)) {
    stop("Interaction.ind is missing. Use t(utils::combn()) to generate interaction matrix.")
  }
  p <- dim(X)[2]

  # RUN EVA TO GET RANKED INITIAL INTERACTION POOL. AND RANKED MAIN EFFECTS.
  EVAoutput <- EVA(X, y, heredity, r1, sigma, varind , interaction.ind, lambda)

  # RUN INITIAL TO GET INITIAL POPULATION
  myParent <- Initial(X, y, EVAoutput, heredity, r1, r2, numElite)

  # CALCULATE INITIAL GENERATION ABC
  G.initial.matrix <- myfamily.eva(myParent, r1, r2, X , y, heredity, sigma, varind, interaction.ind, lambda)

  # OBTAIN THE BEST ABC SCORE OF THE INITIAL GENERATION. AND CALL IT G BEST
  G.best <- as.numeric(G.initial.matrix[1,((r1+r2)+1)])
  EVAoutput <- EVAoutput
  iter <- 0
  G.current <- list()

  # MAIN WHILE LOOP BEGINS HERE
  while (iter < max.iter) {
    iter <- iter + 1
    Offsprings <- Crossover(X, myParent, EVAoutput, heredity, r1, r2, numElite)
    Mutant <- Mutation(myParent, EVAoutput, r1, r2, initial.temp, cooling.rate, X, y, heredity, sigma,varind, interaction.ind, lambda)

    Myfamily <- rbind(myParent, Offsprings, Mutant)

    # CALCULATE CURRENT GENERATION ABC
    G.current.matrix <- myfamily.eva(Myfamily[,1:(r1+r2)], r1, r2, X , y, heredity, sigma, varind, interaction.ind, lambda)

    # OBTAIN THE BEST ABC SCORE OF THE CURRENT GENERATION. AND CALL IT G CURRENT
    G.current[[iter]] <- as.numeric(G.current.matrix[1, ((r1+r2)+1)])
    G.current.numElite <- remove.dup.rows(G.current.matrix, r1, r2)
    myParent <- G.current.numElite[1:numElite,1:(r1+r2)]
    message(paste("Current Iteration:", iter, "Current ABC is:", G.current[[iter]]))
  }
  all.iter.score <- as.numeric(c(G.best, G.current))
  final.model <- sort(c(as.numeric(G.current.numElite[1,1:(r1+r2)])[which(as.numeric(G.current.numElite[1,1:(r1+r2)])!=0)]))
  final.model.names <- predictor_match(final.model, p, r1, r2, interaction.ind)[[1]]
  final.model.score <- as.numeric(G.current.numElite[1,((r1+r2)+1)])

  obj <- list(Final.variable.names = final.model.names,
              Final.variable.idx = final.model,
              Final.model.score = final.model.score,
              All.iter.score = all.iter.score)

  class(obj) <- "hySAINT"
  return(obj)
}
