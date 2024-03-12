#' Evaluating main and interaction effects\cr
#'
#' This function ranks each main and interaction effect. It also calculate the ABC
#'    score for each potential interactions across different heredity structures.
#'    If \code{heredity = "No"} and the the number of potential interactions exceed
#'    \code{choose(1000,2)}, distance correlation between each variable in \code{X}
#'    and \code{y} will be calculated so that it reduces the running time.
#'    This ensures a more efficient evaluation process.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects.
#' @param y Response variable. A \code{n}-dimensional vector.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param r1 At most how many main effects do you want to include in your model?.
#'        For high-dimensional data, \code{r1} cannot be larger than the number of screened main effects.
#' @param sigma The standard deviation of the noise term. In practice, sigma is usually
#'        unknown. Users can estimate sigma from function \code{selectiveInference::estimateSigma},
#'        then use the output as the sigma value.
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
#' @return A list of output. The components are: ranked main effect, \code{ranked.mainpool};
#' and a 4-column matrix contains potential interactions ranked by ABC score, \code{ranked.intermat}.
#' @export
#'
#' @seealso \code{\link{ABC}}, \code{\link{Extract}}
#'
#' @examples # Strong heredity
#' set.seed(0)
#' interaction.ind <- t(combn(10,2))
#' X <- matrix(rnorm(100*10,1,0.1), 100, 10)
#' epl <- rnorm(100,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,3]+X[,1]*X[,2]+X[,1]*X[,3]+epl
#' EVAoutput <- EVA(X, y, r1 = 5, sigma = 0.01, interaction.ind = interaction.ind)

EVA <- function(X, y, heredity = "Strong", r1, sigma, varind = NULL, interaction.ind = NULL,lambda = 10){

  n <- dim(X)[1]
  p <- dim(X)[2]

  # MAIN EFFECT SELECTION
  mainData <- Extract(X, varind = c(1:(dim(X)[2])), interaction.ind)
  screenIID_output <- VariableScreening::screenIID(mainData, y, method = "SIRS")
  mainVec <- as.numeric(order(screenIID_output$measurement, decreasing = TRUE))

  if (p < n){
    # MAIN EFFECT USED FOR INTERACTION EVALUATION
    sis.ix0 <- mainVec[1:r1]
    # MAIN EFFECT RANKED POOL INDEX
    maincandidates.pool <- mainVec
  }else{
    # MAIN EFFECT USED FOR INTERACTION EVALUATION
    sis_output <-  SIS::SIS(mainData, y, penalty = "SCAD")
    sis.ix0 <- sis_output$sis.ix0
    # MAIN EFFECT RANKED POOL INDEX
    mainpool <- mainRank_match(sis.ix0, mainVec)
    maincandidates.pool <- sis.ix0[order(mainpool$ranks)]
  }

  # CREATE INTERACTION POOL FOR ALL POTENTIAL INTERACTIONS
  # THE INTERACTION POOL IS DIFFERENT UNDER DIFFERENT HEREDITY CONDITION
  interpooltemp <- inter_pool_generate(sis.ix0, p, heredity, interaction.ind)
  intercandidates.ind <- as.numeric(unlist(interpooltemp[,3]))

  # ABC EVALUATION FOR INTERACTION
  interscore <- interaction_evaluation(sis.ix0, intercandidates.ind, X, y, heredity, sigma, varind, interaction.ind, lambda)

  # INTERACTION MATRIX. COLUMNS ARE A B INTERACTIONIDX INTERACTIONSCORE
  interaction_matrix <- as.matrix(cbind(interpooltemp, Scroes = interscore))
  if (mydim(interaction_matrix)[1] == 1){
    interaction_matrix_ranked <- interaction_matrix
  }else{
    interaction_matrix_ranked <- interaction_matrix[order(unlist(interaction_matrix[,4]), na.last = TRUE),]
  }
  return(list(
    # INTERACTION MATRIX RANKED
    # FOR P < N CASE THIS INTERACTION MATRIX RANKED ONLY CONTAIN CHOOSE(R1,2) NUMBER
    # OF INTERACTIONS. NEW INTERACTION MATRIX WILL BE GENERATED DURING MUTATION DEPEND
    # ON IF MAIN EFFECTS OUTSIDE OF R1 BEING ADDED
    ranked.mainpool = maincandidates.pool,
    ranked.intermat = interaction_matrix_ranked
  ))
}

