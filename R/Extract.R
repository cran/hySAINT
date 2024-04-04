#' Extracting columns and generating required interaction effects from data\cr
#'
#' This function simplifies the data preparation process by enabling users to
#' extract specific columns from their dataset \code{X}, and automatically
#' generating any necessary interaction effects based on \code{varind}.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects. Note that the interaction
#'        effects should not be included in \code{X} because this function
#'        automatically generates the corresponding interaction effects if needed.
#' @param varind A numeric vector that specifies the indices of variables to be
#'        extracted from \code{X}. Duplicated values are not allowed. See Example
#'        for details.
#' @param interaction.ind A two-column numeric matrix. Each row represents a unique
#'        interaction pair, with the columns indicating the index numbers of the variables
#'        involved in each interaction. Note that \code{interaction.ind} must be generated
#'        outside of this function using \code{t(utils::combn(p,2))}. See Example section for
#'        details.
#'
#' @return A numeric matrix is returned.
#' @export
#'
#' @examples # Generate interaction.ind
#' interaction.ind <- t(combn(4,2))
#'
#' @examples # Generate data
#' set.seed(0)
#' X <- matrix(rnorm(20), ncol = 4)
#' y <- X[, 2] + rnorm(5)
#'
#' @examples # Extract X1 and X1X2 from X1, ..., X4
#' Extract(X, varind = c(1,5), interaction.ind)
#'
#' @examples # Extract X5 from X1, ..., X4
#' Extract(X, varind = 5, interaction.ind)
#'
#' @examples # Extract using duplicated values
#' try(Extract(X, varind = c(1,1), interaction.ind)) # this will not run
#'
Extract <- function(X, varind, interaction.ind = NULL){
  if (is.null(interaction.ind)){
    stop("Interaction.ind is missing. Use t(utils::combn()) to generate interaction matrix.")
  }
  if (as.logical(any(duplicated(varind[which(varind!=0)])))){
    stop("There cannot be duplicated values in varind.")
  }
  ncoln <- length(varind)
  nrown <- nrow(X)
  p <- ncol(X)
  mainind <- varind[which(varind%in%1:p)]
  mainvars.matrix <- X[, mainind]
  interind <- varind[which(varind>p)]
  a <- interaction.ind[interind-p,]
  if (length(a) > 1){
    intervars.matrix <- matrix(0, nrow = nrown, ncol = mydim(a)[1])
    for (i in 1:mydim(a)[1]) {
      if (mydim(a)[1]==1) intervars.matrix[,i] <- X[,a[1]]*X[,a[2]]
      else {
        intervars.matrix[,i] <- X[, a[i,1]]*X[,a[i,2]]
      }
    }
    colnames(intervars.matrix) <- paste0("X", interind)
    data_extract <- cbind(mainvars.matrix, intervars.matrix)
  }
  else{
    data_extract <- mainvars.matrix
  }
  if (length(mainind) == 1){
    data_extract <- as.matrix(data_extract)
    colnames(data_extract)[which(colnames(data_extract)=="mainvars.matrix" )] <- paste0("X", mainind)
  }
  if (length(mainind) == 1 && length(interind) ==0){
    colnames(data_extract) <- paste0("X", mainind)
  }
  if (length(mainind) == 0){
    data_extract <- data_extract
    colnames(data_extract) <- paste0("X", interind)
  }
  if (length(mainind) > 1){
    colnames(data_extract)[1:mydim(mainvars.matrix)[2]] <- paste0("X", mainind)
  }
  return(data_extract)
}

# BELOW ARE THE FUNCTIONS FUNCTIONS THAT DO NOT NEED TO BE DOCUMENTED

# MODIFIED CHOOSE FUNCTION
mychoose <- function(k1){
  if (k1 == 1){
    return(1)
  }else{
    return(choose(k1,2))
  }
}

# MODIFIED DIM FUNCTION
mydim <- function(x){
  if (is.matrix(x)) dim(x)
  else return(1)
}

# MODIFIED SORT FUNCTION
sort_zeros <- function(vec){
  non_zeros <- vec[vec != 0]
  zeros <- vec[vec == 0]
  if (length(zeros) > 0){
    return(c(sort(non_zeros), zeros))
  }else{
    return(c(sort(non_zeros)))
  }
}

# MODIFIED SAMPLE FUNCTION
mysample <- function(x, size, replace = F, prob = NULL){
  if (length(x) == 1) return(x)
  if (length(x) > 1) return(sample(x, size, replace, prob))
}

# FUNCTION USED TO REMOVE DUPLICATED ROWS
remove.dup.rows <- function(Matrix, r1, r2){
  duplicates <- duplicated(Matrix) | duplicated(Matrix, fromLast = TRUE)
  any_duplicates <- any(duplicates)
  if (any_duplicates == TRUE){
    row_dup <- as.numeric(which(duplicates))
    rm_ind <- row_dup[seq(1, length(row_dup), by = 2)]
    no_dup_new <- Matrix[-rm_ind,]
    Matrix <- no_dup_new[order(unlist(no_dup_new[,((r1+r2)+1)]), na.last = TRUE),]
  }else{
    Matrix <- Matrix
  }
  return(Matrix)
}


# FUNCTION USED TO MATCH IDX TO VARIABLE NAMES
predictor_match <- function(candidate.model, p, r1,r2,interaction.ind){
  model.match <- list()
  if (mydim(candidate.model)[1]==1){
    ee <- as.numeric(candidate.model[1:(r1+r2)])
    ee1 <- as.numeric(ee[which(ee%in% 1:p)])
    ee2 <- as.numeric(ee[which(ee > p)])
    if (length(ee2) == 0){
      model.match[[1]] <-c(paste0("X", ee1))
    }
    if (length(ee1) == 0){
      model.match[[1]] <- c( paste0("X", interaction.ind[ee2-p,1], "X", interaction.ind[ee2-p,2]))
    }
    if (!length(ee2) ==0 & !length(ee1) ==0){
      model.match[[1]] <-c(paste0("X", ee1), paste0("X", interaction.ind[ee2-p,1], "X", interaction.ind[ee2-p,2]))
    }
  }else{
    for (i in 1:nrow(candidate.model)) {
      ee <- candidate.model[i, 1:(r1+r2)]
      ee1 <- as.numeric(ee[which(ee%in% 1:p)])
      ee2 <- as.numeric(ee[which(ee > p)])
      if (length(ee2) == 0){
        model.match[[i]] <-c(paste0("X", ee1))
      }
      if (length(ee1) == 0){
        model.match[[i]] <- c( paste0("X", interaction.ind[ee2-p,1], "X", interaction.ind[ee2-p,2]))
      }
      if (!length(ee2) ==0 & !length(ee1) ==0){
        model.match[[i]] <-c(paste0("X", ee1), paste0("X", interaction.ind[ee2-p,1], "X", interaction.ind[ee2-p,2]))
      }
    }
  }
  return(model.match = model.match)
}

# FUNCTION USED TO MATCH MAIN EFFECT
mainRank_match <- function(a,b){

  ranks <- numeric(length(a))
  for (i in 1:length(a)) {
    ranks[i] <- which(b == a[i])
  }
  return(list(mainVec = b,
              mainSIS = a,
              ranks = ranks))
}

# FUNCTION USED DURING CORSSOVER
sample_prob <- function(X, myParent, EVAoutput, heredity = "Strong", r1, r2){
  p <- dim(X)[2]
  InterRank.ind.mat <- EVAoutput$ranked.intermat
  max_model_size <- dim(myParent)[2]
  n_rows <- nrow(myParent)
  top_20_percent <- ceiling(0.2*n_rows)
  row_probabilities <- rep(0.1, n_rows)
  row_probabilities[1:top_20_percent] <- 0.9
  if (n_rows == 3 | n_rows == 2){
    even_numbers <- 2
  }else{
    even_numbers <- seq(2, (n_rows-2), by = 2)
  }

  num_rows_to_sample <- mysample(even_numbers,1)
  sampled_rows <- myParent[sample(n_rows, size = num_rows_to_sample, replace = TRUE, prob = row_probabilities), , drop = FALSE]
  CPMatrix <- sampled_rows
  CMatrix <-  matrix(0, nrow=choose(dim(CPMatrix)[1],2),ncol=max_model_size)
  tempcount <- 0
  for (i in 1:(dim(CPMatrix)[1]-1)) {
    for (j in ((i+1):(dim(CPMatrix)[1]))) {
      tempcount <- tempcount + 1
      crossind <- as.numeric(unlist(c(CPMatrix[i,][which(!CPMatrix[i,]==0)],
                                      CPMatrix[j,][which(!CPMatrix[j,]==0)])))
      counts <- table(crossind)
      unique_numbers <- as.numeric(names(counts[counts == 1]))
      common_numbers <- as.numeric(names(counts[counts > 1]))
      bind <- c(common_numbers,unique_numbers)
      bind_main <- bind[bind <= p]
      bind_inter <- bind[bind > p]
      bothmain <- bind_main[which(bind_main%in%common_numbers)]
      bothinter <- bind_inter[which(bind_inter%in%common_numbers)]
      main_selection_rates <- rep(0.5, length(bind_main))
      inter_selection_rates <- rep(0.5, length(bind_inter))
      if (length(bothmain) != 0) {
        main_selection_rates <- rep(0.3, length(bind_main))
        main_selection_rates[1:length(bothmain)] <- 0.7
      }
      if (length(bothinter) != 0) {
        inter_selection_rates <- rep(0.3, length(bind_inter))
        inter_selection_rates[1:length(bothinter)] <- 0.7
      }
      sample_bind_main <- sort(mysample(bind_main, min(mysample(c(round(r1/2,0):r1),1), length(bind_main)),
                                        replace = FALSE, prob = main_selection_rates))
      sample_bind_main <- as.numeric(sort(unique(sample_bind_main)))
      sample_bind_inter <- sort(mysample(bind_inter, min(mysample(c(1:r2),1), length(bind_inter)),
                                         replace = FALSE, prob = inter_selection_rates))
      sampled_ind <- sort(c(sample_bind_main, sample_bind_inter))
      main_ind <- sort(sampled_ind[sampled_ind <= p])
      maps <- InterRank.ind.mat[InterRank.ind.mat[,3]%in% sampled_ind[which(sampled_ind > p)],]
      if (heredity == "Strong"){
        if (mydim(maps)[1] == 1){
          aaa <- maps[maps[1] %in% main_ind & maps[2] %in% main_ind]
        }else{
          aaa <- maps[maps[,1] %in% main_ind & maps[,2] %in% main_ind,]
        }
        if (mydim(aaa)[1] == 1){
          bbb <- as.numeric(unlist(aaa[3]))
        }else{
          bbb <- as.numeric(aaa[,3])
        }
        if (length(bbb) == 0){
          bbb <- 0
        }else{
          bbb <- bbb
        }
      }else if (heredity == "Weak"){
        if (mydim(maps)[1] == 1){
          aaa <- maps[maps[1] %in% main_ind | maps[2] %in% main_ind]
        }else{
          aaa <- maps[maps[,1] %in% main_ind | maps[,2] %in% main_ind,]
        }
        if (mydim(aaa)[1] == 1){
          bbb <- as.numeric(unlist(aaa[3]))
        }else{
          bbb <- as.numeric(aaa[,3])
        }
        if (length(bbb) == 0){
          bbb <- 0
        }else{
          bbb <- bbb
        }
      }else{
        bbb <- as.numeric(sampled_ind[which(sampled_ind>p)])
        if (length(bbb) == 0){
          bbb <- 0
        }else{
          bbb <- bbb
        }
      }
      fill <- c(main_ind, bbb)
      CMatrix[tempcount,c(1:length(fill))] <- fill
    }
  }
  CMatrix[is.na(CMatrix)] <- 0
  Offsprings <- CMatrix

  return(Offsprings)
}

# FUNCTION USED DURING INITIAL
sampleGpool <- function(X, EVAoutput, heredity, r1, r2){
  p <- dim(X)[2]
  n <- dim(X)[1]
  if (p < n ){
    initial.main.numbs <- EVAoutput$ranked.mainpool[1:r1]
  }else{
    initial.main.numbs <- EVAoutput$ranked.mainpool
  }
  initial.inter.pool <- EVAoutput$ranked.intermat

  initial.main.select <- initial.main.numbs[1:r1]
  if (heredity == "Strong"){
    mapinter.initial.main.select <- initial.inter.pool[initial.inter.pool[, 1] %in% initial.main.select & initial.inter.pool[, 2] %in% initial.main.select, ]
  }else if (heredity == "Weak"){
    mapinter.initial.main.select <- initial.inter.pool[initial.inter.pool[, 1] %in% initial.main.select | initial.inter.pool[, 2] %in% initial.main.select, ]
  }else if (heredity == "No"){
    mapinter.initial.main.select <- initial.inter.pool
  }
  initial.inter.numbs <- as.numeric(mapinter.initial.main.select[,3])
  top_1_percent <- ceiling(0.2*length(initial.inter.numbs))
  # top_1_percent <- ceiling(0.01*length(initial.inter.numbs))
  inter_filler_probabilities <- rep(0.1, length(initial.inter.numbs))
  inter_filler_probabilities[1:top_1_percent] <- 0.9
  initial.inter.select <- mysample(initial.inter.numbs, 1, replace = FALSE, prob = inter_filler_probabilities)

  filler.initial <- c(initial.main.select, initial.inter.select)
  return(filler.initial)
}

# FUNCTION USED TO GENERATE INTERACTION POOL FOR DIFFERENT HEREDITY CONDITIONS
inter_pool_generate <- function (a, p, heredity = "Strong", interaction.ind){
  if (heredity == "No"){
    interpooltemp <- interaction.ind
  }
  else if (heredity == "Weak"){
    df <- rbind(interaction.ind[interaction.ind[,1] %in% a,][order(stats::na.exclude(match(interaction.ind[,1], a))),],
                interaction.ind[interaction.ind[,2] %in% a,][order(stats::na.exclude(match(interaction.ind[,2], a))),])
    interpooltemp <- df[!duplicated(df),]
  }
  else if (heredity == "Strong"){
    interpooltemp <- t(utils::combn(sort(a),2))
  }
  intercandidates.ind <- match(do.call(paste, as.data.frame(interpooltemp)), do.call(paste, as.data.frame(interaction.ind)))+ p
  intercandidates.mat <- as.matrix(cbind(interpooltemp,intercandidates.ind))
  colnames(intercandidates.mat) <- c("a", "b", "idx")
  return(intercandidates.mat)
}

# FUNCTION USED TO EVALUATE INTERACTION
interaction_evaluation <- function(a, b, X, y, heredity = "Strong", sigma = NULL, varind = NULL, interaction.ind = NULL, lambda = 10){
  p <- dim(X)[2]
  if (heredity == "Strong" | heredity == "Weak"){
    interscore <- lapply(1:length(b), function(i) {
      ABC(X, y, heredity, sigma, varind = c(a, b[i]),interaction.ind, lambda)})
  }else if(heredity == "No"){
    if (length(b) <= choose(1000,2)){
      interscore <- lapply(1:length(b), function(i) {
        ABC(X, y, heredity, sigma, varind = c(a, b[i]),interaction.ind, lambda)})
    }else{
      interscore <- lapply(1:length(b), function(i){
        abs(as.numeric(energy::dcor(Extract(X, varind = c(b[i]), interaction.ind), y)))
      })
    }
  }
  return(interscore)
}

# FUNCTION USED TO EVALUATE MY ENTIRE FAMILY
myfamily.eva <- function(Myfamily, r1, r2, X , y, heredity = "Strong",
                         sigma = NULL, varind = NULL, interaction.ind = NULL , lambda = 10){
  G.current <- lapply(1:nrow(Myfamily), function(i) {
    ABC(X , y, heredity, sigma, varind = c(as.numeric(Myfamily[i,])[as.numeric(Myfamily[i,]) != 0]),interaction.ind,lambda)})
  G.current.matrix <- as.matrix(cbind(Myfamily, G.current))
  colnames(G.current.matrix) <- NULL
  if (mydim(G.current.matrix)[1] == 1){
    G.current.matrix.ordered <- G.current.matrix
  }else{
    G.current.matrix.ordered <-  G.current.matrix[order(unlist(G.current.matrix[,((r1+r2)+1)]), na.last = TRUE),]
  }
  return(G.current.matrix.ordered)
}

# FUNCION USED TO PERFORM SIMULATED ANNEALING
Saint <- function(myParent, EVAoutput, r1, r2, initial.temp = 1000, cooling.rate = 0.95, X, y,
                  heredity = "Strong", sigma = NULL, varind = NULL, interaction.ind = NULL, lambda = 10){
  p <- dim(X)[2]

  # CALCULATE CURRENT ABC FOR MYPARENT
  E.best <- lapply(1:nrow(myParent), function(i){
    ABC(X , y, heredity, sigma, varind = c(as.numeric(myParent[i,])[as.numeric(myParent[i,]) != 0]),interaction.ind , lambda)})

  mutant <- matrix(0, nrow = mydim(myParent)[1], ncol = (r1+r2))
  for (i in 1:nrow(myParent)) {
    if (i == 1){current.temp <- initial.temp}

    # DETERMINE CURRENT MAIN
    current.main.numbs <- as.numeric(myParent[i,][which(myParent[i,]%in% 1:p)])
    current.main.size <- length(current.main.numbs)
    current.inter.numbs <- as.numeric(myParent[i,][which(myParent[i,]> p)])

    # DETERMINE WHETHER MAIN EFFECTS IS ADDED OR SUBTRACTED
    # IF 0 THEN DELETE. IF 1 THEN ADD
    subtract.or.add <- mysample(0:1,1)

    # IF 1 IS SELECTED. AND THE MAIN EFFECT IS LESS THAN R1
    # ALLOW ADD
    if (subtract.or.add == 1 & current.main.size < r1){
      not.in.current.main.numbs <- as.numeric(setdiff(EVAoutput$ranked.mainpool, current.main.numbs))

      if (length(not.in.current.main.numbs) == 0){
        new.main.numbs <- current.main.numbs
        new.main.size <- length(new.main.numbs)
      }else{
        new.sample.main.numbers <- mysample(not.in.current.main.numbs, 1, replace = FALSE)
        new.main.numbs <- c(current.main.numbs, new.sample.main.numbers)
        new.main.size <- length(new.main.numbs)
      }

      # IF 1 IS SELECTED. AND THE MAIN EFFECT IS EQUAL TO R1
      # FORCE DELETE
    }else if (subtract.or.add == 1 & current.main.size == r1){
      remove.numbs <- mysample(current.main.numbs, 1, replace = FALSE)
      new.main.numbs <- as.numeric(setdiff(current.main.numbs, remove.numbs))
      new.main.size <- length(new.main.numbs)

      # IF 0 IS SELECED. AND THE MAIN EFFECT WILL BE GREATER THAN ZERO
      # ALLOW DELETE
    }else if (subtract.or.add == 0 & ((current.main.size - 1) >= 1)){
      remove.numbs <- mysample(current.main.numbs, 1, replace = FALSE)
      new.main.numbs <- as.numeric(setdiff(current.main.numbs, remove.numbs))
      new.main.size <- length(new.main.numbs)

      # IF 0 IS SELECTED. AND THE MAIN EFFECT WILL BE LOWER TO ZERO
      # ALLOW DELETE FOR NO HEREDITY
      # FORCE ADD FOR STRONG AND WEAK HEREDITY
    } else if (subtract.or.add == 0 & ((current.main.size - 1) < 1)){
      if (heredity == "No"){
        new.main.numbs <- NULL
      }else{
        not.in.current.main.numbs <- as.numeric(setdiff(EVAoutput$ranked.mainpool, current.main.numbs))
        new.sample.main.numbers <- mysample(not.in.current.main.numbs, 1, replace = FALSE)
        new.main.numbs <- c(current.main.numbs, new.sample.main.numbers)
        new.main.size <- length(new.main.numbs)
      }
    }

    # CALCULATE CURRENT ABC FOR NEW MAIN.
    # THIS IS THE ABC SCORE AFTER NEW MAIN NUMBS ADDED
    E.current.m <- ABC(X, y, heredity, sigma, varind = c(new.main.numbs, current.inter.numbs),interaction.ind, lambda)

    # IF E CURRENT M < E BEST THEN ACCEPT CURRENT MAIN
    # IF E CURRENT M IS NOT BETTER. THEN DO SA
    if (as.numeric(E.current.m) < as.numeric(E.best[[i]])){
      best.current.main.numbs <- new.main.numbs
    }else{
      prob.value <- exp((as.numeric(E.best[[i]]) - as.numeric(E.current.m))/current.temp)
      rand.value <- runif(1)
      if (prob.value <= rand.value) {best.current.main.numbs <- new.main.numbs}else{
        best.current.main.numbs <- current.main.numbs
      }
    }

    # GENERATE NEW INTERACTION POOL. AND CHECK WHETHER IF HEREDITY CONDITION IS
    # SATISFIED FOR CURRENT INTERACTION
    new.interaction.pool <- inter_pool_generate(best.current.main.numbs, p, heredity, interaction.ind)
    new.interaction.pool.idx <- as.numeric(new.interaction.pool[,3])

    if (heredity == "Strong" | heredity == "Weak"){
      current.inter.numbs.h <- intersect(current.inter.numbs, as.numeric(new.interaction.pool[,3]))
      if (length(current.inter.numbs.h) > 0){
        current.inter.numbs <- current.inter.numbs.h
      }else{
        current.inter.numbs <- NULL
      }
    }else{
      current.inter.numbs <- current.inter.numbs
    }
    current.inter.size <- length(current.inter.numbs)

    # DETERMINE WHETHER INTERACTION EFFECTS IS ADDED OR SUBTRACTED
    # IF 0 THEN DELETE. IF 1 THEN ADD
    inter.subtract.or.add <- mysample(0:1,1)

    # IF CURRENT INTER NUMBER IS NULL
    if (is.null(current.inter.numbs)){
      new.inter.numbs <- mysample(new.interaction.pool.idx, 1, replace = FALSE)

      # IF CURRENT INTER NUMBER IS NOT NULL. DO SAMPLE
    }else{

      # IF 1 IS SELECTED. AND THE INTERACTION EFFECT IS LESS THAN R2
      # ALLOW ADD
      if (inter.subtract.or.add == 1 & current.inter.size < r2){
        not.in.current.inter.numbs <- setdiff(new.interaction.pool.idx, current.inter.numbs)
        if (length(not.in.current.inter.numbs) == 0){
          new.inter.numbs <- current.inter.numbs
          new.inter.size <- length(new.inter.numbs)
        }else{
          new.sample.inter.numbers <- mysample(not.in.current.inter.numbs, 1, replace = FALSE)
          new.inter.numbs <- c(current.inter.numbs, new.sample.inter.numbers)
          new.inter.size <- length(new.inter.numbs)
        }

        # IF 1 IS SELECTED. AND THE INTERACTION EFFECT IS EQUAL TO R2
        # FORCE DELETE
      }else if (inter.subtract.or.add == 1 & current.inter.size == r2){
        remove.numbs <- mysample(current.inter.numbs, 1, replace = FALSE)
        new.inter.numbs <- setdiff(current.inter.numbs, remove.numbs)
        new.inter.size <- length(new.inter.numbs)

        # IF 0 IS SELECED. AND THE INTERACTION EFFECT WILL BE GREATER THAN ZERO
        # ALLOW DELETE
      }else if (inter.subtract.or.add == 0 & ((current.inter.size - 1) >= 1)){
        remove.numbs <- mysample(current.inter.numbs, 1, replace = FALSE)
        new.inter.numbs <- setdiff(current.inter.numbs, remove.numbs)
        new.inter.size <- length(new.inter.numbs)

        # IF 0 IS SELECTED. AND THE INTERACTION EFFECT WILL BE LOWER TO ZERO
        # FORCE ADD
      } else if (inter.subtract.or.add == 0 & ((current.inter.size - 1) < 1)){
        not.in.current.inter.numbs <- setdiff(new.interaction.pool.idx, current.inter.numbs)
        if (length(not.in.current.inter.numbs) == 0){
          new.inter.numbs <- current.inter.numbs
          new.inter.size <- length(new.inter.numbs)
        }else{
          new.sample.inter.numbers <- mysample(not.in.current.inter.numbs, 1, replace = FALSE)
          new.inter.numbs <- c(current.inter.numbs, new.sample.inter.numbers)
          new.inter.size <- length(new.inter.numbs)
        }
      }
    }

    # CALCULATE CURRENT ABC FOR NEW MAIN AND NEW INTER.
    # THIS IS THE ABC SCORE AFTER NEW MAIN NUMB AND NEW INTER NUMB ADDED
    E.current <- ABC(X, y, heredity, sigma, varind = c(best.current.main.numbs, new.inter.numbs),interaction.ind, lambda)

    # IF E CURRENT < E.BEST THEN ACCEPT CURRENT INTER
    # IF E CURRENT IS NOT BETTER. THEN SA
    if (as.numeric(E.current) < as.numeric(E.best[[i]])){
      best.current.inter.numbs <- new.inter.numbs
    }else{
      prob.value <- exp((as.numeric(E.best[[i]]) - as.numeric(E.current))/current.temp)
      rand.value <- runif(1)
      if (prob.value <= rand.value) {best.current.inter.numbs <- new.inter.numbs}else{
        best.current.inter.numbs <- current.inter.numbs
      }
    }

    filler <- c(best.current.main.numbs, best.current.inter.numbs)
    mutant[i, 1:length(filler)] <- filler
    current.temp <- current.temp*cooling.rate
  }
  return(mutant)
}


