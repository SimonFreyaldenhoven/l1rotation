utils::globalVariables(c("l1_norm", "non_outlier", ".data", "x"))

collate_solutions <- function(rmat_min, initial_loadings, X) {

  stopifnot(nrow(rmat_min) == ncol(initial_loadings))
  stopifnot(is.matrix(rmat_min))
  stopifnot(is.matrix(initial_loadings))
  stopifnot(is.matrix(X))

  n <- nrow(initial_loadings)
  factorno <- nrow(rmat_min)
  no_randomgrid <- ncol(rmat_min)
  epsilon_rot <- 0.05

  l1_min <- colSums(abs(initial_loadings %*% rmat_min))
  sort_index <- order(l1_min)
  l1_min_sort <- l1_min[sort_index]
  rmat_min_sort <- rmat_min[, sort_index]
  rmat_min_sort <- rmat_min_sort * pracma::repmat(sign(rmat_min_sort[1, ]), factorno, 1)

  distances <- calculate_pairwise_distances(rmat_min_sort, l1_min_sort, epsilon_rot, factorno)

  candidates <- distances$rmat_min_sort %>%
    matrix_to_dataframe() %>%
    dplyr::mutate(l1_norm = distances$l1_min_sort)

  candidates <- stats::aggregate(rep(1, nrow(candidates)), by = as.list(candidates), FUN = sum) %>%
    dplyr::rename(n = x) %>%
    dplyr::arrange(l1_norm) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(non_outlier = n/gridsize(factorno) >= 0.005)

  # Construct R from the candidates
  rmat_min_unique <- candidates %>%
    dplyr::filter(non_outlier) %>%
    dplyr::select(-c(l1_norm, n, non_outlier)) %>%
    dataframe_to_matrix()

  h_n <- 1/log(n)
  amount_sparsity <- colSums(abs(initial_loadings %*% rmat_min_unique) < h_n)

  candidates <- candidates %>%
    dplyr::filter(non_outlier) %>%
    dplyr::mutate(l0_norm = amount_sparsity)

  # up until here, candidates are arranged by l1 norm
  consolidated_mins <- consolidate_local_mins(initial_loadings, candidates, sorting_column = "l0_norm")

  consolidated_mins <- fill_with_pc(consolidated_mins, initial_loadings, X, factorno)

  rotated_loadings <- consolidated_mins$rotated_loadings
  R <- consolidated_mins$R

  loadings <- matrix_to_dataframe(R) %>%
    dplyr::left_join(candidates, by = paste0("V", 1:factorno)) %>%
    dplyr::mutate(
      l1_norm = c(l1_norm[!is.na(l1_norm)], consolidated_mins$l1_norm_update),
      n = ifelse(is.na(n), 0, n)
    )

  return(
    list(
      diagnostics = list(R = R, fval = loadings$l1_norm, sol_frequency = loadings$n),
      rotated_loadings = consolidated_mins$rotated_loadings
    )
  )

}


calculate_pairwise_distances_old <- function(rmat_min_sort, l1_min_sort, epsilon_rot, factorno) {
  no_randomgrid <- ncol(rmat_min_sort)

  # Pre-allocate norms matrix
  norms <- matrix(0, nrow = no_randomgrid, ncol = no_randomgrid)

  # Original nested loop implementation
  for (i in 1:no_randomgrid) {
    for (j in i:no_randomgrid) {
      norms[i,j] <- (pracma::Norm(rmat_min_sort[, i] - rmat_min_sort[, j], p = 2) / sqrt(factorno))
      if (norms[i,j] < epsilon_rot) {
        rmat_min_sort[, j] <- rmat_min_sort[, i]
        l1_min_sort[j] <- l1_min_sort[i]
      }
    }
  }

  list(norms = norms, rmat_min_sort = rmat_min_sort, l1_min_sort = l1_min_sort)
}

calculate_pairwise_distances <- function(rmat_min_sort, l1_min_sort, epsilon_rot, factorno) {
  no_randomgrid <- ncol(rmat_min_sort)

  # Pre-allocate results matrix
  norms <- matrix(0, nrow = no_randomgrid, ncol = no_randomgrid)

  # Calculate cross-product matrix
  cross_prod <- crossprod(rmat_min_sort)
  diag_vals <- diag(cross_prod)

  # Calculate distances using matrix operations
  # ||a - b||^2 = ||a||^2 + ||b||^2 - 2<a,b>
  for(i in 1:no_randomgrid) {
    for(j in i:no_randomgrid) {
      norms[i,j] <- sqrt(round(diag_vals[i] + diag_vals[j] - 2 * cross_prod[i,j], digits = 15)) / sqrt(factorno)
      if(norms[i,j] < epsilon_rot) {
        rmat_min_sort[,j] <- rmat_min_sort[,i]
        l1_min_sort[j] <- l1_min_sort[i]
      }
    }
  }

  list(norms = norms, rmat_min_sort = rmat_min_sort, l1_min_sort = l1_min_sort)
}


dataframe_to_matrix <- function(dataframe) {
  dataframe %>% as.matrix() %>% t()
}

matrix_to_dataframe <- function(matrix){
  matrix %>% t() %>% as.data.frame()
}

consolidate_local_mins <- function(initial_loadings, candidates, sorting_column = "l0_norm") {

  # Sorting column determines what is used to pick local minima
  # n - number of occurences
  # l0_norm - approximate l0-norm

  if (sorting_column != "l1_norm") {

    candidates <- candidates %>%
      dplyr::arrange(dplyr::desc(.data[[sorting_column]]))

    factor_cols <- grepl("^V", names(candidates))

    rmat_min_unique <- candidates[factor_cols] %>%
      dataframe_to_matrix()
  }

  # Consolidate candidates
  rotated_loadings <- initial_loadings %*% rmat_min_unique[, 1] # First candidate
  R <- rmat_min_unique[, 1]

  factorno <- nrow(rmat_min_unique)
  upperK <- ncol(rmat_min_unique)
  n <- nrow(initial_loadings)

   for (kk in 2:upperK) {

     temp <- cbind(rotated_loadings, initial_loadings %*% rmat_min_unique[, kk])
     non_singular <- min(eigen(t(temp) %*% temp)$values / n) > sqrt(1 / factorno) / 3 &&
       min(eigen(t(temp) %*% temp)$values / n) > min(eigen(t(rotated_loadings) %*% rotated_loadings)$values / n) / 4

     if(non_singular){
       rotated_loadings <- temp
       R <- cbind(R, rmat_min_unique[, kk])
     }
   }

  return(list(rotated_loadings = rotated_loadings, R = R))

}


fill_with_pc <- function(consolidated_mins, initial_loadings, X, factorno){

  rotated_loadings <- consolidated_mins$rotated_loadings
  R <- consolidated_mins$R

  candidateno <- ncol(rotated_loadings)
  I <- diag(factorno)
  l1_norm_update <- c()
  n <- nrow(X)
  t <- ncol(X)

  temp_F <-(X %*% initial_loadings/n) %*% solve(t(initial_loadings) %*% initial_loadings/n)
  eig_x <- diag(t(temp_F) %*% temp_F) * (n/t)

  while (candidateno < factorno) {
    message("Supplementing with PCs...")
    min_eig <- rep(0, factorno)
    for (ell in 1:factorno) {
      temp <- cbind(rotated_loadings, initial_loadings[, ell])
      min_eig[ell] <- min(eigen(t(temp) %*% temp)$values)
    }
    index <- which.max(min_eig * sqrt(eig_x[1:factorno]))
    rotated_loadings <- cbind(rotated_loadings, initial_loadings[, index])
    R <- cbind(R, I[, index])
    message(paste("PC used:", index))
    l1_norm_update <- c(l1_norm_update, sum(abs(initial_loadings[, index])))

    candidateno <- ncol(rotated_loadings)
  }

  return(list(rotated_loadings = rotated_loadings, R = R, l1_norm_update = l1_norm_update))

}
