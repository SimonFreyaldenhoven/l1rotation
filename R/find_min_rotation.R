
spherical_to_cartesian <- function(theta){
  if(!is.matrix(theta)) {
    r <- length(theta) + 1
    R <- rep(0, r)
    R[1] <- cos(theta[1])
    if(r > 2){
      for (kk in 2:(r-1)) {
        R[kk] <- prod(sin(theta[1:(kk-1)])) * cos(theta[kk])
      }
    }
    R[r] <- prod(sin(theta))
    return(R)
  }
  stopifnot(nrow(theta) > 0)

  r <- nrow(theta) + 1
  no_draws <- ncol(theta)

  R <- matrix(0, nrow = r, ncol = no_draws)

  R[1, ] <- cos(theta[1, ])

  if(r > 2){
    for (kk in 2:(r - 1)) {
      R[kk, ] <- col_prod(sin(theta[1:(kk - 1), ]))*cos(theta[kk, ])

    }
  }

  if(r > 1){
    R[r, ] <- col_prod(sin(theta))
  }

  return(R)
}

objectivefcn_spherical <- function(theta, initial_loadings) {
  R <- spherical_to_cartesian(theta)
  sum(abs(initial_loadings %*% R))
}

gridsize <- function(factorno) {
  # defines number of random draws to start search for local minma from
  if (factorno == 2) {
    no_randomgrid <- 500
  } else if (factorno == 3) {
    no_randomgrid <- 1000
  } else if (factorno == 4) {
    no_randomgrid <- 2000
  } else if (factorno == 5) {
    no_randomgrid <- 4000
  } else if (factorno > 5 && factorno < 9) {
    no_randomgrid <- 6000
  } else {
    no_randomgrid <- 10000
  }

  return(no_randomgrid)
}


find_min_rotation <- function(initial_loadings, parallel = FALSE, n_cores = NULL) {

  stopifnot(is.matrix(initial_loadings))
  stopifnot(ncol(initial_loadings) > 1)
  if(any(is.na(initial_loadings)) | any(is.infinite(initial_loadings))) stop("initial_loadings contains missing or infinite values.")
  if(!all(is.numeric(initial_loadings))) stop("initial_loadings contains non-numeric values.")

  stopifnot((n_cores %% 1 == 0 & n_cores > 0) | is.null(n_cores))
  stopifnot(is.logical(parallel))
  if(parallel & is.null(n_cores)) stop("parallel set to TRUE but n_cores is NULL. Please specify n_cores for parallel execution.")
  if(!parallel & !is.null(n_cores)) warning("parallel set to FALSE but n_cores is not null. Defaulting to sequential execution.")


  if(parallel) cluster <- setup_cluster(n_cores)

  r <- ncol(initial_loadings)
  no_draws <- gridsize(r)
  l1_norm <- rep(0, no_draws)
  exitflag <- rep(0, no_draws)

  # Create starting points for algorithm
  initial_draws <- matrix(stats::rnorm(r * no_draws), nrow = r)
  initial_draws <- normalize(initial_draws, p = 2)

  theta <- cartesian_to_spherical(initial_draws)

  # Optimization in polar coordinates happens w.r.t. theta
  angles <- theta
  l <- nrow(angles)
  results <- list()
  `%dopar%` <- foreach::`%dopar%`
  bind_rows <- dplyr::bind_rows

  functions_to_keep <- c("col_prod", "spherical_to_cartesian", "objectivefcn_spherical")

  if(parallel) {
    results <- foreach::foreach(rep = 1:no_draws, .combine = "bind_rows", .export = functions_to_keep) %dopar% {

      starting_point <- theta[, rep]
      result <- stats::optim(
        starting_point,
        objectivefcn_spherical, initial_loadings = initial_loadings,
        control = list(maxit = 200 * l, ndeps = 1e-4, reltol = 1e-7, warn.1d.NelderMead = FALSE),
        method = 'Nelder-Mead'
      )

      result_tbl <- data.frame(rep = rep, par = result$par, l1_norm = result$value, exitflag = result$convergence)
      results[rep] <- list(result_tbl)

    }
    parallel::stopCluster(cl = cluster)

    angles <- matrix(data = results$par, nrow = l, byrow = FALSE)

    l1_norm <- results %>%
      dplyr::group_by(rep) %>%
      dplyr::slice(1) %>%
      dplyr::pull(l1_norm)

    exitflag <- results %>%
      dplyr::group_by(rep) %>%
      dplyr::slice(1) %>%
      dplyr::pull(exitflag)

  } else{

    for (rep in cli::cli_progress_along(1:no_draws, "Finding rotations")) {

      starting_point <- theta[, rep]
      result <- stats::optim(
        starting_point,
        objectivefcn_spherical, initial_loadings = initial_loadings,
        control = list(maxit = 200 * l, ndeps = 1e-4, reltol = 1e-7, warn.1d.NelderMead = FALSE),
        method = 'Nelder-Mead'
      )

      angles[, rep] <- result$par
      l1_norm[rep] <- result$value
      exitflag[rep] <- result$convergence
    }
  }

  # Convert back to cartesian coordinates
  R <- spherical_to_cartesian(angles)

  return(list(R = R, l1_norm = l1_norm, exitflag = exitflag))
}

setup_cluster <- function(n_cores){
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk)) {
    # use 2 cores in CRAN
    n_cores <- 1L
  } else {
    # use all cores in devtools::test()
  }
  cluster <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cluster)
  return(cluster)
}

normalize <- function(X, p = 2){
  stopifnot(is.matrix(X))
  norms <- apply(X, p = p, 2, pracma::Norm)
  X_norm <- sweep(X, 2, norms, FUN = "/")
  return(X_norm)
}

# Returns the norm of each column of a matrix
vecnorm <- function(X, p = 2){
  apply(X, p = p, 2, pracma::Norm)
}


col_prod <- function(data){
  if(is.matrix(data)) matrixStats::colProds(data)
  else{
    c(data)
  }
}

# Assumes radius is equal to 1 (that is, X is normalized)
cartesian_to_spherical <- function(X){
  stopifnot(nrow(X) > 1)
  r <- nrow(X)
  no_draws <- ncol(X)

  theta <- matrix(0, nrow = r - 1, ncol = no_draws)
  if(r-2 > 0){
    for (kk in 1:(r - 2)) {
      theta[kk, ] <- atan2( vecnorm(X[(kk + 1):r, ]), X[kk, ])
    }
  }
  theta[r - 1, ] <- atan2( X[r, ], X[(r - 1), ] )

  return(theta)
}

