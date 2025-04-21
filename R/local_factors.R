utils::globalVariables(c("column", "value"))

#' Check whether local factors are present and find the rotation of the loading matrix with the smallest l1-norm.
#'
#' @description
#' `local_factors` tests whether local factors are present and returns both the Principal Component estimate of the loadings and the rotation of the loadings with the smallest l1-norm. It also produces graphical illustrations of the results.
#'
#' @param X A (usually standardized) t by n matrix of observations.
#' @param r An integer denoting the number of factors in X.
#' @param parallel A logical denoting whether the algorithm should be run in parallel.
#' @param n_cores An integer denoting how many cores should be used, if parallel == TRUE.
#'
#' @returns Returns a list with the following components:
#'  * `has_local_factors` A logical equal to `TRUE` if local factors are present.
#'  * `initial_loadings` Principal component estimate of the loading matrix.
#'  * `rotated_loadings` Matrix that is the rotation of the loading matrix that produces the smallest l1-norm.
#'  * `rotation_diagnostics` A list containing 3 components:
#'      * `R` Rotation matrix that when used to rotate `initial_loadings` produces the smallest l1-norm.
#'      * `l1_norm` Vector of length `r` containing the value of the l1 norm each solution generates.
#'      * `sol_frequency` Vector of length `r` containing the frequency in the initial grid of each solution.
#'  * `pc_plot` Tile plot of the Principal Component estimate of the loading matrix.
#'  * `rotated_plot` Tile plot of the l1-rotation of the loading matrix estimate.
#'  * `small_loadings_plot` Plot of the number of small loadings for each column of the l1-rotation of the loading matrix estimate.
#'
#' @export
#'
#'
#' @examples
#' # Minimal example with 2 factors, where X is a 224 by 207 matrix
#' lf <- local_factors(X = example_data, r = 2)
#'
#' # Visualize Principal Component estimate of the loadings
#' lf$pc_plot
#'
#' # Visualize l1-rotation loadings
#' lf$pc_rotated_plot
#'
local_factors <- function(X, r, parallel = FALSE, n_cores = NULL) {

  stopifnot(is.matrix(X) | is.data.frame(X))
  if("data.frame" %in% class(X)) X <- as.matrix(X)
  stopifnot(is.numeric(r))
  stopifnot(r %% 1 == 0 & r > 0)
  stopifnot(ncol(X) >= r)

  if(any(is.na(X)) | any(is.infinite(X))) stop("X cannot contain missing or infinite values.")
  if(!all(is.numeric(X))) stop("X cannot contain non-numeric values.")

  stopifnot(is.numeric(n_cores) | is.null(n_cores))
  if(is.numeric(n_cores)) stopifnot(n_cores %% 1 == 0 & n_cores > 0)
  stopifnot(is.logical(parallel))
  if(parallel & is.null(n_cores)) stop("parallel set to TRUE but n_cores is NULL. Please specify n_cores for parallel execution.")
  if(!parallel & !is.null(n_cores)) warning("parallel set to FALSE but n_cores is not null. Defaulting to sequential execution.")


  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  eig_X <- pca$d^2
  initial_loadings <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  rotn_result <- find_local_factors(X, r, parallel = parallel, n_cores = n_cores)
  test_result <- test_local_factors(X, r, loadings = rotn_result$rotated_loadings)
  has_local_factors <- test_result$has_local_factors
  rotated_loadings <- rotn_result$rotated_loadings

  # Illustrate loading matrices
  pc_plot <- plot_loading_matrix(initial_loadings, xlab = "k", title = "Principal Component estimate")
  rotated_plot <- plot_loading_matrix(rotated_loadings, xlab = "k", title = "Rotated estimate (l1-criterion)")
  small_loadings_plot <- plot_small_loadings(test_result, r)
  return(list(
    has_local_factors = has_local_factors,
    initial_loadings = initial_loadings,
    rotated_loadings = rotated_loadings,
    rotation_diagnostics = rotn_result$diagnostics,
    pc_plot = pc_plot,
    rotated_plot = rotated_plot,
    small_loadings_plot = small_loadings_plot))
}

plot_loading_matrix <- function(data, xlab = "", ylab = "", title = ""){

  if(is.matrix(data)) data <- convert_mat_to_df(data)

  scale_fill_pal <- select_palette(data$value, type = "difference")

  ggplot2::ggplot(data, ggplot2::aes(column, as.numeric(row), fill = value)) +
    ggplot2::geom_tile() +
    scale_fill_pal +
    ggplot2::labs(x = xlab, y = ylab, title = title, fill = "") +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      barwidth = 1,
      barheight = 15,
      label.theme = ggplot2::element_text(size = 12),
      draw.ulim = TRUE,
      draw.llim = TRUE)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
      plot.background = ggplot2::element_rect(fill = "white", color = "white")
    )
}

plot_small_loadings <- function(result, r, xlab = "k", ylab = "", title = ""){

  n_small <- result$n_small
  gamma <- result$gamma
  h_n <- result$h_n

  data.frame(value = result$n_small) %>%
    dplyr::mutate(factor = 1:r) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = factor, y = value), size = 3) +
    ggplot2::geom_hline(yintercept = gamma, linetype = "dashed", linewidth = 1) +
    ggplot2::ylim(c(min(gamma - 10, min(n_small - 5)), max(gamma + 5, max(n_small) + 5))) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::xlim(c(1, r)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"))
}


convert_mat_to_df <- function(mat){

  r <- ncol(mat)
  n <- nrow(mat)

  df <- data.frame(
    value = as.numeric(mat), row = rep(1:n, r),
    column = as.numeric(sapply(1:r, rep, times = n))
  ) %>%
    dplyr::mutate(column = factor(column, levels = 1:ncol(mat)))

  return(df)
}

select_palette <- function(range, type, breaks = NULL){

  # Define midpoint
  mp <- switch(type,
               difference = 0,
               ratio = 1,
               level = stats::median(range)
  )

  min <- min(range, na.rm = TRUE)
  max <- max(range, na.rm = TRUE)

  # Get limits
  d <- max(abs(mp - min), abs(max - mp))
  if (min < mp & max > mp) {
    auto_limits <- c(mp - d, mp + d)
    colors <- c("#008EFF", "#EBF0F4", "#961046")
  }
  if (min >= mp){
    auto_limits <- c(mp, max)
    colors <- c("#EBF0F4", "#961046")

  }
  if (max < mp) {
    auto_limits <- c(mp - d, mp)
    colors <- c("#008EFF", "#EBF0F4")

  }

  palette <- ggplot2::scale_fill_gradient2(
    high = scales::muted("darkblue"), mid = "white", low = "maroon",
    limits = auto_limits, midpoint = mp, oob = scales::squish
  )
  return(palette)
}



