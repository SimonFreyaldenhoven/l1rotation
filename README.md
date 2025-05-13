
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/l1rotation)](https://CRAN.R-project.org/package=l1rotation)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/l1rotation)](https://r-pkg.org/pkg/l1rotation)
<!-- badges: end -->

## Overview

The `l1rotation` package implements the l1-rotation criterion of
[Freyaldenhoven (2025)](https://doi.org/10.21799/frbp.wp.2020.25) to
simplify the loading matrix in factor models.

## Installation

``` r
# Install from CRAN
install.packages("l1rotation")

# Install latest version from GitHub
install.packages("devtools")
devtools::install_github("SimonFreyaldenhoven/l1rotation")
```

## Usage

Find a minimal example below. For more information see the package
[website](https://kobleary.github.io/l1rotation/),
[documentation](https://cran.r-project.org/package=l1rotation/l1rotation.pdf)
and
[vignette](https://cran.r-project.org/package=l1rotation/vignettes/documentation.html).

``` r
library(l1rotation)
set.seed(916) 

# Minimal example with 2 factors, where X is a 224 by 207 matrix
lf <- local_factors(X = example_data, r = 2)

# Rerun with parallel processing
lf_lambdas <- local_factors(X = example_data, r = 2, parallel = TRUE, n_cores = 10)

# Visualize Principal Component estimate of the loadings
lf$pc_plot
```

![](man/figures/README-example-1.png)<!-- -->

``` r

# Visualize l1-rotation loadings
lf$rotated_plot
```

![](man/figures/README-example-2.png)<!-- -->

## Citation

Simon Freyaldenhoven. “Identification Through Sparsity in Factor Models:
the l1-rotation criterion.” [Philadelphia Fed Working Paper
20-25](https://doi.org/10.21799/frbp.wp.2020.25), February 2025.

Simon Freyaldenhoven, Ryan Kobler. “`l1rotation` package.” Code and data
repository at <https://github.com/SimonFreyaldenhoven/l1rotation>, March
2025.
