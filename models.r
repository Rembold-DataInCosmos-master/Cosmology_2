library(tidyverse)

#' @title Calculate Comoving Distance from Redshift
#' @description Computes cosmological distances in Mpc from redshifts based on ΛCDM parameters.
#' @param params A numeric vector of cosmological parameters: `c(H0, Ω_m, Ω_k, Ω_Λ)`
#' @param zs A numeric vector of redshift values
#' @param num_samples Number of Riemann slices to use (default = 5000)
#' @return A numeric vector of distances in Mpc
#' @export
calc_distance <- function(params, zs, num_samples = 5000) {
  H0   <- params[1]
  om_m <- params[2]
  om_k <- params[3]
  om_l <- params[4]

  H <- function(z) {
    H0 * sqrt(om_m * (1 + z)^3 + om_k * (1 + z)^2 + om_l)
  }

  max_z <- max(zs)
  edges <- seq(0, max_z, length.out = num_samples)
  bin_sizes <- diff(edges)
  z_range <- edges[-length(edges)] + bin_sizes / 2

  integrand <- (1 / H(z_range)) * bin_sizes
  cum_integral <- cumsum(integrand)

  int_tbl <- tibble(z = z_range, I = cum_integral)

  tibble(z = zs) %>%
    rowwise() %>%
    mutate(
      idx = max(which(int_tbl$z < z)),
      distance = (1 + z) * 3e5 * int_tbl$I[idx]
    ) %>%
    ungroup() %>%
    pull(distance)
}

