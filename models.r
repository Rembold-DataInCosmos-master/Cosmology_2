library(tidyverse)

#' @title Calculate Orbital Velocity from Density Profile
#' @description Computes circular orbital velocity at specified distances based on a given density function.
#' @param params A numeric vector of density parameters passed to `density_fn`
#' @param dist A numeric vector of distances (in km)
#' @param density_fn A function of form `function(params, r)` returning density in kg/m³ at radii `r`
#' @param num_samples Number of Riemann slices (default = 5000)
#' @param spacing Either `"log"` (default) or `"linear"` for sampling
#' @return A numeric vector of velocities (in m/s)
#' @export
calc_velocity <- function(params, dist, density_fn, num_samples = 5000, spacing = "log") {
  G <- 6.67e-11
  max_dist <- max(dist)

  # Create edges
  if (spacing == "log") {
    edges <- 10^seq(log10(1e-6), log10(max_dist), length.out = num_samples)
  } else if (spacing == "linear") {
    edges <- seq(1e-6, max_dist, length.out = num_samples)
  } else {
    stop("Unknown value in `spacing` parameter. Use 'log' or 'linear'.")
  }

  bin_sizes <- diff(edges)
  radii <- edges[-length(edges)] + bin_sizes / 2

  densities <- density_fn(params, radii)
  masses <- 4 * pi * densities * radii^2 * bin_sizes
  cum_masses <- cumsum(masses)

  mass_tbl <- tibble(r = radii, m = cum_masses)

  tibble(d = dist) %>%
    rowwise() %>%
    mutate(
      idx = max(which(mass_tbl$r < d)),
      velocity = sqrt(G * mass_tbl$m[idx] / d)
    ) %>%
    ungroup() %>%
    pull(velocity)
}

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

