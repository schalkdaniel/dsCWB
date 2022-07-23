#' Get linear predictor from B-spline.
#'
#' @param x [numeric] Vector of x values
#' @param bs_dim [integer(1)] Number of base functions for the spline
#'   (default = 10L). (Corresponds to number of inner knots for
#'   the spline).
#' @param sigma [numeric(1)] Standard deviation for the normally
#'   distributed random variable from which the parameter are
#'   drawn (default = 3).
#' @param offset [numeric(1)] Shift on the y-axis of the linear
#'   predictor (default = 0).
#' @param stz [logical(1)] Sum to zero constraint for the spline.
#' @return A list with \code{x}, \code{y}, and specific information
#'   about the spline such as parameter.

x = sort(runif(100))
bs_dim = 10L
sigma = 3
offset = 0
stz = TRUE

spl = simSpline(x, stz = TRUE)

plot(x = spl$x, y = spl$y, type = "l")

simSpline = function(x, bs_dim = 10L, sigma = 3, offset = 0, stz = FALSE) {
  checkmate::assertNumeric(x = sigma, len = 1L)
  checkmate::assertNumeric(x = offset, len = 1L)
  #if (bs_dim < 7) stop("Need bs_dim >= 7 !")

  nk = bs_dim - 2

  xu = max(x)
  xl = min(x)

  xr = xu - xl
  xl = xl - xr * 0.001
  xu = xu + xr * 0.001

  dx = (xu - xl)/(nk - 1)
  kn = seq(xl - dx * 3, xu + dx * 3, length = nk + 4 + 2)

  # create the spline basis functions
  X = splines::spline.des(kn, x, 4, x * 0)$design
  if (stz)
    X = X %*% cpsp::getSubtractionRotation(X, cbind(rep(1, nrow(X))))

  # multiply with random coefficients to get random functions
  coefs = rnorm(ncol(X), sd = sigma)

  return(list(y = X %*% coefs + offset, x = x, X = X, offset = offset, coefs = coefs, knots = kn))
}

#' Simulate feature with site specific effects
#'
#' @param n       [integer(1)] Number of observations.
#' @param nsites [integer(1)] Number sites (randomly drawn).
#' @param from    [numeric(1)] Lower boundary of the feature.
#' @param up      [numeric(1)] Upper boundary of the feature.
#' @return A list with \code{x}, \code{y}, the site, and specific
#'   information about all splines (for main and site effects).
simNumeric = function(n = 1000L, nsites = 3L, from = NULL, to = NULL, offset = 0, ...) {
  checkmate::assertIntegerish(x = n, len = 1L)
  checkmate::assertIntegerish(x = nsites, len = 1L)

  if (is.null(from)) from = runif(1, -100, 100)
  if (is.null(to)) to = from + runif(1, 0, 100)

  checkmate::assertNumeric(x = from, len = 1L, upper = to)
  checkmate::assertNumeric(x = to, len = 1L, lower = from)

  x = seq(from = from, to = to, length.out = n)
  main_effect = simSpline(x, offset = offset, ...)

  idx_sites = sample(seq_len(nsites), n, TRUE)
  site_effects = lapply(seq_len(nsites), function(i) simSpline(x[idx_sites == i],
    offset = 0, stz = TRUE, ...))

  y = main_effect$y
  for (i in seq_len(nsites)) {
    site_effects[[i]]$y = site_effects[[i]]$y - mean(site_effects[[i]]$y)
    y_site = site_effects[[i]]$y
    y[idx_sites == i] = y[idx_sites == i] + y_site
  }

  out = list()
  out$splines = list(main = main_effect, site_effects = site_effects)
  out$data = data.frame(y = y, x = x, site = idx_sites)

  return(out)
}

simCategorical = function(n = 1000L, nsites = 3L, ncat, from = NULL, to = NULL, offset = 0) {
  if (is.null(from)) from = runif(1, -100, 100)
  if (is.null(to)) to = from + runif(1, 0, 100)
  means = runif(ncat, min = from, max = to)
  means = means - mean(means) + offset
  idx = sample(seq_len(ncat), n, replace = TRUE)
  x = means[idx]

}




# shared_params contains:
#   - the parameter for the spline, knots, etc. for numerical features
#   - the means for categorical features.
# The site effect is calculated at the server if site_effects and the seed returned to the host to reconstruct all data.

simData = function(n, pnum, ncat, snr = 1, site_effects = TRUE, base_seed = NULL, shared_params) {

}
