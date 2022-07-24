#' @title Create a B-spline basis representation
#'
#' @param x (`numeric()`) Data vector.
#' @param knots_min (`numeric(1L`) Minimal value for the knots.
#' @param knots_max (`numeric(1L`) Maximal value for the knots.
#' @param nknots (`integer(1L)`) Number of knots.
#' @param ord (`integer(1L)`) Basis function degree.
#' @return Sparse matrix containing the spline basis.
#' @author Daniel S.
#' @import cpsp
#' @export
trafoSpline = function(x, knots_min, knots_max, nknots = 20L, ord = 3L) {
  checkmate::assertNumeric(x = knots_min, len = 1L, any.missing = FALSE)
  checkmate::assertNumeric(x = knots_max, len = 1L, any.missing = FALSE)
  checkmate::assertNumeric(x = x, lower = knots_min, upper = knots_max, any.missing = FALSE)
  checkmate::assertCount(x = nknots, positive = TRUE)
  checkmate::assertCount(x = ord, positive = TRUE)

  if (any(x > knots_max)) {
    idx_max = x > knots_max
    x[idx_max] = knots_max
  }
  if (any(x < knots_min)) {
    idx_min = x < knots_min
    x[idx_min] = knots_min
  }

  knots = cpsp::createKnots(values = c(knots_min, knots_max), n_knots =  nknots, degree =  ord)
  spb = cpsp::createSparseSplineBasis(values = x, degree = ord, knots = knots)
  return(spb)
}

#' @title Create dummy data matrix
#'
#' @param x (`character()`) Data vector.
#' @param levels (`character()`) Levels of the character vector.
#' @return Sparse dummy matrix..
#' @author Daniel S.
#' @export
trafoCategorical = function(x, levels) {
  checkmate::assertCharacter(x = x, any.missing = FALSE)
  checkmate::assertCharacter(x = levels, any.missing = FALSE)

  n = length(x)
  p = length(levels)

  cidx = vapply(x, function(k) which(k == levels), integer(1L))
  X = Matrix::sparseMatrix(i = seq_along(x), j = cidx, x = 1)
  return(X)
}
