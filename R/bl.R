#' @title P-Spline base learner class
#'
#' @description
#' This class creates a base learner that can be used for calculating
#' predictions, fitting parameters, etc.
#'
#' @export
BlSpline = R6Class("BlSpline",
  public = list(

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param knots_min (`numeric(1L)`)\cr
    #'   Minimum value of inner knots.
    #' @param knots_max (`numeric(1L)`)\cr
    #'   Maximum value of inner knots.
    #' @param nknots (`integer(1L)`)\cr
    #'   Number of inner knots.
    #' @param ord (`integer(1L)`)\cr
    #'   Polynomial degree of the basis functions.
    #' @param derivs (`integer(1L)`)\cr
    #'   Number of penalized differences.
    initialize = function(knots_min, knots_max,
      nknots = 20L, ord = 3L, derivs = 2L) {

      checkmate::assertNumeric(x = knots_min, len = 1L, null.ok = FALSE)
      checkmate::assertNumeric(x = knots_max, len = 1L, null.ok = FALSE)
      checkmate::assertIntegerish(x = nknots, len = 1L, min = 1)
      checkmate::assertIntegerish(x = ord, len = 1L, min = 1)
      checkmate::assertIntegerish(x = derivs, len = 1L, min = 1)

      private$p_knots_min = knots_min
      private$p_knots_max = knots_max
      private$p_ord    = ord
      private$p_nknots = nknots
      private$p_derivs = derivs
    },

    #' @description
    #' Calculates the basis transformation of the base learner.
    #' @param x (`numeric()`)\cr
    #'   Data vector.
    basisTrafo = function(x) {
      checkmate::assertNumeric(x, any.missing = FALSE)
      return(trafoSpline(x, private$p_knots_min, private$p_knots_max,
        private$p_nknots, private$p_ord))
    },

    #' @description
    #' Initializes data based on given data. This function should
    #' be used at the sites.
    #'
    #' @param x (`numeric()`)\cr
    #'   Data vector.
    #' @param feature (`character(1L)`)\cr
    #'   Character containing the name of the feature.
    #' @param df (`numeric(1L)`)\cr
    #'   Value indicating the degrees of freedom for the base learner.
    #' @param val_idx (`integer()`)\cr
    #'   Vector containing the indices used for validation.
    #' @param pen (`numeric(1L)`)\cr
    #'   Penalty term for the base learner. If set, df is ignored.
    initData = function(x, feature, df, val_idx = NULL, pen = NULL) {
      checkmate::assertNumeric(x = x, any.missing = FALSE)
      checkmate::assertCharacter(x = feature, len = 1L)
      checkmate::assertNumeric(x = df, any.missing = FALSE, len = 1L)

      #x = eval(parse(text = paste0(symbol, "[[\"", feature, "\"]]")))
      #if (is.null(x)) stop("Feature \"", feature, "\" was not found in ", symbol)

      private$p_feature = feature
      private$p_df = df
      private$p_val_idx = val_idx

      train_idx = setdiff(seq_along(x), val_idx)
      private$p_train_idx = train_idx
      checkmate::assertIntegerish(x = val_idx, min = 1, max = length(x), null.ok = TRUE)
      private$p_design = trafoSpline(x, private$p_knots_min, private$p_knots_max,
        private$p_nknots, private$p_ord)

      private$p_xtx     = t(private$p_design[train_idx, ]) %*% private$p_design[train_idx, ]
      private$p_penmat  = compboostSplines::penaltyMat(ncol(private$p_design), private$p_derivs)

      if (is.null(pen))
        private$p_penalty = compboostSplines::demmlerReinsch(as.matrix(private$p_xtx), private$p_penmat, private$p_df)
      else
        private$p_penalty = pen

      private$p_xtx_inv = solve(private$p_xtx + private$p_penalty * private$p_penmat)
    },

    #' @description
    #' Initializes data based on the hat matrix. This function should
    #' be used at the host.
    #'
    #' @param feature (`character(1L)`)\cr
    #'   Character containing the name of the feature.
    #' @param df (`numeric(1L)`)\cr
    #'   Value indicating the degrees of freedom for the base learner.
    #' @param xtx (`matrix()`)\cr
    #'   Hat matrix used for parameter estimation.
    initDataXtX = function(feature, df, xtx) {
      checkmate::assertCharacter(x = feature, len = 1L)

      private$p_feature = feature
      private$p_df      = df

      private$p_xtx     = xtx
      private$p_penmat  = compboostSplines::penaltyMat(ncol(xtx), private$p_derivs)
      private$p_penalty = compboostSplines::demmlerReinsch(private$p_xtx, private$p_penmat, private$p_df)
      private$p_xtx_inv = solve(private$p_xtx + private$p_penalty * private$p_penmat)
    },

    #' @description
    #' Initializes an empty base learner. This one can be used for
    #' simple predictions without data in it.
    #'
    #' @param feature (`character(1L)`)\cr
    #'   Character containing the name of the feature.
    #' @param df (`numeric(1L)`)\cr
    #'   Value indicating the degrees of freedom for the base learner.
    #' @param penalty (`numeric(1L)`)\cr
    #'   Set the penalty (just for information, has no effect).
    initEmpty = function(feature, df = NULL, penalty = NULL) {
      checkmate::assertCharacter(x = feature, len = 1L)
      checkmate::assertNumeric(x = df, len = 1L, any.missing = FALSE, null.ok = TRUE)
      checkmate::assertNumeric(x = penalty, len = 1L, any.missing = FALSE, null.ok = TRUE)

      private$p_feature = feature
      private$p_df      = df
      private$p_penalty = penalty
    },

    #' @description
    #' Get XtX.
    getXtX = function() {
      return(private$p_xtx)
    },

    #' @description
    #' Get penalty matrix.
    getPenaltyMat = function() {
      return(private$p_penmat)
    },

    #' @description
    #' Get Xty.
    #' @param y (`numeric()`)\cr
    #'   y variable which is multiplied with X
    #' @param type (`character(1L)`)\cr
    #'   Type which target should be returned. Choices are
    #'   `full` (for the whole vector), `train` (for the vector used for training),
    #'   and `val` (for the validation vector).
    getXty = function(y, type = "full") {
      checkmate::assertChoice(type, c("full", "train", "val"))

      if (type == "full") xt = t(private$p_design)
      if (type == "train") xt = t(private$p_design[private$p_train_idx, ])
      if (type == "val") xt = t(private$p_design[private$p_val_idx, ])

      checkmate::assertNumeric(y, len = ncol(xt), any.missing = FALSE)
      return(as.vector(xt %*% y))
    },


    #' @description
    #' Calculate parameter estimates based on
    #' be used at the host.
    #'
    #' @param y (`numeric(1L)`)\cr
    #'   Response vector.
    #' @param xty (`matrix()|numeric()`)\cr
    #'   X^t y with X design matrix and y response vector.
    train = function(y = NULL, xty = NULL) {
      checkmate::assertNumeric(x = y, len = length(private$p_train_idx), any.missing = FALSE, null.ok = TRUE)
      checkmate::assertNumeric(x = xty, len = ncol(private$p_xty_inv), any.missing = FALSE, null.ok = TRUE)

      if (! xor(is.null(y), is.null(xty)))
        stop("Either `y` or `xty`, but not both, must be given")

      if (! is.null(xty)) {
        out = private$p_xtx_inv %*% xty
      } else {
        out = private$p_xtx_inv %*% t(private$p_design[private$p_train_idx, ]) %*% y
      }
      return(as.vector(out))
    },

    #' @description
    #' Calculate prediction based on a parameter vector.
    #'
    #' @param par (`numeric()`)\cr
    #'   Parameter vector used for predicting on the data set
    #'   in design. If `NULL`, the internal parameter vector is used.
    #' @param type (`character(1L)`)\cr
    #'   Type which predictions should be returned. Choices are
    #'   `full` (for the whole vector), `train` (for the vector used for training),
    #'   and `val` (for the validation vector).
    predict = function(par = NULL, type = "full") {
      checkmate::assertNumeric(par, len = ncol(private$p_xtx), null.ok = TRUE)

      if (is.null(private$p_design)) return(NA)
      if (is.null(par)) par = private$p_param

      if (type == "full") out = private$p_design %*% par
      if (type == "train") out = private$p_design[private$p_train_idx, ] %*% par
      if (type == "val") out = private$p_design[private$p_val_idx, ] %*% par

      return(as.vector(out))
    },

    #' @description
    #' Calculate prediction for the new data.
    #'
    #' @param newdata (`data.frame()`)\cr
    #'   New data. Note that `newdata` must contain the feature.
    #' @param par (`numeric()`)\cr
    #'   Parameter vector used for predicting on new data.
    predictNewdata = function(newdata, par = NULL) {
      checkmate::assertDataFrame(newdata)
      checkmate::assertNumeric(par, len = ncol(private$p_xtx), null.ok = TRUE)

      x = newdata[[private$p_feature]]
      if (is.null(x)) stop("Feature \"", private$p_feature, "\" was not found in ", newdata)

      design = trafoSpline(x, private$p_knots_min, private$p_knots_max,
        private$p_nknots, private$p_ord)

      if (is.null(par)) {
        if (is.null(private$p_param))
          out = rep(0, nrow(design))
        else
          out = design %*% private$p_param
      } else {
        out = design %*% par
      }
      return(as.vector(out))
    },

    #' @description
    #' Add parameter vector to the already stored parameter vector.
    #'
    #' @param par (`numeric()`)\cr
    #'   Parameter vector added to the internal stored parameter.
    updateParam = function(par) {
      checkmate::assertNumeric(par, len = ncol(private$p_xtx))

      if (is.null(private$p_param)) private$p_param = 0
      private$p_param = private$p_param + par
    },

    #' @description
    #' Get knot range.
    getKnotRange = function() {
      return(c(private$p_knots_min, private$p_knots_max))
    },

    #' @description
    #' Re-set the penalty terms.
    #' @param pen (`numeric(1)`)\cr
    #'   New penalty terms.
    updatePenalty = function(pen) {
      checkmate::assertNumeric(pen, len = 1L, lower = 0, any.missing = FALSE)
      private$p_penalty = pen
    },

    #' @description
    #' Get feature name.
    getFeature = function() {
      return(private$p_feature)
    },

    #' @description
    #' Get hyperparameters (df, penalty, order, and derivs).
    getHyperpars = function() {
      return(list(df = private$p_df, penalty = private$p_penalty,
        ord = private$p_ord, derivs = private$p_derivs))
    },

    #' @description
    #' Get parameter vector.
    getParam = function() {
      return(private$p_param)
    },

    #' @description
    #' Get the base learner type.
    getType = function() {
      return("numeric")
    }
  ),
  private = list(
    p_knots_min = NULL,
    p_knots_max = NULL,
    p_nknots = NULL,
    p_ord = NULL,
    p_derivs = NULL,
    p_feature = NULL,
    p_df = NULL,
    p_train_idx = NULL,
    p_val_idx = NULL,
    p_design = NULL,
    p_xtx = NULL,
    p_xtx_inv = NULL,
    p_penmat = NULL,
    p_penalty = NULL,
    p_param = NULL
  )
)

#' @title One-hot encoding base learner class
#'
#' @description
#' This class creates a base learner that can be used for calculating
#' predictions, fitting parameters, etc.
#'
#' @export
BlOneHot = R6Class("BlOneHot",
  public = list(

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param categories (`character()`)\cr
    #'   Categories of the categorical feature.
    initialize = function(categories) {

      checkmate::assertCharacter(x = categories, any.missing = FALSE)

      dict = seq_along(categories)
      names(dict) = categories

      private$p_dict = dict
    },

    #' @description
    #' Calculates the basis transformation of the base learner.
    #' @param x (`numeric()`)\cr
    #'   Data vector.
    basisTrafo = function(x) {
      if (is.factor(x)) x = as.character(x)
      checkmate::assertCharacter(x, any.missing = FALSE)
      ind = private$p_dict[x]

      i = seq_along(ind)[! is.na(ind)]
      j = na.omit(ind)

      return(Matrix::sparseMatrix(i = i, j = j, x = 1,
        dims = c(length(ind), length(private$p_dict))))
    },

    #' @description
    #' Initializes data based on given data. This function should
    #' be used at the sites.
    #'
    #' @param x (`numeric()`)\cr
    #'   Data vector.
    #' @param feature (`character(1L)`)\cr
    #'   Character containing the name of the feature.
    #' @param df (`numeric(1L)`)\cr
    #'   Value indicating the degrees of freedom for the base learner.
    #' @param val_idx (`integer()`)\cr
    #'   Vector containing the indices used for validation.
    #' @param pen (`numeric(1L)`)\cr
    #'   Penalty term for the base learner. If set, df is ignored.
    initData = function(x, feature, df, val_idx = NULL, pen = NULL) {
      checkmate::assertCharacter(x = x, any.missing = FALSE)
      checkmate::assertCharacter(x = feature, len = 1L)
      checkmate::assertNumeric(x = df, any.missing = FALSE, len = 1L)

      private$p_feature = feature
      private$p_df = df
      private$p_val_idx = val_idx

      train_idx = setdiff(seq_along(x), val_idx)
      private$p_train_idx = train_idx
      checkmate::assertIntegerish(x = val_idx, min = 1, max = length(x), null.ok = TRUE)
      private$p_design = self$basisTrafo(x)

      private$p_xtx     = t(private$p_design[train_idx, ]) %*% private$p_design[train_idx, ]
      private$p_penmat  = diag(ncol(private$p_design))

      if (is.null(pen))
        private$p_penalty = compboostSplines::demmlerReinsch(as.matrix(private$p_xtx), private$p_penmat, private$p_df)
      else
        private$p_penalty = pen

      private$p_xtx_inv = diag(1 / diag(private$p_xtx + private$p_penalty * private$p_penmat))
    },

    #' @description
    #' Initializes data based on the hat matrix. This function should
    #' be used at the host.
    #'
    #' @param feature (`character(1L)`)\cr
    #'   Character containing the name of the feature.
    #' @param df (`numeric(1L)`)\cr
    #'   Value indicating the degrees of freedom for the base learner.
    #' @param xtx (`matrix()`)\cr
    #'   Hat matrix used for parameter estimation.
    initDataXtX = function(feature, df, xtx) {
      checkmate::assertCharacter(x = feature, len = 1L)

      private$p_feature = feature
      private$p_df      = df

      private$p_xtx     = xtx
      private$p_penmat  = diag(ncol(xtx))
      private$p_penalty = compboostSplines::demmlerReinsch(private$p_xtx, private$p_penmat, private$p_df)
      private$p_xtx_inv = diag(1 / diag(private$p_xtx + private$p_penalty * private$p_penmat))
    },

    #' @description
    #' Initializes an empty base learner. This one can be used for
    #' simple predictions without data in it.
    #'
    #' @param feature (`character(1L)`)\cr
    #'   Character containing the name of the feature.
    #' @param df (`numeric(1L)`)\cr
    #'   Value indicating the degrees of freedom for the base learner.
    #' @param penalty (`numeric(1L)`)\cr
    #'   Set the penalty (just for information, has no effect).
    initEmpty = function(feature, df = NULL, penalty = NULL) {
      checkmate::assertCharacter(x = feature, len = 1L)
      checkmate::assertNumeric(x = df, len = 1L, any.missing = FALSE, null.ok = TRUE)
      checkmate::assertNumeric(x = penalty, len = 1L, any.missing = FALSE, null.ok = TRUE)

      private$p_feature = feature
      private$p_df      = df
      private$p_penalty = penalty
    },

    #' @description
    #' Get XtX.
    getXtX = function() {
      return(private$p_xtx)
    },

    #' @description
    #' Get penalty matrix.
    getPenaltyMat = function() {
      return(private$p_penmat)
    },

    #' @description
    #' Get Xty.
    #' @param y (`numeric()`)\cr
    #'   y variable which is multiplied with X
    #' @param type (`character(1L)`)\cr
    #'   Type which target should be returned. Choices are
    #'   `full` (for the whole vector), `train` (for the vector used for training),
    #'   and `val` (for the validation vector).
    getXty = function(y, type = "full") {
      checkmate::assertChoice(type, c("full", "train", "val"))

      if (type == "full") xt = t(private$p_design)
      if (type == "train") xt = t(private$p_design[private$p_train_idx, ])
      if (type == "val") xt = t(private$p_design[private$p_val_idx, ])

      checkmate::assertNumeric(y, len = ncol(xt), any.missing = FALSE)
      return(as.vector(xt %*% y))
    },


    #' @description
    #' Calculate parameter estimates based on
    #' be used at the host.
    #'
    #' @param y (`numeric(1L)`)\cr
    #'   Response vector.
    #' @param xty (`matrix()|numeric()`)\cr
    #'   X^t y with X design matrix and y response vector.
    train = function(y = NULL, xty = NULL) {
      checkmate::assertNumeric(x = y, len = length(private$p_train_idx), any.missing = FALSE, null.ok = TRUE)
      checkmate::assertNumeric(x = xty, len = ncol(private$p_xty_inv), any.missing = FALSE, null.ok = TRUE)

      if (! xor(is.null(y), is.null(xty)))
        stop("Either `y` or `xty`, but not both, must be given")

      if (! is.null(xty)) {
        out = private$p_xtx_inv %*% xty
      } else {
        out = private$p_xtx_inv %*% t(private$p_design[private$p_train_idx, ]) %*% y
      }
      return(as.vector(out))
    },

    #' @description
    #' Calculate prediction based on a parameter vector.
    #'
    #' @param par (`numeric()`)\cr
    #'   Parameter vector used for predicting on the data set
    #'   in design. If `NULL`, the internal parameter vector is used.
    #' @param type (`character(1L)`)\cr
    #'   Type which predictions should be returned. Choices are
    #'   `full` (for the whole vector), `train` (for the vector used for training),
    #'   and `val` (for the validation vector).
    predict = function(par = NULL, type = "full") {
      checkmate::assertNumeric(par, len = ncol(private$p_xtx), null.ok = TRUE)

      if (is.null(private$p_design)) return(NA)
      if (is.null(par)) par = private$p_param

      if (type == "full") out = private$p_design %*% par
      if (type == "train") out = private$p_design[private$p_train_idx, ] %*% par
      if (type == "val") out = private$p_design[private$p_val_idx, ] %*% par

      return(as.vector(out))
    },

    #' @description
    #' Calculate prediction for the new data.
    #'
    #' @param newdata (`data.frame()`)\cr
    #'   New data. Note that `newdata` must contain the feature.
    #' @param par (`numeric()`)\cr
    #'   Parameter vector used for predicting on new data.
    predictNewdata = function(newdata, par = NULL) {
      checkmate::assertDataFrame(newdata)
      checkmate::assertNumeric(par, len = ncol(private$p_xtx), null.ok = TRUE)

      x = newdata[[private$p_feature]]
      if (is.null(x)) stop("Feature \"", private$p_feature, "\" was not found in ", newdata)

      design = self$basisTrafo(x)

      if (is.null(par)) {
        if (is.null(private$p_param))
          out = rep(0, nrow(design))
        else
          out = design %*% private$p_param
      } else {
        out = design %*% par
      }
      return(as.vector(out))
    },

    #' @description
    #' Add parameter vector to the already stored parameter vector.
    #'
    #' @param par (`numeric()`)\cr
    #'   Parameter vector added to the internal stored parameter.
    updateParam = function(par) {
      checkmate::assertNumeric(par, len = ncol(private$p_xtx))

      if (is.null(private$p_param)) private$p_param = 0
      private$p_param = private$p_param + par
    },

    #' @description
    #' Get knot range.
    getDictionary = function() {
      return(private$p_dict)
    },

    #' @description
    #' Re-set the penalty terms.
    #' @param pen (`numeric(1)`)\cr
    #'   New penalty terms.
    updatePenalty = function(pen) {
      checkmate::assertNumeric(pen, len = 1L, lower = 0, any.missing = FALSE)
      private$p_penalty = pen
    },

    #' @description
    #' Get feature name.
    getFeature = function() {
      return(private$p_feature)
    },

    #' @description
    #' Get hyperparameters (df, penalty, order, and derivs).
    getHyperpars = function() {
      return(list(df = private$p_df, penalty = private$p_penalty))
    },

    #' @description
    #' Get parameter vector.
    getParam = function() {
      return(private$p_param)
    },

    #' @description
    #' Get the base learner type.
    getType = function() {
      return("categorical")
    }
  ),
  private = list(
    p_feature = NULL,
    p_dict = NULL,
    p_df = NULL,
    p_train_idx = NULL,
    p_val_idx = NULL,
    p_design = NULL,
    p_xtx = NULL,
    p_xtx_inv = NULL,
    p_penmat = NULL,
    p_penalty = NULL,
    p_param = NULL
  )
)

if (FALSE) {
  devtools::document()

  nsim = 100L

  x = runif(nsim, 0, 10)
  y = sin(x) + rnorm(nsim, 0, 0.2)

  bls = BlSpline$new(min(x), max(x))
  bls$initData(x, "x", 10)
  par = bls$train(y)
  pred = bls$predict(par)

  plot(x, y)

  ord = order(x)
  lines(x = x[ord], y = pred[ord], col = "red")

  bls$updateParam(par)
  bls$getParam()
  bls$predict()
  bls$getHyperpars()



  bls = BlOneHot$new(LETTERS[1:4])
  bls$getDictionary()
  x = LETTERS[1:4][sample(4, 10, TRUE)]
  x[4] = "bla"
  bls$basisTrafo(x)
  bls$initData(x, "x", df = 2)
  bls$getXtX()
  bls$getHyperpars()
  bls$getXty(1:10)
  param = bls$train(1:10)
  bls$getParam()
  bls$updateParam(param)
  bls$getParam()

  newdat = rev(x)
  bls$predictNewdata(data.frame(x = newdat))
}
