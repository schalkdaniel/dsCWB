#' @title Client model
#'
#' @description
#' This class creates a CWB model container for the sites.
#' It collects all base learner and also the data (data frame,
#' pseudo residuals, loss function, etc.) per site.
#'
#' @export
ClientModel = R6Class("ClientModel",
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param symbol (`character(1L)`)\cr
    #'   Character containing the name of the data.
    #' @param target (`character(1L)`)\cr
    #'   Character containing the name of the target variable.
    #' @param feature_names (`character()`)\cr
    #'   Character vector of all target variables.
    #' @param learning_rate (`numeric(1L)`)\cr
    #'   Learning rate.
    #' @param df (`numeric(1L)`)\cr
    #'   Degrees of freedom.
    #' @param nknots (`integer(1L)`)\cr
    #'   Number of inner knots.
    #' @param ord (`integer(1L)`)\cr
    #'   Polynomial degree of basis functions.
    #' @param derivs (`integer(1L)`)\cr
    #'   Number of penalized differences.
    #' @param val_fraction (`numeric(1L)`)\cr
    #'   Fraction of observations used for validation.
    #' @param positive (`integer(1L)`)\cr
    #'   Character indicating the positive class in the binary classification setting.
    #' @param seed (`integer(1L)`)\cr
    #'   Seed used for drawing the validation data.
    initialize = function(symbol = "D", target, feature_names = NULL, learning_rate = 0.1, df = 5L,
      nknots = 20L, ord = 3L, derivs = 2L, val_fraction = NULL, positive = NULL, seed = NULL) {

      checkSymbol(symbol)

      checkmate::assertCharacter(x = symbol, len = 1L, any.missing = FALSE)
      checkmate::assertCharacter(x = target, len = 1L, any.missing = FALSE)
      checkmate::assertCharacter(x = feature_names, any.missing = FALSE, null.ok = TRUE)

      checkmate::assertNumeric(x = learning_rate, len = 1L, any.missing = FALSE)
      checkmate::assertNumeric(x = df, len = 1L, any.missing = FALSE)

      checkmate::assertIntegerish(x = nknots, len = 1L, any.missing = FALSE)
      checkmate::assertIntegerish(x = ord, len = 1L, any.missing = FALSE)
      checkmate::assertIntegerish(x = derivs, len = 1L, any.missing = FALSE)
      checkmate::assertNumeric(x = val_fraction, len = 1L, any.missing = FALSE, null.ok = TRUE)
      checkmate::assertCharacter(x = positive, len = 1L, any.missing = FALSE, null.ok = TRUE)
      checkmate::assertIntegerish(x = seed, len = 1L, any.missing = FALSE, null.ok = TRUE)

      n = eval(parse(text = paste0("nrow(", symbol, ")")), envir = .GlobalEnv)
      if (! is.null(val_fraction)) {
        if (! is.null(seed)) set.seed(seed)
        val_idx = sample(seq_len(n), ceiling(n * val_fraction))
      } else {
        val_idx = NULL
      }

      private$p_n = n
      private$p_symbol = symbol
      private$p_target = target
      private$p_feature_names = feature_names
      private$p_learning_rate = learning_rate
      private$p_df = df
      private$p_nknots = nknots
      private$p_ord = ord
      private$p_derivs = derivs
      private$p_val_idx = val_idx
      private$p_train_idx = setdiff(seq_len(n), val_idx)

      tt = private$getTargetRaw()
      private$p_task = getTaskType(tt)

      if (private$p_task == "bin-class") {
        if (! is.null(positive)) {
          if (! positive %in% unique(tt)) stop("Positive class is not in the target variable!")
          private$p_positive = positive
        } else {
          p_positive = unique(tt)[1]
        }
        private$p_loss = LossBinomial$new()
      }
      if (private$p_task == "regression") {
        private$p_loss = LossQuadratic$new()
      }
      if (is.null(feature_names)) {
        fn = eval(parse(text = paste0("names(", private$p_symbol, ")")), envir = .GlobalEnv)
        private$p_feature_names = setdiff(fn, target)
      }
    },

    #' @description
    #' Get loss optimal constant.
    getLossOptimalConstant = function() {
      return(private$p_loss$init(self$getTarget("train")))
    },

    #' @description
    #' Initialize the loss optimal constant.
    #' @param const (`numeric(1L)`)\cr
    #'    Loss optimal constant.
    initConstantModel = function(const) {
      if (! is.null(private$p_prediction)) stop("Model already initialized")
      checkmate::assertNumeric(x = const, len = 1L, any.missing = FALSE)
      private$p_offset = const
      private$p_prediction = const
      private$p_pseudo_resids = private$p_loss$pseudoResids(self$getTarget(), const)
    },

    #' @description
    #' Get the SSE for each base learner.
    #' @param par_list (`list()`)\cr
    #'   List of parameter values. Each entry must have a name
    #'   that corresponds to a base learner.
    #' @param par_list_binary (`character(1L)`)\cr
    #'   Parameter list as encoded binary object.
    getSSE = function(par_list = NULL, par_list_binary = NULL) {

      checkmate::assertCharacter(par_list_binary, len = 1L, null.ok = TRUE)

      if ((! is.null(par_list)) && (! is.null(par_list_binary)))
        stop("Use either par_list or par_list_binary but not both")

      if (! is.null(par_list_binary))
        par_list = decodeBinary(par_list_binary)

      checkmate::assertList(par_list, null.ok = TRUE)
      nuisance = lapply(names(par_list), function(pln) checkmate::assertChoice(pln, self$getBlNames()))

      out_pl = NA
      if (! is.null(par_list)) {
        out_pl = lapply(names(par_list), function(pln) {
          pred = private$p_bls[[pln]]$predict(par_list[[pln]], "train")
          return(sum((self$getPseudoResids("train") - pred)^2))
        })
        names(out_pl) = names(par_list)
      }

      out_bl = lapply(private$p_bls, function(bl) {
        par = bl$train(self$getPseudoResids("train"))
        pred = bl$predict(par, "train")
        return(sum((self$getPseudoResids("train") - pred)^2))
      })

      return(list(from_list = out_pl, from_bl = out_bl))
    },

    #' @description
    #' Get the risk.
    #' @description
    #' Get target from the data from the symbol.
    #' @param type (`character(1L)`)\cr
    #'   Type which target shold be returned. Choices are
    #'   `full` (for the whole vector), `train` (for the vector used for training),
    #'   and `val` (for the validation vector).
    getRisk = function(type = "train") {
      checkmate::assertChoice(type, c("full", "train", "val"))
      if (is.null(self$getPrediction(type)))
        stop("Prediction must be initialized to calculate risk.")
      return(private$p_loss$risk(self$getTarget(type), self$getPrediction(type)))
    },

    #' @description
    #' Return XtX per base learner.
    getXtX = function() {
      return(lapply(private$p_bls, function(bl) bl$getXtX()))
    },

    #' @description
    #' Get Xty per base learner.
    getXty = function() {
      return(lapply(private$p_bls, function(bl) bl$getXty(self$getPseudoResids("train"), "train")))
    },

    #' @description
    #' Get pseudo residuals.
    #' @param type (`character(1L)`)\cr
    #'   Type which target shold be returned. Choices are
    #'   `full` (for the whole vector), `train` (for the vector used for training),
    #'   and `val` (for the validation vector).
    getPseudoResids = function(type = "full") {
      checkmate::assertChoice(type, c("full", "train", "val"))
      if (type == "full") return(private$p_pseudo_resids)
      if (type == "train") return(private$p_pseudo_resids[private$p_train_idx])
      if (type == "val") return(private$p_pseudo_resids[private$p_val_idx])
    },

    #' @description
    #' Get base learner names.
    getBlNames = function() {
      return(names(private$p_bls))
    },

    #' @description
    #' Get base learner parameters.
    getBlParams = function() {
      return(lapply(private$p_bls, function(bl) bl$getParam()))
    },

    #' @description
    #' Update the model by adding a base learner.
    #' @param blname (`character(1L)`)\cr
    #'   Character containing the base learner name.
    #' @param par (`numeric()`)\cr
    #'   Parameter used to update pseudo residuals and predictions.
    #' @param par_binary (`character(1L)`)\cr
    #'   Binary encoded parameter used to update pseudo residuals and predictions.
    updateCWBParts = function(blname, par = NULL, par_binary = NULL) {
      checkmate::assertCharacter(blname, len = 1L, any.missing = FALSE)
      checkmate::assertNumeric(par, any.missing = FALSE, null.ok = TRUE)
      checkmate::assertCharacter(par_binary, len = 1L, any.missing = FALSE, null.ok = TRUE)

      if (! xor(is.null(par), is.null(par_binary)))
        stop("Just one of par or par_binary must be specified")

      if (! is.null(par_binary)) {
        par = decodeBinary(par_binary)
        checkmate::assertNumeric(par, any.missing = FALSE)
      }

      private$p_prediction = private$p_prediction + private$p_learning_rate * private$p_bls[[blname]]$predict(par)
      private$p_pseudo_resids = private$p_loss$pseudoResids(self$getTarget(), self$getPrediction())

    },

    #' @description
    #' Update the model by adding a base learner.
    #' @param blname (`character(1L)`)\cr
    #'   Character containing the base learner name.
    update = function(blname) {
      checkmate::assertChoice(blname, self$getBlNames())
      par = private$p_bls[[blname]]$train(self$getPseudoResids("train"))
      private$p_bls[[blname]]$updateParam(par * private$p_learning_rate)

      self$updateCWBParts(blname, par)
    },

    #' @description
    #' Re-set the penalty terms.
    #' @param ll_pen (`list()`)\cr
    #'   Named list with new penalty terms.
    #' @param anistrop (`logical(1L)`)\cr
    #'   Flag indicating whether the penalty should be done anistrop or isotrop.
    updatePenalty = function(ll_pen, anistrop = TRUE) {
      if (! anistrop) {
        checkmate::assertList(ll_pen)

        blnames = names(ll_pen)
        blnames0 = names(private$p_bls)

        nuisance = lapply(blnames, checkmate::assertChoice, choices = blnames0)

        for (bln in blnames) {
          private$p_bls[[bln]]$updatePenalty(ll_pen[[bln]])
        }
      } else {
        checkmate::assertNumeric(ll_pen)

        blp = private$bls[[bln]]$getPenalty()
        blpmat = private$bls[[bln]]$getPenaltyMat()
        pnew = blp * blpmat + ll_pen * diag(ncol(blpmat))
        private$bls[[bln]]$updatePenalty(pnew)
      }
    },

    #' @description
    #' Add base learners to the model.
    #' @param ll_init (`list()`)\cr
    #'   List containing the initialization parameters.
    addBaselearners = function(ll_init) {
      checkmate::assertList(ll_init)
      for (ff in private$p_feature_names) {
        x = private$getFeatureFromData(ff)
        if (is.numeric(x)) {
          e = try({
            private$addBlNumeric(ff, paste0(ff, "-spline"), ll_init[[ff]]$min, ll_init[[ff]]$max)
          }, silent = TRUE)
          control = "numeric"
        }
        if (is.character(x) || is.factor(x)) {
          e = try({
            private$addBlCategorical(ff, paste0(ff, "-onehot"), ll_init[[ff]]$table)
          }, silent = TRUE)
          control = "categorical"
        }
        if (inherits(e, "try-error")) {
          stop("Error in base learner creation of the ", control, " feature ", ff,
            ": \n ", attr(e, "condition")$message)
        }
      }
    },

    #' @description
    #' Get feature names
    getFeatureNames = function() {
      return(private$p_feature_names)
    },

    #' @description
    #' Get prediction.
    #' @param type (`character(1L)`)\cr
    #'   Type which target shold be returned. Choices are
    #'   `full` (for the whole vector), `train` (for the vector used for training),
    #'   and `val` (for the validation vector).
    getPrediction = function(type = "full") {
      checkmate::assertChoice(type, c("full", "train", "val"))
      if (length(private$p_prediction) == 1L)
        pred = rep(private$p_prediction, private$p_n)
      else
        pred = private$p_prediction

      if (type == "full") return(pred)
      if (type == "train") return(pred[private$p_train_idx])
      if (type == "val") return(pred[private$p_val_idx])
    },

    #' @description
    #' Return offset of the model.
    getOffset = function() {
      return(private$p_offset)
    },

    #' @description
    #' Get target from the data from the symbol.
    #' @param type (`character(1L)`)\cr
    #'   Type which target shold be returned. Choices are
    #'   `full` (for the whole vector), `train` (for the vector used for training),
    #'   and `val` (for the validation vector).
    getTarget = function(type = "full") {
      if (is.null(private$p_task)) stop("No task is set!")

      if (private$p_task == "regression") return(private$getTargetRaw(type))

      if (private$p_task == "bin-class")
        return(ifelse(private$getTargetRaw(type) == private$p_positive, 1, -1))
    },

    #' @description
    #' Get the number of rows per train and validation set.
    getTrainValObs = function() {
      return(list(
        ntrain = length(private$p_train_idx),
        nval = length(private$p_val_idx)))
    }
  ),
  active = list(
    #' @field loss (`Loss`)\cr
    #'   Loss object used for modelling.
    loss = function(x) {
      if (! missing(x)) stop("`loss` is read only.")
      return(private$p_loss)
    },

    #' @field bls (`list(Baselearner)`)\cr
    #'   List of base learner.
    bls = function(x) {
      if (! missing(x)) stop("`bls` is read only.")
      return(private$p_bls)
    }
  ),
  private = list(
    # @description
    # Get feature from the data from the symbol.
    #
    # @param feature (`character(1L)`)\cr
    #   Character containing the feature name.
    getFeatureFromData = function(feature) {
      checkmate::assertCharacter(x = feature, len = 1L, any.missing = FALSE)
      x = eval(parse(text = paste0(private$p_symbol, "[[\"", feature, "\"]]")), envir = .GlobalEnv)
      if (is.null(x)) stop("Feature \"", feature, "\" was not found in ", symbol)
      return(x)
    },

    # @description
    # Get target from the data from the symbol.
    # @param type (`character(1L)`)\cr
    #   Type which target shold be returned. Choices are
    #   `full` (for the whole vector), `train` (for the vector used for training),
    #   and `val` (for the validation vector).
    getTargetRaw = function(type = "full") {
      checkmate::assertChoice(type, c("full", "train", "val"))
      if (type == "full") return(private$getFeatureFromData(private$p_target))
      if (type == "train") return(private$getFeatureFromData(private$p_target)[private$p_train_idx])
      if (type == "val") return(private$getFeatureFromData(private$p_target)[private$p_val_idx])
    },


    # @description
    # Add a numeric base learner.
    #
    # @param feature (`character(1L)`)\cr
    #   Character containing the feature name.
    # @param blname (`character(1L)`)\cr
    #   Base learner name (id used for references).
    # @param knots_min (`numeric(1L)`)\cr
    #   Minimum value of inner knots.
    # @param knots_max (`numeric(1L)`)\cr
    #   Maximum value of inner knots.
    addBlNumeric = function(feature, blname, knots_min, knots_max) {
      checkmate::assertCharacter(x = feature, len = 1L, any.missing = FALSE)
      checkmate::assertNumeric(x = knots_min, len = 1L, null.ok = FALSE)
      checkmate::assertNumeric(x = knots_max, len = 1L, null.ok = FALSE)

      x = private$getFeatureFromData(feature)
      bl = BlSpline$new(knots_min, knots_max, private$p_nknots, private$p_ord, private$p_derivs)
      bl$initData(x, feature, private$p_df, private$p_val_idx)

      private$p_bls[[blname]] = bl
    },

    # @description
    # Add a categorical base learner.
    #
    # @param feature (`character(1L)`)\cr
    #   Character containing the feature name.
    # @param blname (`character(1L)`)\cr
    #   Base learner name (id used for references).
    # @param dict (`character()`)\cr
    #   Unique categories of the feature.
    addBlCategorical = function(feature, blname, dict) {
      checkmate::assertCharacter(x = feature, len = 1L, any.missing = FALSE)
      checkmate::assertCharacter(x = blname, len = 1L, any.missing = FALSE)
      checkmate::assertCharacter(x = dict, any.missing = FALSE)

      x = as.character(private$getFeatureFromData(feature))
      bl = BlOneHot$new(dict)
      if (private$p_df > length(dict))
        df = length(dict)
      else
        df = private$p_df

      bl$initData(x, feature, df, private$p_val_idx)

      private$p_bls[[blname]] = bl
    },

    p_n = NULL,
    p_symbol = NULL,
    p_target = NULL,
    p_feature_names = NULL,
    p_learning_rate = NULL,
    p_df = NULL,
    p_nknots = NULL,
    p_ord = NULL,
    p_derivs = NULL,
    p_train_idx = NULL,
    p_val_idx = NULL,
    p_bls = NULL,
    p_pseudo_resids = NULL,
    p_prediction = NULL,
    p_offset = NULL,
    p_task = NULL,
    p_positive = NULL,
    p_loss = NULL
  )
)
