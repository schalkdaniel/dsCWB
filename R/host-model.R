#' @title Host model
#'
#' @description
#' This class creates a CWB model container for the host.
#' It collects all base learners.
#'
#' @export
HostModel = R6Class("HostModel",
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param symbol (`character(1L)`)\cr
    #'   Character containing the name of the data.
    #' @param target (`character(1L)`)\cr
    #'   Character containing the name of the target variable.
    #' @param target_type (`character(1L)`)\cr
    #'   The target type, must be `regression` or `bin-class`.
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
    #' @param positive (`integer(1L)`)\cr
    #'   Character indicating the positive class in the binary classification setting.
    initialize = function(symbol = "D", target, target_type = NULL, feature_names = NULL, learning_rate = 0.1, df = 5L,
      nknots = 20L, ord = 3L, derivs = 2L, positive = NULL) {

      checkmate::assertCharacter(x = symbol, len = 1L, any.missing = FALSE)
      checkmate::assertCharacter(x = target, len = 1L, any.missing = FALSE)
      checkmate::assertChoice(x = target_type, choices = c("regression", "bin-class"))
      checkmate::assertCharacter(x = feature_names, any.missing = FALSE)

      checkmate::assertNumeric(x = learning_rate, len = 1L, any.missing = FALSE)
      checkmate::assertNumeric(x = df, len = 1L, any.missing = FALSE)

      checkmate::assertIntegerish(x = nknots, len = 1L, any.missing = FALSE)
      checkmate::assertIntegerish(x = ord, len = 1L, any.missing = FALSE)
      checkmate::assertIntegerish(x = derivs, len = 1L, any.missing = FALSE)
      checkmate::assertCharacter(x = positive, len = 1L, any.missing = FALSE, null.ok = TRUE)

      private$p_symbol = symbol
      private$p_target = target
      private$p_feature_names = feature_names
      private$p_learning_rate = learning_rate
      private$p_df = df
      private$p_nknots = nknots
      private$p_ord = ord
      private$p_derivs = derivs
      private$p_positive = positive

      if (target_type == "regression") private$p_loss = LossQuadratic$new()
    },

    #' @description
    #' Add base learners to the model.
    #' @param ll_init (`list()`)\cr
    #'   List containing the initialization parameters.
    #' @param ll_xtx (`list()`)\cr
    #'   List containing the XtX matrices.
    addBaselearners = function(ll_init, ll_xtx) {
      checkmate::assertList(ll_init)
      checkmate::assertList(ll_xtx)

      private$p_vinit = ll_init

      ninit = names(ll_init)
      nxtx  = names(ll_xtx)

      check = vapply(X = nxtx, FUN = function(nx) {
        any(vapply(ninit, function(ni) grepl(ni, nx), logical(1L)))
      }, FUN.VALUE = logical(1L))
      if (! all(check))
        stop("Names of `ll_xtx` must be names of `ll_init`")

      for (ff in ninit) {
        checkmate::assertChoice(ff, private$p_feature_names)
        if (ll_init[[ff]]$class == "numeric") {
          blname = paste0(ff, "-spline")

          bl = BlSpline$new(ll_init[[ff]]$min, ll_init[[ff]]$max, private$p_nknots, private$p_ord, private$p_derivs)
          bl$initDataXtX(ff, private$p_df, as.matrix(ll_xtx[[blname]]))

          private$p_bls[[blname]] = bl
        }
      }
    },

    #' @description
    #'   Get the data symbol
    getDataSymbol = function() {
      return(private$p_symbol)
    },

    #' @description
    #'   Get feature names
    getFeatureNames = function() {
      return(private$p_feature_names)
    },

    #' @description
    #'   Get parameter based on xty.
    #' @param ll_xty (`list()`)\cr
    #'   List containing the Xty matrices.
    getParam = function(ll_xty) {
      checkmate::assertList(x = ll_xty)
      nuisance = lapply(names(ll_xty), function(name_xty) checkmate::assertChoice(name_xty, names(private$p_bls)))
      ll_pars = list()
      for (fn in names(ll_xty)) {
        ll_pars[[fn]] = private$p_bls[[fn]]$train(xty = ll_xty[[fn]])
      }
      return(ll_pars)
    },

    #' @description
    #' Conduct an update step.
    #' @param blname (`character(1L)`)\cr
    #'   Base learner name.
    #' @param par (`numeric()`)\cr
    #'   Parameter vector corresponding to blname.
    update = function(blname, par) {
      checkmate::assertCharacter(blname, len = 1L, any.missing = FALSE)
      private$p_bls[[blname]]$updateParam(par * private$p_learning_rate)
      #if (grepl("site-", blname)) {
        #datashield.assign(private$p_connections, "cm", paste0("updateClientBaselearner(\"cm\", ", blname, ")"))
      #} else {
        #private$p_bls[[blname]]$updateParam(par * private$p_learning_rate)
        #ds_cl = paste0("updateClientBaselearner(\"cm\", ", blname, ", ", decodeBinary(par), ")")
        #datashield.assign(private$p_connections, "cm", ds_cl)
      #}
    },

    #' @description
    #' Log an updated step.
    #' @param blname (`character(1L)`)\cr
    #'   Base learner name.
    #' @param effect_type (`character(1L)`)\cr
    #'   Type (shared `shared` or site-specific `site`) of the base learner.
    #' @param sse (`numeric(1L)`)\cr
    #'   SSE of the added base learner.
    #' @param risk_train (`numeric(1L)`)\cr
    #'   Training risk.
    #' @param risk_val (`numeric(1L)`)\cr
    #'   Validation risk.
    log = function(blname, effect_type, sse, risk_train = NA, risk_val = NA) {
      checkmate::assertCharacter(blname, len = 1L, any.missing = FALSE)
      checkmate::assertChoice(effect_type, c("shared", "site"))
      checkmate::assertNumeric(sse, len = 1L, any.missing = FALSE)
      checkmate::assertNumeric(risk_train, len = 1L)
      checkmate::assertNumeric(risk_val, len = 1L)

      if (is.null(private$p_log)) {
        private$p_log = data.frame(iteration = 1, bl = blname, effect_type = effect_type, sse = sse,
          risk_train = risk_train, risk_val = risk_val)
      } else {
        lnew = data.frame(iteration = nrow(private$p_log) + 1, bl = blname, effect_type = effect_type, sse = sse,
          risk_train = risk_train, risk_val = risk_val)
        private$p_log = rbind(private$p_log, lnew)
      }
    }
  ),
  active = list(
    #' @field loss (`Loss`)\cr
    #'   Loss object used for modelling.
    loss = function(x) {
      if (! missing(x)) stop("loss is read only.")
      return(private$p_loss)
    },

    #' @field lr (numeric(1L)\cr
    #'   Learning rate.
    lr = function(x) {
      if (! missing(x)) stop("`lr` is read only.")
      return(private$p_learning_rate)
    },

    #' @field bls (`list(Baselearner)`)\cr
    #'   List of base learner.
    bls = function(x) {
      if (! missing(x)) stop("`bls` is read only.")
      return(private$p_bls)
    }
  ),
  private = list(
    p_symbol = NULL,
    p_target = NULL,
    p_feature_names = NULL,
    p_learning_rate = NULL,
    p_df = NULL,
    p_nknots = NULL,
    p_ord = NULL,
    p_derivs = NULL,
    p_bls = NULL,
    p_params = NULL,
    p_bltrace = NULL,
    p_log = NULL,
    p_loss = NULL,
    p_vinit = NULL,
    p_positive = NULL
  )
)
