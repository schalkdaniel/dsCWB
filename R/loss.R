# Init function e.g. sum/mean for l2 loss
# Aggregate function which takes the init parts and aggregates it.

#' @title Get the task type of a vector
#' @param x (`vector()`)\cr
#'   Vector from which the task is extracted. Just regression and binary classification is supported.
#' @return Character containing the task type.
#' @export
getTaskType = function(x) {
  if (is.numeric(x)) return("regression")
  if (is.character(x) || is.factor(x)) {
    nclasses = length(unique(x))
    if (nclasses != 2) stop("Just binary classification is supported! Your target vector has ", nclasses, " classes!")
    return("bin-class")
  }
  stop("Cannot identify the task!")
}

#' @title Quadratic loss class
#'
#' @description
#' This class defines the loss and the specific functions required for boosting.
#'
#' @export
LossQuadratic = R6Class("LossQuadratic",
  public = list(

    #' @description
    #' Get the constant initialization of the loss.
    #' @param truth (`numeric()`)\cr
    #'   Response vector.
    init = function(truth) {
      checkmate::assertNumeric(truth, any.missing = FALSE)
      return(mean(truth))
    },

    #' @description
    #' Calculate the pointwise loss.
    #' @param truth (`numeric()`)\cr
    #'   The response variable.
    #' @param pred (`numeric()`)\cr
    #'   The predicted values for the response.
    loss = function(truth, pred) {
      checkmate::assertNumeric(truth, any.missing = FALSE)
      checkmate::assertNumeric(pred, any.missing = FALSE)
      if (! length(pred) %in% c(1, length(truth)))
        stop("Assertion on 'pred' failed: Must have length ", length(truth), " or 1, but has length ", length(pred))

      return(0.5 * (truth - pred)^2)
    },

    #' @description
    #' Calculate the pseudo residuals.
    #' @param truth (`numeric()`)\cr
    #'   The response variable.
    #' @param pred (`numeric()`)\cr
    #'   The predicted values for the response.
    pseudoResids = function(truth, pred) {
      checkmate::assertNumeric(truth, any.missing = FALSE)
      checkmate::assertNumeric(pred, any.missing = FALSE)
      if (! length(pred) %in% c(1, length(truth)))
        stop("Assertion on 'pred' failed: Must have length ", length(truth), " or 1, but has length ", length(pred))

      return(truth - pred)
    },

    #' @description
    #' Calculate the risk.
    #' @param truth (`numeric()`)\cr
    #'   The response variable.
    #' @param pred (`numeric()`)\cr
    #'   The predicted values for the response.
    risk = function(truth, pred) {
      checkmate::assertNumeric(truth, any.missing = FALSE)
      checkmate::assertNumeric(pred, any.missing = FALSE)
      if (! length(pred) %in% c(1, length(truth)))
        stop("Assertion on 'pred' failed: Must have length ", length(truth), " or 1, but has length ", length(pred))

      return(sum(self$loss(truth, pred)))
    },

    #' @description
    #' Aggregate initializations.
    #' @param inits (`numeric()`)\cr
    #'   Initialization values.
    #' @param w (`numeric()`)\cr
    #'   Weights.
    aggregateInit = function(inits, w) {
      checkmate::assertNumeric(inits, any.missing = FALSE)
      checkmate::assertNumeric(w, any.missing = FALSE, len = length(inits))
      return(stats::weighted.mean(x = inits, w = w))
    }
  )
)

#' @title Binomial loss class
#'
#' @description
#' This class defines the loss and the specific functions required for boosting.
#'
#' @export
LossBinomial = R6Class("LossBinomial",
  public = list(

    #' @description
    #' Get the constant initialization of the loss.
    #' @param truth (`numeric()`)\cr
    #'   Response vector.
    init = function(truth) {
      checkmate::assertIntegerish(truth, any.missing = FALSE, lower = -1, upper = 1)
      if (any(truth == 0)) stop("Target must be encoded as -1 and 1")
      return(mean(truth == 1))
      #p = mean(truth)
      #return(0.5 * log(p / (1 - p)))
    },

    #' @description
    #' Calculate the pointwise loss.
    #' @param truth (`numeric()`)\cr
    #'   The response variable.
    #' @param pred (`numeric()`)\cr
    #'   The predicted values for the response.
    loss = function(truth, pred) {
      checkmate::assertIntegerish(truth, any.missing = FALSE, lower = -1, upper = 1)
      if (any(truth == 0)) stop("Target must be encoded as -1 and 1")
      checkmate::assertNumeric(pred, any.missing = FALSE)
      if (! length(pred) %in% c(1, length(truth)))
        stop("Assertion on 'pred' failed: Must have length ", length(truth), " or 1, but has length ", length(pred))

      return(log(1 + exp(-truth * pred)))
    },

    #' @description
    #' Calculate the pseudo residuals.
    #' @param truth (`numeric()`)\cr
    #'   The response variable.
    #' @param pred (`numeric()`)\cr
    #'   The predicted values for the response.
    pseudoResids = function(truth, pred) {
      checkmate::assertIntegerish(truth, any.missing = FALSE, lower = -1, upper = 1)
      if (any(truth == 0)) stop("Target must be encoded as -1 and 1")
      checkmate::assertNumeric(pred, any.missing = FALSE)
      if (! length(pred) %in% c(1, length(truth)))
        stop("Assertion on 'pred' failed: Must have length ", length(truth), " or 1, but has length ", length(pred))

      return(truth / (1 + exp(truth * pred)))
    },

    #' @description
    #' Calculate the risk.
    #' @param truth (`numeric()`)\cr
    #'   The response variable.
    #' @param pred (`numeric()`)\cr
    #'   The predicted values for the response.
    risk = function(truth, pred) {
      checkmate::assertNumeric(truth, any.missing = FALSE)
      checkmate::assertNumeric(pred, any.missing = FALSE)
      if (! length(pred) %in% c(1, length(truth)))
        stop("Assertion on 'pred' failed: Must have length ", length(truth), " or 1, but has length ", length(pred))

      return(sum(self$loss(truth, pred)))
    },

    #' @description
    #' Aggregate initializations.
    #' @param inits (`numeric()`)\cr
    #'   Initialization values.
    #' @param w (`numeric()`)\cr
    #'   Weights.
    aggregateInit = function(inits, w) {
      checkmate::assertNumeric(inits, any.missing = FALSE)
      checkmate::assertNumeric(w, any.missing = FALSE, len = length(inits))
      p = stats::weighted.mean(x = inits, w = w)
      return(log(p / (1 - p)))
    }
  )
)

