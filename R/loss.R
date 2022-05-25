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
    #' @param ll_inits (`list()`)\cr
    #'   List with initialization values.
    aggregateInit = function(ll_inits) {
      checkmate::assertList(ll_inits)
      n = length(ll_inits[[1]])
      nuisance = lapply(ll_inits, checkmate::assertNumeric, len = n, any.missing = FALSE)
      return(Reduce("+", ll_inits))
    }
  )
)
