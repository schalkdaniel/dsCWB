#' @title Initialize a client model
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
#'   Fraction of samples used for validation.
#' @param positive (`character(1L)`)\cr
#'   Name of the positive class in binary classification.
#' @param seed (`numeric(1L)`)\cr
#'   Seed for generating validation data (only applies when val_fraction is set).
#' @return Client model of R6 class ClientModel.
#' @export
createClientModel = function(symbol, target, feature_names = NULL, learning_rate = 0.1, df = 5, nknots = 20L,
  ord = 3L, derivs = 2L, val_fraction = NULL, positive = NULL, seed = NULL) {

  cm = ClientModel$new(symbol, target, feature_names, learning_rate, df, nknots, ord, derivs,
    val_fraction, positive, seed)
  return(cm)
}


#' @title Initialize the site base learner.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @param ll_init_binary (`character(1L)`)\cr
#'   Binary encoded list of initialization values.
#' @export
initClientModel = function(model_symbol, ll_init_binary) {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(ll_init_binary, len = 1L, any.missing = FALSE)
  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)

  ll_init = decodeBinary(ll_init_binary)
  checkmate::assertList(ll_init)

  #list(feature = c(min = 1, max = 2))
  cm$addBaselearner(ll_init)

  return(cm)
}

#' @title Get feature names from sites.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @export
getFeatureNames = function(model_symbol) {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)
  return(cm$getFeatureNames())
}

#' @title Get risk from model.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @param risk_type (`character(1L)`)\cr
#'   The risk type, must be one of "full", "train" (default), or "val".
#' @export
getRisk = function(model_symbol, risk_type = "train") {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  checkmate::assertChoice(model_symbol, choices = c("full", "train", "val"))
  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)
  return(cm$getRisk(risk_type))
}


#' @title Initialize client model with constant.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @param init (`numeric(1L)`)\cr
#'   Constant initialization.
#' @export
initSiteConstant = function(model_symbol, init) {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)
  cm$initConstantModel(init)
  return(cm)
}


#' @title Get loss optimal constants per site.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @export
getOptimalConstant = function(model_symbol) {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)
  return(cm$getLossOptimalConstant())
}

#' @title Update model base learner
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @param blname (`character(1L)`)\cr
#'   Name of the base learner.
#' @param par_binary (`character(1L)`)\cr
#'   Binary encoded parameter used to update pseudo residuals and predictions.
#' @export
updateClientBaselearner = function(model_symbol, blname, par_binary = NULL) {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(blname, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(par_binary, len = 1L, any.missing = FALSE, null.ok = TRUE)

  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)

  if (is.null(par_binary))
    cm$update(blname)
  else
    cm$updateCWBParts(blname, par_binary = par_binary)

  return(cm)
}

#' @title Get base learner initialization
#' @param symbol (`character(1L)`)\cr
#'   Character containing the name of the data.
#' @param fn_binary (`character(1L)`)\cr
#'   Binary encoded character vector of all target variables.
#' @return List with initialization values.
#' @export
getClientInit = function(symbol, fn_binary) {

  checkmate::assertCharacter(symbol, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(fn_binary, len = 1L, any.missing = FALSE)

  feature_names = decodeBinary(fn_binary)
  checkmate::assertCharacter(feature_names, any.missing = FALSE)

  dat = eval(parse(text = symbol), envir = .GlobalEnv)
  ll_out = lapply(feature_names, function(ff) {
    if (is.numeric(dat[[ff]])) {
      return(list(class = "numeric", min = min(dat[[ff]]), max = max(dat[[ff]])))
    }
    if (is.factor(dat[[ff]]) || is.character(dat[[ff]])) {
      return(list(class = "categorical", table = as.character(unique(dat[[ff]]))))
    }
  })
  #ll_out = list()
  #for (ff in feature_names) {
    #if (is.numeric(dat[[ff]])) {
      #ll_out[[ff]] = list(class = "numeric", min = min(dat[[ff]]), max = max(dat[[ff]]))
    #}
    #if (is.factor(dat[[ff]]) || is.character(dat[[ff]])) {
      #ll_out[[ff]] = list(class = "categorical", table = as.character(unique(dat[[ff]])))
    #}
  #}
  names(ll_out) = feature_names
  return(ll_out)
}

#' @title Get XtX from the client models.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @return XtX matrices of the registered base learners.
#' @export
getClientXtX = function(model_symbol = "cwb") {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)
  return(cm$getXtX())
}

#' @title Get Xty from the client models.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @return Xty matrices of the registered base learners.
#' @export
getClientXty = function(model_symbol = "cwb") {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)
  return(cm$getXty())
}

#' @title Get the SSE from the client models.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @param par_list_binary (`list()`)\cr
#'   Encoded list of parameter values. Each entry must have a name
#'   that corresponds to a base learner.
#' @return SSE values of the of the registered base learners.
#' @export
getClientSSE = function(model_symbol = "cwb", par_list_binary = NULL) {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(par_list_binary, len = 1L, null.ok = TRUE)

  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)

  return(cm$getSSE(par_list_binary = par_list_binary))
}

#' @title Get the model coefficients.
#' @param model_symbol (`character(1L)`)\cr
#'   The name of the model at the server.
#' @return Model coefficients per base learner.
#' @export
getClientModelCoefficients = function(model_symbol = "cwb") {
  checkmate::assertCharacter(model_symbol, len = 1L, any.missing = FALSE)

  cm = eval(parse(text = model_symbol), envir = .GlobalEnv)
  return(cm$getBlParams())
}

if (FALSE) {
  q()
  R
  devtools::load_all()
  cm = createClientModel("iris", "Sepal.Length", val_fraction = 0.2, seed = 31415L)
  fn = getFeatureNames("cm")
  ll_init = getClientInit("iris", fn)
  cm = initClientModel("cm", encodeObject(ll_init))
  cm$getBlNames()

  cm$initConstantModel(mean(iris$Sepal.Length))
  sses = cm$getSSE()

  p1 = cm$getPrediction()
  cm$getPseudoResids()

  cm$update(names(which.min(sses$from_bl)))
  cm$getBlParams()

  p2 = cm$getPrediction()
  cm$getPseudoResids()

  ll_xty = cm$getXty()
  ll_xtx = cm$getXtX()

  p_bl = p2 - p1
  plt = data.frame(y = iris$Sepal.Length, x = iris$Petal.Length, pred = p1 + 10*p_bl)
  #library(ggplot2)
  #ggplot(plt, aes(x = x)) + geom_point(aes(y = y)) + geom_line(aes(y = pred))

  cm$getBlNames()
  par_list = list("Sepal.Width-spline" = rnorm(24))

  sses2 = getClientSSE("cm", encodeObject(par_list))
  getClientXtX("cm")
  getClientXty("cm")

  hm = HostModel$new("iris", "Sepal.Length", "regression", fn)
  hm$addBaselearner(ll_init, ll_xtx)
  ll_par = hm$getParam(ll_xty)

  sses3 = getClientSSE("cm", encodeObject(ll_par))
  ll_sses = list(ds1 = sses3, ds2 = sses3, ds3 = sses3)
  getMinimalSSE(ll_sses)
}
