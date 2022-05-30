aggregateInit = function(ll_init) {
  ss = names(ll_init)
  fn = names(ll_init[[1]])
  ll_out = list()
  for (f in fn) {
    for (s in ss) {
      mm = c()
      ma = c()
      cls = c()
      if (ll_init[[s]][[f]]$class == "numeric") {
        mm = min(mm, ll_init[[s]][[f]]$min)
        ma = max(ma, ll_init[[s]][[f]]$max)
      }
      if (ll_init[[s]][[f]]$class == "categorical") {
        cls = unique(c(cls, ll_init[[s]][[f]]$table))
      }
    }
    if (ll_init[[1]][[f]]$class == "numeric") ll_out[[f]] = list(min = mm, max = ma)
    if (ll_init[[1]][[f]]$class == "categorical") ll_out[[f]] = list(table = cls)
  }
  return(ll_out)
}

aggregateXtX = function(xtxs) {
  ss = names(xtxs)
 fn = names(xtxs[[1]])
  ll_out = list()
  for (f in fn) {
    xtx = 0
    for (s in ss) {
      xtx = xtx + xtxs[[s]][[f]]
    }
    ll_out[[f]] = xtx
  }
  return(ll_out)
}

getMinimalSSE = function(ll_sses) {
  checkmate::assertList(ll_sses)
  nuisance = lapply(ll_sses, function(lls) {
    checkmate::assertList(lls)
    nuisance = lapply(names(lls), function(nll) {
      checkmate::assertChoice(nll, c("from_list", "from_bl"))})
  })

  sses_aggr = ll_sses[[1]]
  for (n1 in names(sses_aggr)) {
    for (n2 in names(sses_aggr[[n1]])) {
      sses_aggr[[n1]][[n2]] = 0
    }
  }
  for (n1 in names(ll_sses)) {
    for (n2 in names(ll_sses[[n1]])) {
      for (n3 in names(ll_sses[[n1]][[n2]])) {
        sses_aggr[[n2]][[n3]] = sses_aggr[[n2]][[n3]] + ll_sses[[n1]][[n2]][[n3]]
      }
    }
  }
  min_list_idx = which.min(unlist(sses_aggr$from_list))
  min_bl_idx = which.min(unlist(sses_aggr$from_bl))

  from_list = sses_aggr$from_list[[min_list_idx]]
  from_bl = sses_aggr$from_bl[[min_bl_idx]]

  combined_min = c(from_list, from_bl)
  names(combined_min) = c(names(min_list_idx), names(min_bl_idx))

  min_idx = which.min(combined_min)
  sse_min = min(c(from_list, from_bl))
  names(sse_min) = names(min_idx)

  attr(sse_min, "effect_type") = ifelse(unname(min_idx) == 1, "shared", "site")

  return(sse_min)
}


if (FALSE) {
  devtools::load_all()

  ll_init = list(ds1 = getClientInit("iris", names(iris)),
                 ds2 = getClientInit("iris", names(iris)))

  ll_new = aggregateInit(ll_init)
}

transChar = function(x) {
  if (is.numeric(x)) {
    checkmate::assertNumeric(x, any.missing = FALSE, null.ok = TRUE)
    qchar = ""
  } else {
    checkmate::assertCharacter(x, any.missing = FALSE, null.ok = TRUE)
    qchar = "\""
  }

  if (is.null(x))
    cchar = paste0(qchar, "NULL", qchar)
  else
    cchar = paste0(qchar, x, qchar)

  return(paste(cchar, collapse = ", "))
}

#' @title Train a distributed CWB model.
#' @param connections (``)\cr
#'   Connections to the DataSHIELD servers.
#' @param symbol (`character(1L)`)\cr
#'   Character containing the name of the data.
#' @param target (`character(1L)`)\cr
#'   Character containing the name of the target variable.
#' @param feature_names (`character()`)\cr
#'   Character vector of all target variables.
#' @param mstop (`integer(1L)`)\cr
#'   Number of boosting iterations.
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
#' @param patience (`integer(1L)`)\cr
#'   Number of consecutive iterations without improvement to stop the algorithm.
#' @param eps_for_break (`numeric(1L)`)\cr
#'   Minimal relative improvement without stopping the algorithm.
#' @param positive (`character(1L)`)\cr
#'   Name of the positive class in binary classification.
#' @param seed (`numeric(1L)`)\cr
#'   Seed for generating validation data (only applies when val_fraction is set).
#' @return Client model of R6 class ClientModel.
#' @importFrom DSI datashield.aggregate datashield.assign
#' @export
dsCWB = function(connections, symbol, target = NULL, feature_names, mstop = 100L,
  learning_rate = 0.1, df = 5, nknots = 20L, ord = 3L, derivs = 2L, val_fraction = NULL,
  patience = NULL, eps_for_break = NULL, positive = NULL, seed = NULL) {

  checkmate::assertCharacter(x = symbol, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(x = target, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(x = feature_names, any.missing = FALSE)
  checkmate::assertCount(x = mstop, positive = TRUE)
  checkmate::assertNumeric(x = learning_rate, len = 1L, any.missing = FALSE)
  checkmate::assertNumeric(x = df, len = 1L, any.missing = FALSE)
  checkmate::assertIntegerish(x = nknots, len = 1L, any.missing = FALSE)
  checkmate::assertIntegerish(x = ord, len = 1L, any.missing = FALSE)
  checkmate::assertIntegerish(x = derivs, len = 1L, any.missing = FALSE)
  checkmate::assertNumeric(val_fraction, lower = 0, upper = 1, len = 1L, any.missing = FALSE)
  checkmate::assertCount(x = patience)
  checkmate::assertNumeric(x = eps_for_break, lower = 0, upper = 1, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(x = positive, len = 1L, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertCount(x = seed)

  hm = HostModel$new(symbol = symbol, target = target, feature_names = feature_names,
    learning_rate = learning_rate, df = df, nknots = nknots, ord = ord, derivs = derivs,
    positive = positive, connections = connections)

  tchar = transChar(target)
  fn = transChar(feature_names)
  oobchar = transChar(val_fraction)
  pchar = transChar(positive)
  schar = transChar(seed)

  cl_init = paste0("getClientInit(\"", symbol, "\", c(", fn, "))")
  ll_init = datashield.aggregate(connections, cl_init)

  ll_init = aggregateInit(ll_init)
  # Create client models:


  # TODO! Adjust degrees of freedom for the client models, must be smaller than
  # the one used for the shared effects.
  cl_mod = paste0("createClientModel(\"", symbol, "\", ", tchar, ", c(", fn, "), ",
    learning_rate, ",", df, ",", nknots, ",", ord, ",", derivs, ",", oobchar, ",",
    pchar, ",", schar, ")")

  datashield.assign(connections, "cm", cl_mod)
  datashield.assign(connections, "cm", paste0(initClientModel("cm", encodeObject(ll_init))))

  cl_const_init = paste0("getOptimalConstant(\"cm\")")
  inits = datashield.aggregate(connections, cl_const_init)

  # Host model aggregates constant initializations:
  const_init = hm$loss$aggregateInit(inits)

  cl_site_const_init = paste0("initSiteConstant(\"cm\", ", const_init, ")")
  datashield.assign(connections, "cm", cl_site_const_init)

  cl_site_xtx = paste0("getClientXtX(\"cm\")")
  xtxs = datashield.aggregate(connections, cl_site_xtx)

  ll_xtx = aggregateXtX(xtxs)

  # Init host model:
  hm$setOffset(const_init)
  hm$addBaselearner(ll_init, ll_xtx)
  #hm$connect(connections)

  # Training the model:
  train_iter = TRUE
  risk_old = Inf
  k = 1
  while (train_iter) {
    # Get Xty and SSEs from fitted site-specific effects from the ds servers:
    ll_xty = datashield.aggregate(connections, quote(getClientXty("cm")))
    ll_shared_effects_param = hm$getParam(ll_xty)

    cl_sses = paste0("getClientSSE(\"cm\", ", encodeObject(ll_shared_effects_param), ")")
    ll_sses = datashield.aggregate(connections, paste0(getClientSSE("cm", cl_sses)))

    min_sse = getMinimalSSE(ll_sses)

    # Get risk:
    ll_rt = datashield.aggregate(connections, quote(getRisk("cm", "train")))
    risk_train = Reduce("+", ll_rt)

    risk_val = NA
    if (! is.null(val_fraction)) {
      ll_rv = datashield.aggregate(connections, quote(getRisk("cm", "val")))
      risk_val = Reduce("+", ll_rv)
    }

    # Add to log:
    hm$log(names(min_sse), attr(min_sse, "effect_type"), min_sse, risk_train, risk_val)

    # Update model on ds servers if it was selected:
    bl_param = NULL
    if (names(min_sse) %in% names(ll_shared_effects_param))
      bl_param = ll_shared_effects_param[[names(min_sse)]]

    hm$update(names(min_sse), bl_param)

    # Determine stopping criteria:
    if (is.infinite(risk_old))
      val_eps = 1
    else
      val_eps = (risk_old - risk_train) / risk_old

    if ((k >= mstop) || (val_eps <= eps_for_break)) train_iter = FALSE

    risk_old = risk_train
    k = k + 1
  }
  return(hm)

  ## Initialize client models:
  #datashield.assign(createClientModel())
  # Get XtX
  #datashield.aggregate()


  # Initialize host and client model

  # Initialize clients:
  # - Initialze prediction (loss optimal constant)
  # - Initialize base learner (e.g. share min/max vlaues for knots calculation)
  # - share XtX
  # Initialize host by aggregating XtX

  #while (is_not_finished) {

  # fit shared effects (host model)
  # fit client models

  # share SSE
  # select SSE
  # update host or client models
  # share validation loss
  # Check if algo should be stopped

  # After stop:
  # - Return params from clinet models
  # -
  # }
}

