aggregateInit = function(ll_init) {
  ss = names(ll_init)
  fn = names(ll_init[[1]])
  ll_out = list()
  for (f in fn) {
    mm = c()
    ma = c()
    cls = c()
    for (s in ss) {
      if (ll_init[[s]][[f]]$class == "numeric") {
        mm = min(mm, ll_init[[s]][[f]]$min)
        ma = max(ma, ll_init[[s]][[f]]$max)
      }
      if (ll_init[[s]][[f]]$class == "categorical") {
        cls = unique(c(cls, ll_init[[s]][[f]]$table))
      }
    }
    if (ll_init[[1]][[f]]$class == "numeric") ll_out[[f]] = list(min = mm, max = ma, class = "numeric")
    if (ll_init[[1]][[f]]$class == "categorical") ll_out[[f]] = list(table = cls, class = "categorical")
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

getMinimalSSE = function(ll_sses, force_shared = FALSE) {
  checkmate::assertList(ll_sses)
  nuisance = lapply(ll_sses, function(lls) {
    checkmate::assertList(lls)
    nuisance = lapply(names(lls), function(nll) {
      checkmate::assertChoice(nll, c("from_list", "from_bl"))})
  })

  sses_aggr = ll_sses[[1]]
  out = list()
  for (n1 in names(sses_aggr)) {
    for (n2 in names(sses_aggr[[n1]])) {
      out[[n1]][[n2]] = 0
    }
  }
  for (n1 in names(ll_sses)) {
    for (n2 in names(ll_sses[[n1]])) {
      for (n3 in names(ll_sses[[n1]][[n2]])) {
        out[[n2]][[n3]] = out[[n2]][[n3]] + ll_sses[[n1]][[n2]][[n3]]
      }
    }
  }
  min_list_idx = which.min(unlist(out$from_list))
  min_bl_idx = which.min(unlist(out$from_bl))

  if (length(min_list_idx) == 0)
    from_list = NULL
  else
    from_list = out$from_list[[min_list_idx]]

  if (length(min_bl_idx) == 0)
    from_bl = NULL
  else
    from_bl = out$from_bl[[min_bl_idx]]


  combined_min = c(from_list, from_bl)
  names(combined_min) = c(names(min_list_idx), names(min_bl_idx))

  min_idx = which.min(combined_min)

  # It may happen due to the transmission of the sse values that there are small
  # differences even though there aren't. Hence, a check with the diff is done
  # instead of checking `combined_min[[1]] == combined_min[[2]]`.
  if ((abs(diff(combined_min)) < 1e-9) || force_shared) {
    min_idx = 1
    names(min_idx) = names(combined_min)[1]
  }

  sse_min = min(c(from_list, from_bl))
  names(sse_min) = names(min_idx)

  attr(sse_min, "effect_type") = ifelse(min_idx == 1, "shared", "site")

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
    cchar = "NULL"
  else
    cchar = paste0(qchar, x, qchar)

  return(paste(cchar, collapse = ", "))
}

calculateDF = function(xtxs, hm, df) {
  checkmate::assertList(xtxs)
  checkmate::assertR6(hm, "HostModel")
  checkmate::assertNumeric(df, lower = 1, len = 1L, any.missing = FALSE)

  bls = hm$bls
  blnames = names(bls)
  pens = list()

  for (bln in blnames) {
    llmat = list()
    llpen = list()
    for (sn in names(xtxs)) {
      llmat[[sn]] = xtxs[[sn]][[bln]]
      llpen[[sn]] = bls[[bln]]$getPenaltyMat()
    }
    xtx_tensor = Matrix::bdiag(llmat)
    penmat_tensor = Matrix::bdiag(llpen)
    pens[[bln]] = cpsp::demmlerReinsch(as.matrix(xtx_tensor), as.matrix(penmat_tensor), df)
  }
  return(pens)
}

#' @title Train a distributed CWB model.
#' @param connections (`list(OpalConnection)`)\cr
#'   Connections to the DataSHIELD servers (see `?DSI::newDSLoginBuilder`).
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
#' @param trace (`logical(1L)`)\cr
#'   Indicator if the fitting trace should be printed or not.
#' @param force_shared_iters (`integer(1L)`)\cr
#'   Number of iterations in the beginning that are just training the shared effect.
#' @param random_intercept (`logical(1L)`)\cr
#'   Indicator whether random intercepts should be added or not.
#' @param random_intercept_df (`integer(1L)`)\cr
#'   Penalty of the random intercept.
#' @return Client model of R6 class ClientModel.
#' @importFrom DSI datashield.aggregate datashield.assign
#' @examples
#' \dontrun{
#' library(DSI)
#' library(DSOpal)
#'
#' # Establish connection to the DataSHIELD servers:
#' surl     = "https://opal-demo.obiba.org/"
#' username = "administrator"
#' password = "password"
#'
#' builder = newDSLoginBuilder()
#'
#' builder$append(
#'   server   = "ds-test-server-dummy1",
#'   url      = surl,
#'   user     = username,
#'   password = password
#' )
#' builder$append(
#'   server   = "ds-test-server-dummy2",
#'   url      = surl,
#'   user     = username,
#'   password = password
#' )
#' connections = datashield.login(logins = builder$build(), assign = TRUE)
#' datashield.assign(connections, "dat", quote(iris))
#'
#' feature_names = setdiff(names(iris), "Sepal.Length")
#' cwb = dsCWB(connections, "dat", "Sepal.Length", feature_names)
#' cwb$getLog()
#'
#' datashield.logout(connections)
#' }
#' @export
dsCWB = function(connections, symbol, target = NULL, feature_names, mstop = 100L,
  learning_rate = 0.1, df = 5, nknots = 20L, ord = 3L, derivs = 2L, val_fraction = NULL,
  patience = NULL, eps_for_break = 0, positive = NULL, seed = NULL, trace = TRUE,
  force_shared_iters = NULL, random_intercept = TRUE, random_intercept_df = 2) {

  checkConnection(connections)

  checkmate::assertCharacter(x = symbol, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(x = target, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(x = feature_names, any.missing = FALSE)
  checkmate::assertCount(x = mstop, positive = TRUE)
  checkmate::assertNumeric(x = learning_rate, len = 1L, any.missing = FALSE)
  checkmate::assertNumeric(x = df, any.missing = FALSE, lower = 0)
  checkmate::assertIntegerish(x = nknots, len = 1L, any.missing = FALSE)
  checkmate::assertIntegerish(x = ord, len = 1L, any.missing = FALSE)
  checkmate::assertIntegerish(x = derivs, len = 1L, any.missing = FALSE)
  checkmate::assertNumeric(val_fraction, lower = 0, upper = 1, len = 1L, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertCount(x = patience, null.ok = TRUE, positive = TRUE)
  checkmate::assertNumeric(x = eps_for_break, lower = 0, upper = 1, len = 1L, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertCharacter(x = positive, len = 1L, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertCount(x = seed, null.ok = TRUE)
  checkmate::assertCount(x = force_shared_iters, null.ok = TRUE)

  if (random_intercept) {
    cnms = names(dsBaseClient::ds.names(symbol, connections))
    for (i in seq_along(connections)) {
      dsBaseClient::ds.make(transChar(cnms[i]), "SNAME", connections[i])
    }
  }

  if (length(df) > 2) stop("df must be of length <= 2")

  symchar = transChar(symbol)
  tchar = transChar(target)
  fn = transChar(feature_names)
  oobchar = transChar(val_fraction)
  pchar = transChar(positive)
  schar = transChar(seed)

  model_symbol = "cm"
  mchar = transChar(model_symbol)

  tt = datashield.aggregate(connections, paste0("getClientTaskType(", symchar, ",", tchar, ")"))

  if (length(unique(unlist(tt))) > 1)
    stop("Clients have returned different task types.")

  tt = tt[[1]]

  if (length(df) == 2)
    dfs = df[1]
  else
    dfs = df

  hm = HostModel$new(symbol = symbol, target = target, target_type = tt, feature_names = feature_names,
    learning_rate = learning_rate, df = dfs, nknots = nknots, ord = ord, derivs = derivs,
    positive = positive)

  ## Create client models:
  ## ======================================

  call_init_client_model = NULL
  eval(parse(text = paste0("call_init_client_model = quote(createClientModel(",
    symchar, ",", tchar, ", c(", fn, ") ,", learning_rate, ", ", dfs, ", ", nknots, ", ", ord, ", ",
    derivs, ", ", oobchar, ", ", pchar, ", ", schar, ", ", random_intercept, ", ", random_intercept_df,
  "))")))
  datashield.assign(connections, model_symbol, call_init_client_model)

  call_get_init = paste0("getClientInit(", symchar, ", \"", encodeObject(feature_names), "\")")
  cl_init = datashield.aggregate(connections, call_get_init)
  cli = aggregateInit(cl_init)

  if (random_intercept)
    cli$random_intercept = list(table = cnms, class = "categorical")


  cq = NULL
  call_init = paste0("initClientModel(", mchar, ", \"", encodeObject(cli), "\")")
  eval(parse(text = paste0("cq = quote(", call_init, ")")))
  datashield.assign(connections, model_symbol, cq)

  cl_opt_const = NULL
  eval(parse(text = paste0("cl_opt_const = quote(getOptimalConstant(", mchar ,"))")))
  cinit = datashield.aggregate(connections, cl_opt_const)

  w_opt_const = NULL
  eval(parse(text = paste0("w_opt_const = quote(getClientTrainValObs(", mchar ,"))")))
  winit = datashield.aggregate(connections, w_opt_const)

  ntrain = vapply(winit, function(w) w$ntrain, integer(1L))

  co = hm$loss$aggregateInit(unlist(cinit), unlist(ntrain) / Reduce("+", ntrain))

  cq_cinit = NULL
  cl_site_const_init = paste0("initSiteConstant(", mchar, ", ", co, ")")
  eval(parse(text = paste0("cq_cinit = quote(", cl_site_const_init, ")")))
  datashield.assign(connections, model_symbol, cq_cinit)

  cl_site_xtx = paste0("getClientXtX(\"cm\")")
  xtxs = datashield.aggregate(connections, cl_site_xtx)
  ll_xtx = aggregateXtX(xtxs)

  # Init host model:
  hm$setOffset(co)
  if (random_intercept) {
    ll_tmp = ll_xtx[which(! names(ll_xtx) == "random_intercept")]
    cli_tmp = cli[which(! names(cli) == "random_intercept")]
  } else {
    ll_tmp = ll_xtx
    cli_tmp = cli
  }

  hm$addBaselearners(cli_tmp, ll_tmp)

  cl_pen_update = NULL
  if (length(df) == 2)
    dfa = df[2]
  else
    dfa = df

  penalty_global = lapply(hm$bls, FUN = function(bl) bl$getPenalty())

  ns = ncol(ll_xtx$random_intercept)
  pri = cpsp::demmlerReinsch(as.matrix(ll_xtx$random_intercept), diag(ns), random_intercept_df)
  if (random_intercept) {
    penalty_global$random_intercept = pri
  }
  eval(parse(text = paste0("cl_pen_update = quote(updateClientPenalty(", mchar, ", ",
    transChar(encodeObject(penalty_global)), ", anisotrop = FALSE, simple = TRUE, update_random_intercept = TRUE))")))
  datashield.assign(connections, model_symbol, cl_pen_update)

  eval(parse(text = paste0("cl_pen_update = quote(updateClientPenalty(", mchar, ", ",
    transChar(encodeObject(pri)), "))")))
  datashield.assign(connections, model_symbol, cl_pen_update)

  # Training the model:
  train_iter = TRUE
  risk_old = Inf
  early_stop = FALSE
  k = 1
  if (all(! is.null(patience), ! is.null(eps_for_break), ! is.null(val_fraction))) {
    p0 = 0
    early_stop = TRUE
    risk_val_old = Inf
  }

  # Get offset risk:
  # Get risk:
  ll_rt = datashield.aggregate(connections, eval(parse(text = paste0("quote(getRisk(", mchar, ", \"train\"))"))))
  risk_train_offset = Reduce("+", ll_rt) / sum(unlist(ntrain))

  risk_val_offset = NA
  if (! is.null(val_fraction)) {
    nval = vapply(winit, function(w) w$nval, integer(1L))
    ll_rv = datashield.aggregate(connections, eval(parse(text = paste0("quote(getRisk(", mchar, ", \"val\"))"))))
    risk_val_offset = Reduce("+", ll_rv) / sum(unlist(nval))
  }

  hm$log("_intercept", "shared", NA, risk_train_offset, risk_val_offset)
  if (is.null(force_shared_iters)) force_shared_iters = 0

  while (train_iter) {
    # Get Xty and SSEs from fitted site-specific effects from the ds servers:
    ll_xty = datashield.aggregate(connections, paste0("getClientXty(", mchar, ")"))

    if ("random_intercept" %in% names(ll_xty[[1]])) {
      for (sn in names(ll_xty)) ll_xty[[sn]]$random_intercept = NULL
    }

    ll_shared_effects_param = hm$getParam(aggregateXtX(ll_xty))

    cl_sses = paste0("getClientSSE(", mchar, ", \"", encodeObject(ll_shared_effects_param), "\")")
    ll_sses = datashield.aggregate(connections, cl_sses)

    if (k <= force_shared_iters)
      force_shared = TRUE
    else
      force_shared = FALSE

    min_sse = getMinimalSSE(ll_sses, force_shared)

    blname = transChar(names(min_sse))
    if (attr(min_sse, "effect_type") == "site") {
      # update client
      cl_update = paste0("quote(updateClientBaselearner(", mchar, ", ", blname, "))")
      datashield.assign(connections, model_symbol, eval(parse(text = cl_update)))
    } else {
      par = ll_shared_effects_param[[names(min_sse)]]
      hm$update(names(min_sse), par)
      #  update client parts.
      cl_update = paste0("quote(updateClientBaselearner(", mchar, ", ", transChar(names(min_sse)), ", ",
        transChar(encodeObject(par)), "))")
      datashield.assign(connections, model_symbol, eval(parse(text = cl_update)))
    }

    # Get risk:
    ll_rt = datashield.aggregate(connections, eval(parse(text = paste0("quote(getRisk(", mchar, ", \"train\"))"))))
    risk_train = Reduce("+", ll_rt) / sum(unlist(ntrain))
    risk_val = NA
    if (! is.null(val_fraction)) {
      nval = vapply(winit, function(w) w$nval, integer(1L))
      ll_rv = datashield.aggregate(connections, eval(parse(text = paste0("quote(getRisk(", mchar, ", \"val\"))"))))
      risk_val = Reduce("+", ll_rv) / sum(unlist(nval))
    }

    # Add to log:
    hm$log(names(min_sse), unname(attr(min_sse, "effect_type")), unname(min_sse), risk_train, risk_val)



    # Stop if maximal number of iterations is reached:
    if (k >= mstop) train_iter = FALSE

    # Stop if early stopping hits:
    if (early_stop) {
      if (is.infinite(risk_val_old))
        val_eps = 1
      else
        val_eps = (risk_val_old - risk_val) / risk_val_old

      if (val_eps <= eps_for_break)
        p0 = p0 + 1
      else
        p0 = 0

      if (p0 == patience) train_iter = FALSE
    }

    if (trace) {
      lline = hm$getLog(k)
      vrisk = ""
      pstring = ""

      ctime = paste0("[", Sys.time(), "] ")
      if (! is.na(lline$risk_val)) vrisk = paste0("risk (val) = ", round(lline$risk_val, 4), ", ")
      if (early_stop) pstring = paste0("patience = ", p0)
      cat(ctime, " ", lline$iteration, ": risk (train) = ", lline$risk_train, ", ", vrisk,
        "bl = ", lline$bl, " (", lline$effect_type, ", ", round(lline$sse, 4), "), ",
        pstring, "\n", sep = "")
    }

    risk_val_old = risk_val
    k = k + 1
  }
  site_params = datashield.aggregate(connections, paste0("getClientModelCoefficients(", mchar, ")"))
  hm$setSiteCoefficients(site_params)
  return(hm)
}

