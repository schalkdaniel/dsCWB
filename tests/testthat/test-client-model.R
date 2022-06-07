context("Client model")

#library(testthat)
#devtools::load_all()

test_that("client model can be initialized correctly",  {
  symbol = "iris"
  target = "Sepal.Width"
  learning_rate = 0.1
  cm = expect_silent(ClientModel$new(symbol, target, learning_rate = learning_rate))

  expect_equal(cm$getTarget(), iris[[target]])
  expect_equal(cm$getFeatureNames(), setdiff(names(iris), target))
  expect_true(checkmate::checkClass(cm$loss, "LossQuadratic"))
  expect_null(cm$getPrediction())
  expect_equal(cm$getBlParams(), list())
  expect_null(cm$getBlNames())
  expect_null(cm$getPseudoResids())
  expect_equal(cm$getXty(), list())
  expect_equal(cm$getXtX(), list())
  expect_error(cm$getRisk())
  expect_equal(cm$getLossOptimalConstant(), mean(iris[[target]]))
})

test_that("client model initializes correctly", {
  symbol = "iris"
  target = "Sepal.Width"
  learning_rate = 0.1
  cm = expect_silent(ClientModel$new(symbol, target, learning_rate = learning_rate))

  expect_silent(cm$initConstantModel(cm$getLossOptimalConstant()))
  expect_equal(cm$getPrediction(), rep(mean(iris[[target]]), nrow(iris)))
  expect_equal(cm$getPseudoResids(), cm$loss$pseudoResids(cm$getTarget(), cm$getPrediction()))
  expect_equal(cm$getRisk(), cm$loss$risk(cm$getTarget(), cm$getPrediction()))

  init = expect_silent(getClientInit(symbol, cm$getFeatureNames()))
  expect_silent(cm$addBaselearner(init))

  ll_xtx = expect_silent(cm$getXtX())
  ll_xty = expect_silent(cm$getXty())

  for (bl in cm$bls) {
    x = iris[[bl$getFeature()]]
    X = bl$basisTrafo(x)
    XtX = t(X) %*% X
    Xty = as.vector(t(X) %*% cm$getPseudoResids())
    par = as.vector(solve(XtX + bl$getHyperpars()$penalty * bl$getPenaltyMat()) %*% Xty)
    pred = as.vector(X %*% par)

    ## Visual inspection:
    #ord = order(x)
    #plot(y = cm$getTarget(), x = x)
    #lines(y = pred[ord] + cm$getOffset(), x = x[ord])

    expect_equal(pred, bl$predict(par))
    expect_equal(par, bl$train(cm$getPseudoResids()))
    sse = sum((cm$getPseudoResids() - pred)^2)
    expect_equal(ll_xtx[[paste0(bl$getFeature(), "-spline")]], XtX)
    expect_equal(ll_xty[[paste0(bl$getFeature(), "-spline")]], Xty)
    expect_equal(cm$getSSE()$from_bl[[paste0(bl$getFeature(), "-spline")]], sse)
  }
})

test_that("client model generation works on DataSHIELD server", {
  surl     = "https://opal-demo.obiba.org/"
  username = "administrator"
  password = "password"

  pkg = "dsCWB"

  opal = expect_silent(dsCWB:::.tryOPALConnection(opalr::opal.login(username = username, password = password, url = surl)))

  if (inherits(opal, "opal")) {
    # Check if package can be installed:
    expect_true(opalr::dsadmin.install_github_package(opal = opal, pkg = pkg, username = "schalkdaniel", ref = "main"))
    expect_true(opalr::dsadmin.publish_package(opal = opal, pkg = pkg))

    opalr::opal.logout(opal, save = FALSE)

    library(DSI)
    library(DSOpal)

    builder = newDSLoginBuilder()

    builder$append(
      server   = "ds-test-server-dummy1",
      url      = surl,
      user     = username,
      password = password
    )
    builder$append(
      server   = "ds-test-server-dummy2",
      url      = surl,
      user     = username,
      password = password
    )
    connections <<- datashield.login(logins = builder$build(), assign = TRUE)

    datashield.assign(connections, "dat", quote(iris))
    expect_silent(suppressMessages(datashield.assign(connections, "cm", quote(createClientModel("dat", "Petal.Length")))))

    datashield.symbols(connections)

    fn = datashield.aggregate(connections, quote(getFeatureNames("cm")))
    fns = c("Species", "Sepal.Length", "Sepal.Width", "Petal.Width")
    call_get_init = paste0("getClientInit(\"dat\", ", encodeObject(fns), ")")
    cl_init = datashield.aggregate(connections, call_get_init)

    call_init = paste0("initClientModel(\"cm\", ", encodeObject(cl_init), ")")
    datashield.assign(connections, call_init)

    cinit = datashield.aggregate(connections, quote(getOptimalConstant("cm")))
    call_cinit = paste0("initSiteConstant(\"cm\", ", Reduce("+", cinit), ")")

    datashield.assign(connections, call_cinit)

    datashield.aggregate(connections, quote(getClientXtX("cm")))
    datashield.aggregate(connections, quote(getClientXty("cm")))
    datashield.aggregate(connections, quote(getClientSSE("cm")))

    expect_equal(dsBinVal:::.dsLength(connections, "valid"), length(valid) * 2)

    p_cls <<- ifelse(p > 0.5, 1, 0)
    conf_local = table(truth = valid, predicted = p_cls)
    expect_equal(confusion("valid", "p_cls"), conf_local)
    conf = expect_silent(suppressMessages(dsConfusion(connections, "valid", "pred")))
    expect_equal(nrow(conf$confusion), 2)
    expect_equal(nrow(conf$confusion), 2)

    expect_equal(l2sens("iris", "p", nbreaks = 30L)$l2sens, dsL2Sens(connections, "dat", "pred", nbreaks = 30L))
    expect_silent(suppressMessages({
      roc_glm = dsROCGLM(connections, "valid", "pred", dat_name = "iris",
        seed_object = "pred")
    }))
    expect_equal(class(roc_glm), "ROC.GLM")
    expect_output(print(roc_glm))

    expect_silent(suppressMessages({
      roc_glm2 = dsROCGLM(connections, "valid", "pred", dat_name = "iris",
        seed_object = "pred")
    }))
    expect_equal(roc_glm, roc_glm2)

    # Suppress ggplot warnings:
    gg = expect_silent(suppressWarnings(suppressMessages(plot(roc_glm))))
    expect_true(inherits(gg, "ggplot"))
    expect_output(print(roc_glm))


    datashield.assign(connections, "dat_no_na", quote(removeMissings("dat")))
    nuisance = lapply(DSI::datashield.symbols(connections), function(s) {
       expect_true("dat_no_na" %in% s)
    })

    ri = datashield.aggregate(connections, quote(getDataSHIELDInfo()))
    expect_equal(class(ri), "list")
    nuisance = lapply(ri, function(r) {
      expect_equal(names(r), c("session_info", "pcks"))
    })

    # Weird, sometimes it complains that message is printed and sometimes that it does
    # not produce messages ...
    cc = expect_silent(suppressMessages(dsCalibrationCurve(connections, "valid", "pred", 10, 3)))
    expect_output(print(cc))

    expect_error(brierScore(connections, 1, 2))
    bs = expect_silent(suppressMessages(dsBrierScore(connections, "valid", "pred")))
    expect_true(is.numeric(bs))

    gg_cc = expect_silent(suppressMessages(plot(cc)))
    expect_true(inherits(gg_cc, "ggplot"))

    datashield.logout(connections)

  }
})
