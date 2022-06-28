context("Client model")

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

  init = expect_silent(getClientInit(symbol, encodeObject(cm$getFeatureNames())))
  expect_silent(cm$addBaselearners(init))

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
    if (bl$getType() == "numeric") {
      addon = "-spline"
    } else {
      addon = "-onehot"
    }

    expect_equal(ll_xtx[[paste0(bl$getFeature(), addon)]], XtX)
    expect_equal(ll_xty[[paste0(bl$getFeature(), addon)]], Xty)
    expect_equal(cm$getSSE()$from_bl[[paste0(bl$getFeature(), addon)]], sse)
  }

  pennew = lapply(cm$bls, function(bl) bl$getHyperpars()$penalty)
  pennew[[1]] = 1
  pennew[[2]] = 2
  pennew[[3]] = 3
  pennew[[4]] = 1
  expect_silent(cm$updatePenalty(pennew))
  expect_equal(cm$bls[[1]]$getHyperpars()$penalty, 1)
  expect_equal(cm$bls[[2]]$getHyperpars()$penalty, 2)
  expect_equal(cm$bls[[3]]$getHyperpars()$penalty, 3)
  expect_equal(cm$bls[[4]]$getHyperpars()$penalty, 1)

  pennew[["bla"]] = 4
  expect_error(cm$updatePenalty(pennew))
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

    expect_true("cm" %in% datashield.symbols(connections)[[1]])

    fn = datashield.aggregate(connections, quote(getFeatureNames("cm")))
    call_get_init = paste0("getClientInit(\"dat\", \"", encodeObject(fn[[1]]), "\")")
    cl_init = datashield.aggregate(connections, call_get_init)
    cli = expect_silent(dsCWB:::aggregateInit(cl_init))

    for (f in c("Sepal.Length", "Sepal.Width", "Petal.Width")) {
      expect_equal(cli[[f]]$min, min(iris[[f]]))
      expect_equal(cli[[f]]$max, max(iris[[f]]))
    }
    expect_equal(cli$Species[[1]], names(table(iris$Species)))

    call_init = paste0("initClientModel(\"cm\", \"", encodeObject(cli), "\")")
    eval(parse(text = paste0("cq = quote(", call_init, ")")))
    expect_silent(suppressMessages(datashield.assign(connections, "cm", cq)))

    cinit = expect_silent(suppressMessages(datashield.aggregate(connections, quote(getOptimalConstant("cm")))))
    co = Reduce("mean", cinit)
    expect_equal(co, mean(iris$Petal.Length))

    call_cinit = paste0("initSiteConstant(\"cm\", ", co, ")")
    eval(parse(text = paste0("cq = quote(", call_cinit, ")")))
    datashield.assign(connections, "cm", cq)

    xtxs = expect_silent(suppressMessages(datashield.aggregate(connections, quote(getClientXtX("cm")))))
    xtys = expect_silent(suppressMessages(datashield.aggregate(connections, quote(getClientXty("cm")))))
    sses = expect_silent(suppressMessages(datashield.aggregate(connections, quote(getClientSSE("cm")))))

    expect_equal(class(xtxs), "list")
    expect_equal(class(xtys), "list")
    expect_equal(class(sses), "list")

    for (sn in names(xtxs)) {
      for (fn in names(xtxs[[1]])) {
        expect_s4_class(xtxs[[sn]][[fn]], "Matrix")
        expect_type(xtys[[sn]][[fn]], "double")
        expect_type(sses[[sn]][["from_bl"]][[fn]], "double")
      }
    }
    datashield.logout(connections)
  }
})
