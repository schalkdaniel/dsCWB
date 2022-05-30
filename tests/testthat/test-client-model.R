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
