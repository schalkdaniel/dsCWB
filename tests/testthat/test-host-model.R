context("Host model")

#git add --all
#git commit -m "updates [ci skip]"
#git push origin main

#q()
#R
#library(testthat)
#devtools::load_all()

test_that("host model can be initialized correctly",  {
  symbol = "iris"
  target = "Sepal.Width"
  fns = setdiff(names(iris), target)

  ## Define client model to get XtX, Xty, etc.:
  cm      = ClientModel$new(symbol, target)
  ll_init = getClientInit(symbol, encodeObject(cm$getFeatureNames()))
  cm$addBaselearners(ll_init)
  ll_xtx = cm$getXtX()
  cm$initConstantModel(mean(iris[[target]]))
  ll_xty = cm$getXty()
  fupdate = names(ll_xty)[1]
  cm$update(fupdate)
  ## -- end client model

  hm = expect_silent(HostModel$new(symbol, target, "regression", fns))

  expect_equal(hm$getDataSymbol(), symbol)
  expect_equal(hm$getFeatureNames(), fns)
  expect_null(hm$bls)
  expect_error({hm$bls = 2})

  expect_silent(hm$addBaselearners(ll_init, ll_xtx))
  pars = expect_silent(hm$getParam(ll_xty))
})

test_that("update step is working", {
  symbol = "iris"
  target = "Sepal.Width"
  fns = setdiff(names(iris), target)

  ## Define client model to get XtX, Xty, etc.:
  cm      = ClientModel$new(symbol, target)
  ll_init = getClientInit(symbol, encodeObject(cm$getFeatureNames()))
  cm$addBaselearners(ll_init)
  ll_xtx = cm$getXtX()
  cm$initConstantModel(mean(iris[[target]]))
  ll_xty = cm$getXty()
  fupdate = names(ll_xty)[1]
  cm$update(fupdate)
  ## -- end client model

  hm = HostModel$new(symbol, target, "regression", fns)
  hm$addBaselearners(ll_init, ll_xtx)
  pars = hm$getParam(ll_xty)

  expect_silent(hm$update(fupdate, pars[[fupdate]]))
  expect_equal(hm$bls[[fupdate]]$getParam(), pars[[fupdate]] * hm$lr)

  expect_null(hm$bls[[2]]$getParam())
  expect_null(hm$bls[[3]]$getParam())
})
