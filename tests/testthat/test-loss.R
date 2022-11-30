context("Loss objects")

test_that("Quadratic loss works correctly", {
  nsim = 100L
  pr = runif(nsim, -10, 100)
  tr = pr + rnorm(nsim)

  l = expect_silent(LossQuadratic$new())

  expect_equal(l$loss(tr, pr), 0.5 * (pr - tr)^2)
  expect_equal(l$risk(tr, pr), sum(0.5 * (pr - tr)^2))
  expect_equal(l$pseudoResids(tr, pr), tr - pr)

  oc = mean(tr)
  ol = l$init(tr)
  expect_equal(round(oc, 4), round(ol, 4))

  trs = list(tr[1:34], tr[35:80], tr[81:100])
  inits = sapply(trs, l$init)
  w = sapply(trs, length)

  expect_equal(l$aggregateInit(inits, w), ol)
})

test_that("Binomial loss works correctly", {
  nsim = 100L
  tr = sample(c(-1, 1), nsim, TRUE)
  pr = runif(nsim, 0, 2) * tr * sample(c(-1, 1), prob = c(0.2, 0.8))

  l = expect_silent(LossBinomial$new())

  expect_error(l$init(-2, 2))
  expect_error(l$init(0, 2))
  expect_error(l$init(-1, "bla"))
  expect_error(l$loss(-2, 2))
  expect_error(l$loss(0, 2))
  expect_error(l$loss(-1, "bla"))
  expect_error(l$risk(-2, 2))
  expect_error(l$risk(0, 2))
  expect_error(l$risk(-1, "bla"))

  expect_equal(l$loss(tr, pr), log(1 + exp(-tr * pr)))
  expect_equal(l$risk(tr, pr), sum(log(1 + exp(-tr * pr))))
  expect_equal(l$pseudoResids(tr, pr), tr / (1 + exp(tr * pr)))

  oc = optimize(function(x) l$risk(tr, x), lower = -2, upper = 2)$minimum
  ol = log(l$init(tr) / (1 - l$init(tr)))
  expect_equal(round(oc, 4), round(ol, 4))

  trs = list(tr[1:34], tr[35:80], tr[81:100])
  inits = sapply(trs, l$init)
  w = sapply(trs, length)

  expect_equal(l$aggregateInit(inits, w), ol)
})
