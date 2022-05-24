context("Check if l2-sensitivity can be calculated")

test_that("l2-sensitivity works", {
  scores <<- rnorm(n = nrow(iris))

  expect_equal(l2sens(iris, scores), l2sens("iris", "scores"))
  expect_equal(class(l2sens("iris", "scores")), "list")
  expect_equal(class(l2sens("iris", "scores", nbreaks = 10)), "list")
  expect_equal(l2sens(iris, scores), l2sens("iris", "scores", col_names = colnames(iris)))
  expect_equal(l2sens(iris, scores), l2sens("iris", "scores", nbreaks = 50))

  expect_error(l2sens("bla", 1:10))
})
