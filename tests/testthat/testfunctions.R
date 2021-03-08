library(persist)
context("Basic function values")

colonies1 <- matrix(c(1,1,10,10,1,10,1,10),ncol=2)


test_that("mnnd is zero for two identical set of coords", {
  expect_equal(mnnd(colonies1,colonies1), 0)
})

test_that("smmd is invariant to area if using same set of coors", {
  expect_equal(smnnd(colonies1,colonies1,1), smnnd(colonies1,colonies1,100))
})

test_that("cspscaled is equal to 1 if colonies don't move", {
  expect_equal(cspscaled(colonies1,colonies1), 1)
})
