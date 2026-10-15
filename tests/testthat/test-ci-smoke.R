test_that("main unitregTMB functions are exported", {
  expect_true(is.function(unitregTMB::unitregTMB))
  expect_true(is.function(unitregTMB::unitregTMB.control))
})

test_that("default control configures serial execution", {
  ctrl <- unitregTMB::unitregTMB.control()
  expect_identical(ctrl$parallel$n, 1L)
})
