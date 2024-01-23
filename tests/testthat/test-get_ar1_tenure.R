



test_that("ex_gaussian_density works",{

  expect_equal(
    ex_gaussian_density(
    x      = 0,
    sigma  = 0.1,
    lambda = 0.5
  ),
  gamlss.dist::dexGAUS(
    x     = 0,
    mu    = 0,
    sigma = 0.1,
    nu    = 0.5
  )
  )

})
