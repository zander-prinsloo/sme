


test_that("sme_checks() works", {

  # don't include t3
  expect_error(
    sme_checks(
      St1_observed = 1,
      St2_observed = 1,
      St2_true     = 1,
      St1_true     = 1,
      theta_01     = 1.5,
      theta_02     = -1.5,
      err          = 1.65,
      mu           = 0.2
    )
  )
  expect_error(
    sme_checks(
      St1_observed = 1,
      St2_observed = 1,
      St2_true     = 1,
      St1_true     = 1,
      theta_01     = 1.5,
      theta_02     = -1.5,
      err          = 1.65,
      mu           = 0.2,
      only_2       = T
    )
  )


  # Status not 0-1
  expect_error(
    sme_checks(
      1, 1, 2, 1, 1, 1,
      1.5, -1.5, 1.65, 0.2
    )
  )
  # Status T or F
  expect_no_error(
    sme_checks(
      T, T, T, F, T, T,
      1.5, -1.5, 1.65, 0.2
    )
  )
  # Paramater T or F
  expect_error(
    sme_checks(
      1, 1, 1, 1, 1, 1,
      F, -1.5, 1.65, 0.2
    )
  )

})
