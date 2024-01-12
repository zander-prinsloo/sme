# test the component functions



test_that("get_ar1_sym_joint_probability works", {

  expect_equal(
    round(get_ar1_sym_joint_probability(
    st2_observed = 1,
    st2_true     = 1,
    st1_observed = 1,
    st1_true     = 1,
    theta_01     = 5,
    theta_02     = 5,
    err          = 5,
    mu           = 0
  ), 1),
  0.5
  )

  expect_equal(
    round(get_ar1_sym_joint_probability(
      st2_observed = 1,
      st2_true     = 1,
      st1_observed = 1,
      st1_true     = 1,
      theta_01     = 5,
      theta_02     = 5,
      err          = 5,
      mu           = 5
    ), 1),
    1
  )


  expect_equal(
    round(get_ar1_sym_joint_probability(
      st2_observed = 1,
      st2_true     = 1,
      st1_observed = 1,
      st1_true     = 1,
      theta_01     = 5,
      theta_02     = 5,
      err          = 5,
      mu           = -5
    ), 1),
    0
  )
})





