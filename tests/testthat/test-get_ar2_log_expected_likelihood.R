
t11 <- 0.12
t12 <- 0.6
t21 <- 0.35
t22 <- 0.05
mu  <- 0.25

test_that("p_num_misclass works", {

  # edge cases----------------------
  expect_equal(
    p_num_misclass(err          = 0,
                   num_misclass = 0),
    1)
  expect_equal(
    p_num_misclass(err          = 1,
                   num_misclass = 0),
    0)

  # normal cases--------------------
  expect_equal(
    p_num_misclass(err          = 0.25,
                   num_misclass = 1),
    dbinom(x = 1,
           size = 4,
           prob = 0.25))
  expect_equal(
    p_num_misclass(err          = 0.1,
                   num_misclass = 0),
    0.9*0.9*0.9*0.9)
  expect_equal(
    p_num_misclass(err          = 0.1,
                   num_misclass = 0) +
      p_num_misclass(err          = 0.1,
                     num_misclass = 1) +
      p_num_misclass(err          = 0.1,
                     num_misclass = 2) +
      p_num_misclass(err          = 0.1,
                     num_misclass = 3) +
      p_num_misclass(err          = 0.1,
                     num_misclass = 4),
    1)

  # vectorize
  expect_equal(p_num_misclass(err          = 0.1,
                              num_misclass = c(0, 1, 2, 3, 4)),
               c(p_num_misclass(err          = 0.1,
                                num_misclass = 0),
                 p_num_misclass(err          = 0.1,
                                num_misclass = 1),
                 p_num_misclass(err          = 0.1,
                                num_misclass = 2),
                 p_num_misclass(err          = 0.1,
                                num_misclass = 3),
                 p_num_misclass(err          = 0.1,
                                num_misclass = 4)))


})


test_that("p_misclass works", {

  # edge cases----------------------
  expect_equal(
    p_misclass(s     = 1,
               strue = 1,
               err   = 0),
    1)
  expect_equal(
    p_misclass(s = 1,
               strue = 1,
               err = 1),
    0)
  expect_equal(
    p_misclass(s = 1,
               strue = 2,
               err = 0),
    0)
  expect_equal(
    p_misclass(s = 1,
               strue = 2,
               err = 1),
    1)

  # normal cases--------------------
  expect_equal(
    p_misclass(s = 1,
               strue = 1,
               err = 0.25),
    0.75)
  expect_equal(
    p_misclass(s = 1,
               strue = 0,
               err = 0.25),
    0.25)
  expect_equal(
    p_misclass(s = 1,
               strue = 1,
               err = 0.1),
    0.9)
  expect_equal(
    p_misclass(s = 1,
               strue = 0,
               err = 0.1),
    0.1)

 e1  <- 0.1
 s1  <- 1
 st1 <- 0
 s2  <- 1
 st2 <- 1
 expect_equal(
   p_misclass(s     = s1,
              strue = st1,
              err   = e1),
   e1)
 expect_equal(
   p_misclass(s     = s2,
              strue = st2,
              err   = e1),
   1 - e1)

 # vectorize
 expect_equal(p_misclass(s     = c(s1, s2),
                         strue = c(st1, st2),
                         err   = e1),
              c(p_misclass(s     = s1,
                           strue = st1,
                           err   = e1),
                p_misclass(s     = s2,
                           strue = st2,
                           err   = e1)))


})

# TRUE 2 PERIOD
#--------------------------------

test_that("p_true_cond_2period works", {

  t11 <- 0.12
  t12 <- 0.6
  t21 <- 0.35
  t22 <- 0.05

  expect_equal(
    p_true_cond_2period(s3true = 1,
                        s2true = 1,
                        s1true = 1,
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    t12 + t21 + t11 - 2*t22)
  expect_equal(
    p_true_cond_2period(s3true = 0,
                        s2true = 1,
                        s1true = 1,
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    1 - (t12 + t21 + t11 - 2*t22))

  expect_equal(
    p_true_cond_2period(s3true = 1,
                        s2true = 1,
                        s1true = 0,
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    t12)
  expect_equal(
    p_true_cond_2period(s3true = 0,
                        s2true = 1,
                        s1true = 0,
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    1 - t12)

  expect_equal(
    p_true_cond_2period(s3true = 1,
                        s2true = 0,
                        s1true = 1,
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    t21)
  expect_equal(
    p_true_cond_2period(s3true = 0,
                        s2true = 0,
                        s1true = 1,
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    1 - t21)

  expect_equal(
    p_true_cond_2period(s3true = 1,
                        s2true = 0,
                        s1true = 0,
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    t22)
  expect_equal(
    p_true_cond_2period(s3true = 0,
                        s2true = 0,
                        s1true = 0,
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    1 - t22)

  # vectorize
  expect_equal(
    p_true_cond_2period(s3true = c(1, 0, 1, 0, 1, 0, 1, 0),
                        s2true = c(1, 1, 0, 0, 1, 1, 0, 0),
                        s1true = c(1, 1, 1, 1, 0, 0, 0, 0),
                        theta11 = t11,
                        theta12 = t12,
                        theta21 = t21,
                        theta22 = t22,
                        mu = 0),
    c(t12 + t21 + t11 - 2*t22,
      1 - (t12 + t21 + t11 - 2*t22),
      t21,
      1 - t21,
      t12,
      1 - t12,
      t22,
      1 - t22))


})






# TRUE 1 PERIOD CONDITIONAL PROBABILITY


test_that("p_true_cond_1period works", {
  t11 <- 0.12
  t12 <- 0.6
  t21 <- 0.35
  t22 <- 0.05
  mu  <- 0.25

  # 1,1
  p1 <- p_true_cond_1period(s2true = 1,
                             s1true = 1,
                             theta11 = t11,
                             theta12 = t12,
                             theta21 = t21,
                             theta22 = t22,
                             mu = mu)
  p2 <- p_true_cond_2period(s3true = 1,
                             s2true = 1,
                             s1true = 1,
                             theta11 = t11,
                             theta12 = t12,
                             theta21 = t21,
                             theta22 = t22,
                             mu = mu)*mu +
     p_true_cond_2period(s3true = 1,
                         s2true = 1,
                         s1true = 0,
                         theta11 = t11,
                         theta12 = t12,
                         theta21 = t21,
                         theta22 = t22,
                         mu = mu)*(1 - mu)
   expect_equal(p1,
                p2)

 # 1, 0
 p1 <- p_true_cond_1period(s2true = 1,
                           s1true = 0,
                           theta11 = t11,
                           theta12 = t12,
                           theta21 = t21,
                           theta22 = t22,
                           mu = mu)
 p2 <- p_true_cond_2period(s3true = 1,
                           s2true = 0,
                           s1true = 1,
                           theta11 = t11,
                           theta12 = t12,
                           theta21 = t21,
                           theta22 = t22,
                           mu = mu)*mu +
   p_true_cond_2period(s3true = 1,
                       s2true = 0,
                       s1true = 0,
                       theta11 = t11,
                       theta12 = t12,
                       theta21 = t21,
                       theta22 = t22,
                       mu = mu)*(1 - mu)
 expect_equal(p1,
              p2)

 # 0, 1
 p1 <- p_true_cond_1period(s2true = 0,
                           s1true = 1,
                           theta11 = t11,
                           theta12 = t12,
                           theta21 = t21,
                           theta22 = t22,
                           mu = mu)
 p2 <- p_true_cond_2period(s3true = 0,
                           s2true = 1,
                           s1true = 1,
                           theta11 = t11,
                           theta12 = t12,
                           theta21 = t21,
                           theta22 = t22,
                           mu = mu)*mu +
   p_true_cond_2period(s3true = 0,
                       s2true = 1,
                       s1true = 0,
                       theta11 = t11,
                       theta12 = t12,
                       theta21 = t21,
                       theta22 = t22,
                       mu = mu)*(1 - mu)
 expect_equal(p1,
              p2)

 # 0, 0
 p1 <- p_true_cond_1period(s2true = 0,
                           s1true = 0,
                           theta11 = t11,
                           theta12 = t12,
                           theta21 = t21,
                           theta22 = t22,
                           mu = mu)
 p2 <- p_true_cond_2period(s3true = 0,
                           s2true = 0,
                           s1true = 1,
                           theta11 = t11,
                           theta12 = t12,
                           theta21 = t21,
                           theta22 = t22,
                           mu = mu)*mu +
   p_true_cond_2period(s3true = 0,
                       s2true = 0,
                       s1true = 0,
                       theta11 = t11,
                       theta12 = t12,
                       theta21 = t21,
                       theta22 = t22,
                       mu = mu)*(1 - mu)
 expect_equal(p1,
              p2)

 # vectorize
 expect_equal(
   p_true_cond_1period(s2true  = c(1, 1, 0, 0),
                       s1true  = c(1, 0, 1, 0),
                       theta11 = t11,
                       theta12 = t12,
                       theta21 = t21,
                       theta22 = t22,
                       mu      = mu),
   c(p_true_cond_1period(s2true  = 1,
                         s1true  = 1,
                         theta11 = t11,
                         theta12 = t12,
                         theta21 = t21,
                         theta22 = t22,
                         mu      = mu),
     p_true_cond_1period(s2true  = 1,
                         s1true  = 0,
                         theta11 = t11,
                         theta12 = t12,
                         theta21 = t21,
                         theta22 = t22,
                         mu      = mu),
     p_true_cond_1period(s2true  = 0,
                         s1true  = 1,
                         theta11 = t11,
                         theta12 = t12,
                         theta21 = t21,
                         theta22 = t22,
                         mu      = mu),
     p_true_cond_1period(s2true  = 0,
                         s1true  = 0,
                         theta11 = t11,
                         theta12 = t12,
                         theta21 = t21,
                         theta22 = t22,
                         mu      = mu))
  )



})


test_that("true unconditional probability works", {

  mu    <- 0.25
  expect_equal(p_true_uncond(strue = 1,
                             mu = mu),
               mu)
  expect_equal(p_true_uncond(strue = 0,
                             mu = mu),
               mu^strue*(1 - mu)^(1 - strue))


  # vectorize
  expect_equal(
    p_true_uncond(strue = c(1, 0),
                  mu = mu),
    c(p_true_uncond(strue = 1,
                    mu = mu),
      p_true_uncond(strue = 0,
                    mu = mu))
  )

})

test_that("true joint 4wave works", {


  d <- expand.grid(s4true = c(1, 0),
                   s3true = c(1, 0),
                   s2true = c(1, 0),
                   s1true = c(1, 0))
  # 1, 1, 1, 1
  p1 <- p_true_joint_4wave(s4true = d$s4true,
                           s3true = d$s3true,
                           s2true = d$s2true,
                           s1true = d$s1true,
                           theta11 = t11,
                           theta12 = t12,
                           theta21 = t21,
                           theta22 = t22,
                           mu = mu)

  expect_equal(sum(p1),
               1)



})


test_that("observed joint prob 4wave works", {



  d <- expand.grid(s4 = c(1, 0),
                   s3 = c(1, 0),
                   s2 = c(1, 0),
                   s1 = c(1, 0),
                   s4true = c(1, 0),
                   s3true = c(1, 0),
                   s2true = c(1, 0),
                   s1true = c(1, 0))

  # 1, 1, 1, 1
  p1 <- p_obs_joint_4wave(s4 = d$s4,
                          s3 = d$s3,
                          s2 = d$s2,
                          s1 = d$s1,
                          s4true = d$s4true,
                          s3true = d$s3true,
                          s2true = d$s2true,
                          s1true = d$s1true,
                          theta11 = t11,
                          theta12 = t12,
                          theta21 = t21,
                          theta22 = t22,
                          err = 0.1,
                          mu = mu)

  expect_equal(sum(p1),
               1)



})







test_that("get_ar2_expected_log_likelihood works", {

  d <- expand.grid(s4 = c(1, 0),
                   s3 = c(1, 0),
                   s2 = c(1, 0),
                   s1 = c(1, 0))

  # vectorize
  expect_equal(
    get_ar2_expected_log_likelihood(s4 = c(1, 0),
                                    s3 = c(1, 1),
                                    s2 = c(1, 1),
                                    s1 = c(1, 0),
                                    theta11 = theta11,
                                    theta12 = theta12,
                                    theta21 = theta21,
                                    theta22 = theta22,
                                    err = err,
                                    mu = mu),
    c(get_ar2_expected_log_likelihood(s4 = 1,
                                     s3 = 1,
                                     s2 = 1,
                                     s1 = 1,
                                     theta11 = theta11,
                                     theta12 = theta12,
                                     theta21 = theta21,
                                     theta22 = theta22,
                                     err = err,
                                     mu = mu),
      get_ar2_expected_log_likelihood(s4 = 0,
                                     s3 = 1,
                                     s2 = 1,
                                     s1 = 0,
                                     theta11 = theta11,
                                     theta12 = theta12,
                                     theta21 = theta21,
                                     theta22 = theta22,
                                     err = err,
                                     mu = mu)) |>
      sum()
  )

})


