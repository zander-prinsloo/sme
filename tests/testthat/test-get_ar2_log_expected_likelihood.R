


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



})

