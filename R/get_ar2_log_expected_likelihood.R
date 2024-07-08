#' AR(2) misclassification expected log likelihood
#'
#' Used as the E-step in the EM algorithm.
#'
#' @inheritParams p_full_misclass
#' @inheritParams estimate_ar2
#'
#' @return numeric
#' @export
get_ar2_expected_log_likelihood <-
  function(s4,
           s3,
           s2,
           s1,
           weights = NULL,
           theta11,
           theta12,
           theta21,
           theta22,
           err,
           mu) {

  if (is.null(weights)) {
    weights <- rep(1, length(s4))
  }

  # expected likelihood
  lik <- lapply(list(c(1, 1, 1, 1),
                     c(0, 1, 1, 1),
                     c(1, 0, 1, 1),
                     c(0, 0, 1, 1),
                     c(1, 1, 0, 1),
                     c(0, 1, 0, 1),
                     c(1, 0, 0, 1),
                     c(0, 0, 0, 1),
                     c(1, 1, 1, 0),
                     c(0, 1, 1, 0),
                     c(1, 0, 1, 0),
                     c(0, 0, 1, 0),
                     c(1, 1, 0, 0),
                     c(0, 1, 0, 0),
                     c(1, 0, 0, 0),
                     c(0, 0, 0, 0)),
                \(x) {
                       loglik <-
                         p_true_joint_4wave(
                                  s4true  = x[1],
                                  s3true  = x[2],
                                  s2true  = x[3],
                                  s1true  = x[4],
                                  theta11 = theta11,
                                  theta12 = theta12,
                                  theta21 = theta21,
                                  theta22 = theta22,
                                  mu      = mu)

                       loglik <- loglik |>
                         log()
                       loglik <- loglik*weights

                       pweights <- p_full_misclass(s4 = s4,
                                                   s3 = s3,
                                                   s2 = s2,
                                                   s1 = s1,
                                                   s4true = x[1],
                                                   s3true = x[2],
                                                   s2true = x[3],
                                                   s1true = x[4],
                                                   err    = err)
                       ll <- loglik*pweights

                     })

  lik <- lik |>
    unlist() |>
    sum()

  lik

}



p_obs_joint_4wave <- function(s4,
                              s3,
                              s2,
                              s1,
                              s4true,
                              s3true,
                              s2true,
                              s1true,
                              theta11,
                              theta12,
                              theta21,
                              theta22,
                              err,
                              mu) {

  # misclass-----------------------------------------------------
  p_total_misclass <-
    p_full_misclass(s4 = s4,
                    s3 = s3,
                    s2 = s2,
                    s1 = s1,
                    s4true = s4true,
                    s3true = s3true,
                    s2true = s2true,
                    s1true = s1true,
                    err = err)

  # true
  p_true <- p_true_joint_4wave(s4true  = s4true,
                               s3true  = s3true,
                               s2true  = s2true,
                               s1true  = s1true,
                               theta11 = theta11,
                               theta12 = theta12,
                               theta21 = theta21,
                               theta22 = theta22,
                               mu      = mu)

  # total
  prob <- p_total_misclass*p_true

  prob

}



p_true_joint_4wave <- function(s4true,
                               s3true,
                               s2true,
                               s1true,
                               theta11,
                               theta12,
                               theta21,
                               theta22,
                               mu) {

  # Arguments


  # 1
  p1 <- p_true_uncond(strue  = s1true,
                      mu     = mu)

  # 1    --> 2
  p2 <- p_true_cond_1period(s2true  = s2true,
                            s1true  = s1true,
                            theta11 = theta11,
                            theta12 = theta12,
                            theta21 = theta21,
                            theta22 = theta22,
                            mu      = mu)

  # 1, 2 --> 3
  p3 <- p_true_cond_2period(s3true  = s3true,
                            s2true  = s2true,
                            s1true  = s1true,
                            theta11 = theta11,
                            theta12 = theta12,
                            theta21 = theta21,
                            theta22 = theta22,
                            mu      = mu)

  # 2, 3 --> 4
  p4 <- p_true_cond_2period(s3true  = s4true,
                            s2true  = s3true,
                            s1true  = s2true,
                            theta11 = theta11,
                            theta12 = theta12,
                            theta21 = theta21,
                            theta22 = theta22,
                            mu      = mu)

  prob <- p1*p2*p3*p4

  prob
}

p_true_uncond <- function(strue, mu) {

  prob <- mu^strue
  prob <- prob*(1 - mu)^(1 - strue)

  prob

}

p_true_cond_1period <- function(s2true,
                                s1true,
                                theta11,
                                theta12,
                                theta21,
                                theta22,
                                mu) {

  prob <- 2*s2true - 1
  prob <- prob*(mu*theta21 +
                  (1 - mu)*theta22 +
                  (mu*(theta11 - theta22) + theta12 - theta22)*s1true)
  prob <- prob + (1 - s2true)

  prob

}

#' Conditional probability of observing `s3true` given `s2true` and `s1true`
#' and parameters
#'
#' @inheritParams get_ar2_expected_log_likelihood
#' @param mu Unconditional employment parameter
#'
#' @return numeric
#' @export
p_true_cond_2period <- function(s3true,
                                s2true,
                                s1true,
                                theta11,
                                theta12,
                                theta21,
                                theta22,
                                mu) {

  prob <- theta22 +
    (theta12 - theta22)*s2true +
    (theta21 - theta22)*s1true +
    (theta11 - theta22)*s1true*s2true

  prob[s3true == 0] <- 1 - prob[s3true == 0]

  prob

}



#' Likelihood of observing the `s` and `strue` given the probability of misclassification `err`
#'
#' @param s observed status
#' @param strue true status
#' @param err probability of misclassification
#'
#' @return numeric
#' @export
p_misclass <- function(s, strue, err) {


  p <- apply(data.frame(s     = s,
                         strue = strue),
             MARGIN = 1,
             FUN = \(x) {
                     prob <- if (x[1] == x[2])
                       1 - err else
                         err
                     prob
  })

  p |>
    unname()

}

#' Probability of getting `num_misclass` misclassifications when
#' the probability of misclassification is `err`
#'
#' @param err probability of misclassification
#' @param num_misclass number of misclassifications
#' @param num_correct number of correct classifications: typicall 4 - num_misclass
#'
#' @return numeric
#' @export
p_num_misclass <- function(err,
                           num_misclass,
                           num_correct = 4 - num_misclass) {

  prob <- dbinom(x    = num_misclass,
                 size = num_misclass + num_correct,
                 prob = err)
  prob

}


#' Probablity of observing the `s` given the
#' true status `strue` and the
#' probability of misclassification `err`
#'
#' @param s4 Binary 0-1 **observed** status vector in the 4th wave
#' @param s3 Binary 0-1 **observed** status vector in the 4th wave
#' @param s2 Binary 0-1 **observed** status vector in the 4th wave
#' @param s1 Binary 0-1 **observed** status vector in the 4th wave
#' @param s4true Binary 0-1 **true** status vector in the 4th wave
#' @param s3true Binary 0-1 **true** status vector in the 4th wave
#' @param s2true Binary 0-1 **true** status vector in the 4th wave
#' @param s1true Binary 0-1 **true** status vector in the 4th wave
#' @param err Misclassification probability - scalar between 0 and 1
#'
#' @return vector of same length as `s4`
#' @export
p_full_misclass <- function(s4,
                            s3,
                            s2,
                            s1,
                            s4true,
                            s3true,
                            s2true,
                            s1true,
                            err) {


  p4 <- p_misclass(s     = s4,
                   strue = s4true,
                   err   = err)
  p3 <- p_misclass(s     = s3,
                   strue = s3true,
                   err   = err)
  p2 <- p_misclass(s     = s2,
                   strue = s2true,
                   err   = err)
  p1 <- p_misclass(s     = s1,
                   strue = s1true,
                   err   = err)

  p <- p1*p2*p3*p4

  p
}
