get_ar2_log_expected_likelihood <-
  function(s4,
           s3,
           s2,
           s1,
           theta11,
           theta12,
           theta21,
           theta22,
           err,
           mu) {

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
                     c(0, 0, 0, 0)), \(x) {
                       p_obs_joint_4wave(s4      = s4,
                                         s3      = s3,
                                         s2      = s2,
                                         s1      = s1,
                                         s4true  = x[1],
                                         s3true  = x[2],
                                         s2true  = x[3],
                                         s1true  = x[4],
                                         theta11 = theta11,
                                         theta12 = theta11,
                                         theta21 = theta21,
                                         theta22 = theta22,
                                         err     = err,
                                         mu      = mu)
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

  # misclass------------------------------------------------------
  number_error <- sum(c(s4 == s4true,
                        s3 == s3true,
                        s2 == s2true,
                        s1 == s1true))
  p_total_misclass <- p_num_misclass(err          = err,
                                     num_misclass = number_error,
                                     num_correct  = 4 - number_error)

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
#' @inheritParams get_ar2_log_expected_likelihood
#' @param mu Unconditional employment parameter
#'
#' @return
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

  prob <- if (s3true == 1) prob else 1 - prob

  prob

}



#' Likelihood of observing the `s` and `strue` given the probability of misclassification `err`
#'
#' @param s observed status
#' @param strue true status
#' @param err probability of misclassification
#'
#' @return
#' @export
p_misclass <- function(s, strue, err) {


  p <- apply(data.frame(s     = s,
                         strue = strue),
             MARGIN = 1,
             FUN = \(x) {
               print(x)
                     prob <- if (x[1] == x[2])
                       1 - err else
                         err
                     prob
  })

  p |> unname()

}

#' Probability of getting `num_misclass` misclassifications when
#' the probability of misclassification is `err`
#'
#' @param err probability of misclassification
#' @param num_misclass number of misclassifications
#' @param num_correct number of correct classifications: typicall 4 - num_misclass
#'
#' @return
#' @export
p_num_misclass <- function(err,
                           num_misclass,
                           num_correct = 4 - num_misclass) {

  prob <- dbinom(x    = num_misclass,
                 size = 4,
                 prob = err)
  prob

}

