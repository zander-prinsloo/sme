get_ar1_tenure_individual_likelihood <- function(
    st1_observed, st2_observed, st3_observed,
    g_1, g_2, g_3,
    h_1, h_2, h_3,
    err,
    mu,
    theta_1,
    theta_2,
    sigma,
    lambda_h,
    lambda_g
){

  #_____________________________________________________________________________
  # Arguments-------------------------------------------------------------------
  if (lambda_h < 0) {
    cli::cli_abort("Exponential parameter `lambda_h` must be non-negative")
  }
  if (lambda_g < 0) {
    cli::cli_abort("Exponential parameter `lambda_g` must be non-negative")
  }
  if (sigma < 0) {
    cli::cli_abort("Gaussian parameter `sigma` must be non-negative")
  }

  #_____________________________________________________________________________
  # Vectorize-------------------------------------------------------------------
  get_ar1_sym_joint_probability_v <- Vectorize(
    FUN = get_ar1_sym_joint_probability,
    vectorize.args = c("st1_observed", "st2_observed", "st3_observed",
      "st1_true", "st2_true", "st3_true",
      "g_1", "g_2", "g_3",
      "h_1", "h_2", "h_3",
      "err",
      "mu",
      "theta_1",
      "theta_2",
      "sigma",
      "lambda_h",
      "lambda_g")
  )

  #_____________________________________________________________________________
  # Computations----------------------------------------------------------------
  lik <- get_ar1_sym_joint_probability_v(
    st1_observed = st1_observed,
    st2_observed = st2_observed,
    st3_observed = st3_observed,
    st1_true     = c(1, 1, 1, 1, 0, 0, 0, 0),
    st2_true     = c(1, 1, 0, 0, 1, 1, 0, 0),
    st3_true     = c(0, 1, 0, 1, 0, 1, 0, 1),
    g_1          = g_1,
    g_2          = g_2,
    g_3          = g_3,
    h_1          = h_1,
    h_2          = h_2,
    h_3          = h_3,
    err          = err,
    mu           = mu,
    theta_1      = theta_1,
    theta_2      = theta_2,
    sigma        = sigma,
    lambda_h     = lambda_h,
    lambda_g     = lambda_g
  )

  #_____________________________________________________________________________
  # Return----------------------------------------------------------------------
  lik

}





#' Joint probability - tenure and unemployment duration
#'
#' This function gives the full joint probability of a single individual
#' (i.e. non-vectorized) having the given observed potentially
#' misclassified status (s),
#' true status (s_star), observed potentially mismeasured or
#' misclassified tenure (g), observed potentially mismeasured or misclassified
#' unemployment duration. This probability is conditioned on the
#' parameters -- mu, theta_1, theta_2, sigma, lambda_h, lambda_g --
#' these specify the underlying distributions.
#'
#' @inheritParams get_ar1_sym_joint_probability
#' @param g_1 numeric: observed tenure t1
#' @param g_2 numeric: observed tenure t2
#' @param g_3 numeric: observed tenure t3
#' @param h_1 numeric: observed unemployment duration t1
#' @param h_2 numeric: observed unemployment duration t2
#' @param h_3 numeric: observed unemployment duration t3
#' @param sigma numeric: standard deviation parameter of normally
#' distributed measurement error term
#' @param lambda_h numeric: parameter of exponentially distributed
#' unemployment duration
#' @param lambda_g numeric: parameter of exponentially distributed
#' tenure
#'
#' @return numeric: joint probability
#' @export
get_ar1_tenure_joint_probability <- function(
    st1_observed, st2_observed, st3_observed,
    st1_true, st2_true, st3_true,
    g_1, g_2, g_3,
    h_1, h_2, h_3,
    err,
    mu,
    theta_1,
    theta_2,
    sigma,
    lambda_h,
    lambda_g
  ) {

  #_____________________________________________________________________________
  # Arguments-------------------------------------------------------------------
  p_err   <- pnorm(err)
  theta_1 <- pnorm(theta_1)
  theta_2 <- pnorm(theta_2)
  mu      <- pnorm(mu)
  if (lambda_h < 0) {
    cli::cli_abort("Exponential parameter `lambda_h` must be non-negative")
  }
  if (lambda_g < 0) {
    cli::cli_abort("Exponential parameter `lambda_g` must be non-negative")
  }
  if (sigma < 0) {
    cli::cli_abort("Gaussian parameter `sigma` must be non-negative")
  }

  #_____________________________________________________________________________
  # Misclassify Q1--------------------------------------------------------------
  p_misclass_1 <- p_err^(st1_observed * st1_true + (1 - st1_observed) * (1 - st1_true)) * # no error
    (1 - p_err)^((1 - st1_observed) * st1_true + st1_observed * (1 - st1_true))           # error

  #_____________________________________________________________________________
  # Misclassify Q2--------------------------------------------------------------
  p_misclass_2 <- p_err^(st2_observed * st2_true + (1 - st2_observed) * (1 - st2_true)) * # no error
    (1 - p_err)^((1 - st2_observed) * st2_true + st2_observed * (1 - st2_true))           # error

  #_____________________________________________________________________________
  # Misclassify Q3--------------------------------------------------------------
  p_misclass_3 <- p_err^(st3_observed * st3_true + (1 - st3_observed) * (1 - st3_true)) * # no error
    (1 - p_err)^((1 - st3_observed) * st3_true + st3_observed * (1 - st3_true))           # error

  #_____________________________________________________________________________
  # S1* in Q1 ----------------------------------------------------------------
  p_employed_q1 <- mu^st1_true * # empl
    (1 - mu)^(1 - st1_true)      # unempl

  #_____________________________________________________________________________
  # S2*|S1*---------------------------------------------------------------------
  p_trans_q2 <- theta_1^(st1_true * st2_true) *      # S*2 = 1 | S*1 = 1
    (1 - theta_1)^(st1_true * (1 - st2_true)) *      # S*2 = 0 | S*1 = 1
    theta_2^((1 - st1_true) * st2_true) *            # S*2 = 1 | S*1 = 0
    (1 - theta_2)^((1 - st1_true) * (1 - st2_true))  # S*2 = 0 | S*1 = 0

  #_____________________________________________________________________________
  # S3*|S2*---------------------------------------------------------------------
  p_trans_q3 <- theta_1^(st2_true * st3_true) *      # S*3 = 1 | S*2 = 1
    (1 - theta_1)^(st2_true * (1 - st3_true)) *      # S*3 = 0 | S*2 = 1
    theta_2^((1 - st2_true) * st3_true) *            # S*3 = 1 | S*2 = 0
    (1 - theta_2)^((1 - st2_true) * (1 - st3_true))  # S*3 = 0 | S*2 = 0

  #_____________________________________________________________________________
  # tenure Q1-------------------------------------------------------------------
  p_tenure_q1 <- ex_gaussian_density(g_1,
                                  sigma,
                                  lambda_g)^st1_observed           # s1 = 1

  #_____________________________________________________________________________
  # unemployment duration Q1----------------------------------------------------
  p_unempl_dur_q1 <- ex_gaussian_density(h_1,
                                      sigma,
                                      lambda_h)^(1 - st1_observed) # s1 = 0

  #_____________________________________________________________________________
  # tenure Q2-------------------------------------------------------------------
  p_tenure_q2 <-
    dnorm(
      x    = g_2 - g_1 - 0.25,
      mean = 0,
      sd   = sqrt(2 * sigma^2)
    )^(st2_observed*st1_observed*st1_true*st2_true) *         # S*2 = 1, S*1 = 1, S2 = 1, S1 = 1
    dnorm(
      x    = g_2 - 0.125,
      mean = 0,
      sd   = sigma
    )^(st2_observed*st2_true*(1 - st1_true)) *       # S*2 = 1, S*1 = 0, S2 = 1, S1 = 0/1
    ex_gaussian_density(
      g_2,
      sigma,
      lambda_g
    )^st2_observed*((1 - st2_true) +                 # S*2 = 0, S2 = 1
               st2_true*(1 - st1_observed)*st1_true) # S*2 = 1, S*1 = 1, S2 = 1, S1 = 0

  #_____________________________________________________________________________
  # unemployment duration Q2----------------------------------------------------
  p_unempl_dur_q2 <-
    dnorm(
      x    = h_2 - h_1 - 0.25,
      mean = 0,
      sd   = sqrt(2 * sigma^2)
    )^((1 - st2_observed)*(1 - st1_observed)*(1 - st1_true)*(1 - st2_true)) * # S*2=0, S*1=0, S2=0, S1=0
    dnorm(
      x    = h_2 - 0.125,
      mean = 0,
      sd   = sigma
    )^((1 - st2_observed) * (1 - st2_true) * st1_true) *            #
    ex_gaussian_density(
      h_2,
      sigma,
      lambda_h
    )^((1 - st2_observed)*(st2_true +
                    (1 - st2_true)*st1_observed*(1 - st1_true)))

  #_____________________________________________________________________________
  # tenure Q3-------------------------------------------------------------------
  p_tenure_q3 <-
    dnorm(
      x    = g_3 - g_1 - 0.5,
      mean = 0,
      sd   = sqrt(2*sigma^2)
    )^(st3_observed*st2_observed*st1_observed*st3_true*(1 - st2_true)*st1_true) *
    dnorm(
      x    = g_3 - g_2 - 0.25,
      mean = 0,
      sd   = sqrt(2 * sigma^2)
    )^(st3_observed*st2_observed*st2_true*st3_true) *
    dnorm(
      x    = g_3 - 0.125,
      mean = 0,
      sd   = sigma
    )^(st3_observed*st3_true*(1 - st2_true)) *
    ex_gaussian_density(
      g_3,
      sigma,
      lambda_g
    )^st3_observed*((1 - st3_true) + st3_true*(1 - st2_observed)*st2_true*((1 - st1_observed) + (1 - st1_true)))

  #_____________________________________________________________________________
  # unemployment duration Q3----------------------------------------------------
  p_unempl_dur_q3 <-
    dnorm(
      x    = h_3 - h_1 - 0.5,
      mean = 0,
      sd   = sqrt(2 * sigma^2)
    )^((1 - st3_observed) * (1 - st2_observed) * (1 - st1_observed) * (1 - st3_true) * st2_true * (1 - st1_true)) *
    dnorm(
      x    = h_3 - h_2 - 0.25,
      mean = 0,
      sd   = sqrt(2 * sigma^2)
    )^((1 - st3_observed) * (1 - st2_observed) * (1 - st2_true) * (1 - st3_true)) *
    dnorm(
      x    = h_3 - 0.125,
      mean = 0,
      sd   = sigma
    )^((1 - st3_observed) * (1 - st3_true) * st2_true) *
    ex_gaussian_density(
      h_3,
      sigma,
      lambda_h
    )^((1 - st3_observed) * (st3_true +
                      (1 - st3_true) * st2_observed * (1 - st2_true) * (st1_observed + st1_true)))


  #_____________________________________________________________________________
  # Likelihood------------------------------------------------------------------
  lik <- p_misclass_1 *
    p_misclass_2 *
    p_misclass_3 *
    p_employed_q1 *
    p_trans_q2 *
    p_trans_q3 *
    p_tenure_q1 *
    p_unempl_dur_q1 *
    p_tenure_q2 *
    p_unempl_dur_q2 *
    p_tenure_q3 *
    p_unempl_dur_q3

  #_____________________________________________________________________________
  # Return------------------------------------------------------------------
  lik

  #
  # # Compute the probability
  # probability <- err^(st1_observed * st1_true + (1 - st1_observed) * (1 - st1_true)) *
  #   (1 - err)^((1 - st1_observed) * st1_true + st1_observed * (1 - st1_true)) *
  #   err^(st2_observed * st2_true + (1 - st2_observed) * (1 - st2_true)) *
  #   (1 - err)^((1 - st2_observed) * st2_true + st2_observed * (1 - st2_true)) *
  #   err^(st3_observed * st3_true + (1 - st3_observed) * (1 - st3_true)) *
  #   (1 - err)^((1 - st3_observed) * st3_true + st3_observed * (1 - st3_true)) *
  #   mu^st1_true * (1 - mu)^(1 - st1_true) *
  #   theta_1^(st1_true * st2_true) * (1 - theta_1)^(st1_true * (1 - st2_true)) *
  #   theta_2^((1 - st1_true) * st2_true) * (1 - theta_2)^((1 - st1_true) * (1 - st2_true)) *
  #   theta_1^(st2_true * st3_true) * (1 - theta_1)^(st2_true * (1 - st3_true)) *
  #   theta_2^((1 - st2_true) * st3_true) * (1 - theta_2)^((1 - st2_true) * (1 - st3_true)) *
  #   ex_gaussian_density(g_1, sigma, lambda_g)^st1_observed *
  #   ex_gaussian_density(h_1, sigma, lambda_h)^(1 - st1_observed) *
  #   dnorm(g_2 - g_1 - 0.25, mean = 0, sd = sqrt(2 * sigma^2))^(st2_observed * st1_observed * st1_true * st2_true) *
  #   dnorm(g_2 - 0.125, mean = 0, sd = sigma)^(st2_observed * st2_true * (1 - st1_true)) *
  #   ex_gaussian_density(g_2, sigma, lambda_g)^st2_observed * ((1 - st2_true) + st2_true * (1 - st1_observed) * st1_true) *
  #   dnorm(h_2 - h_1 - 0.25, mean = 0, sd = sqrt(2 * sigma^2))^((1 - st2_observed) * (1 - st1_observed) * (1 - st1_true) * (1 - st2_true)) *
  #   dnorm(h_2 - 0.125, mean = 0, sd = sigma)^((1 - st2_observed) * (1 - st2_true) * st1_true) *
  #   ex_gaussian_density(h_2, sigma, lambda_h)^((1 - st2_observed) * (st2_true + (1 - st2_true) * st1_observed * (1 - st1_true))) *
  #   dnorm(g_3 - g_1 - 0.5, mean = 0, sd = sqrt(2 * sigma^2))^(st3_observed * st2_observed * st1_observed * st3_true * (1 - st2_true) * st1_true) *
  #   dnorm(g_3 - g_2 - 0.25, mean = 0, sd = sqrt(2 * sigma^2))^(st3_observed * st2_observed * st2_true * st3_true) *
  #   dnorm(g_3 - 0.125, mean = 0, sd = sigma)^(st3_observed * st3_true * (1 - st2_true)) *
  #   ex_gaussian_density(g_3, sigma, lambda_g)^st3_observed * ((1 - st3_true) + st3_true * (1 - st2_observed) * st2_true * ((1 - st1_observed) + (1 - st1_true))) *
  #   dnorm(h_3 - h_1 - 0.5, mean = 0, sd = sqrt(2 * sigma^2))^((1 - st3_observed) * (1 - st2_observed) * (1 - st1_observed) * (1 - st3_true) * st2_true * (1 - st1_true)) *
  #   dnorm(h_3 - h_2 - 0.25, mean = 0, sd = sqrt(2 * sigma^2))^((1 - st3_observed) * (1 - st2_observed) * (1 - st2_true) * (1 - st3_true)) *
  #   dnorm(h_3 - 0.125, mean = 0, sd = sigma)^((1 - st3_observed) * (1 - st3_true) * st2_true) *
  #   ex_gaussian_density(h_3, sigma, lambda_h)^((1 - st3_observed) * (st3_true + (1 - st3_true) * st2_observed * (1 - st2_true) * (st1_observed + st1_true)))
  #
  # Return the computed probability
  #return(probability)
}









#' Density of exGaussian distribution
#'
#' Sum of a gaussian and exponential random variable.
#' Former distributed N(0,sigma), latter distributed exp(lambda)
#'
#' @param x numeric value
#' @param sigma numeric: standard deviation for gaussian component
#' @param lambda numeric: parameter of exponential component
#'
#' @return numeric atomic vector length 1
#' @keywords internal
ex_gaussian_density <- function(x, sigma, lambda) {

  #_____________________________________________________________________________
  # Arguments-------------------------------------------------------------------
  if (lambda < 0) {
    cli::cli_abort("`lambda` parameter must be non-negative, as parameter of exponential distribution")
  }
  if (sigma < 0) {
    cli::cli_abort("`sigma` parameter must be non-negative, as parameter of gaussian distribution")
  }
  #_____________________________________________________________________________
  # Integrand-------------------------------------------------------------------
  integrand <- function(z) {
    dnorm(z, mean = 0, sd = sigma) * dexp(x - z, rate = lambda)
  }
  #_____________________________________________________________________________
  # Integrate--------------------------------------------------------------------
  result <- integrate(integrand, lower = -Inf, upper = x)

  #_____________________________________________________________________________
  # Return----------------------------------------------------------------------
  result$value

}




















