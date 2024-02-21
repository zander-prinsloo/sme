
#' Log-likelihood of full observed sample - tenure/unemployment duration
#'
#' @inheritParams get_ar1_tenure_individual_likelihood
#' @inheritParams get_ar1_sym_full_likelihood
#'
#' @return numeric: log-likelihood
#' @export
get_ar1_tenure_full_likelihood <- function(
    st1_observed,
    st2_observed,
    st3_observed,
    g_1,
    g_2,
    g_3,
    h_1,
    h_2,
    h_3,
    theta_1,
    theta_2,
    sigma,
    lambda_g,
    lambda_h,
    err,
    #mu,
    skip_checks = FALSE
){
  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------
  if (!skip_checks) {
    if (
      length(unique(lengths(
        list(
          st1_observed,
          st2_observed,
          st3_observed,
          g_1, g_2, g_3,
          h_1, h_2, h_3
        )
      ))) > 1
    ) {
      cli::cli_abort(
        "The observed status (s_observed), tenure (g), and unemployment duration (h)
      vectors for all three periods should all have the same length."
      )
    }
  }

  #_____________________________________________________________________________
  # Computations ---------------------------------------------------------------
  lik <- get_ar1_tenure_individual_likelihood(
      st1_observed   = st1_observed,
      st2_observed   = st2_observed,
      st3_observed   = st3_observed,
      g_1            = g_1,
      g_2            = g_2,
      g_3            = g_3,
      h_1            = h_1,
      h_2            = h_2,
      h_3            = h_3,
      theta_1        = theta_1,
      theta_2        = theta_2,
      sigma          = sigma,
      lambda_g       = lambda_g,
      lambda_h       = lambda_h,
      err            = err,
      #mu             = mu,
      by_true_status = FALSE,
      skip_checks    = skip_checks
  )

  #_____________________________________________________________________________
  # Return --------------------------------------------------------------------
  lik

}


#' Likelihood vector - tenure and unemployment duration
#'
#' Vectorized version of [get_ar1_tenure_joint_probability], but integrating
#' out the true statuses. Gives the
#'
#' @inheritParams get_ar1_tenure_joint_probability
#' @param by_true_status logical: TRUE breaks the likelihood down by each of the
#' 8 combinations of true statuses over three periods while FALSE integrates over
#' the true statuses, giving a single likelihood for each element
#'
#' @return numeric data frame with 1 row if `by_true_status = FALSE` else 8 rows,
#' and columns equal to the length of the input vectors
#' @export
#'
#' @examples
#' get_ar1_tenure_individual_likelihood(
#' st1_observed = 1,
#' st2_observed = 1,
#' st3_observed = 1,
#' g_1          = 0.5,
#' g_2          = 0.75,
#' g_3          = 1.1,
#' h_1          = 0,
#' h_2          = 0,
#' h_3          = 0,
#' err          = 1.65,
#' mu           = 0.1,
#' theta_1      = 1.65,
#' theta_2      = -1.65,
#' sigma        = 0.1,
#' lambda_h     = 2,
#' lambda_g     = 2
#' )
get_ar1_tenure_individual_likelihood <- function(
    st1_observed,
    st2_observed,
    st3_observed,
    g_1,
    g_2,
    g_3,
    h_1,
    h_2,
    h_3,
    theta_1,
    theta_2,
    sigma,
    lambda_g,
    lambda_h,
    err,
    #mu,
    by_true_status = FALSE,
    skip_checks    = FALSE
){

  #_____________________________________________________________________________
  # Arguments-------------------------------------------------------------------
  if (!skip_checks) {
    if (lambda_h < 0) {
      cli::cli_abort("Exponential parameter `lambda_h` must be non-negative")
    }
    if (lambda_g < 0) {
      cli::cli_abort("Exponential parameter `lambda_g` must be non-negative")
    }
    if (sigma < 0) {
      cli::cli_alert_info("Gaussian parameter `sigma` must be non-negative - not incl exp transform")
    }
    if (
      length(unique(lengths(
        list(
          st1_observed,
          st2_observed,
          st3_observed,
          g_1, g_2, g_3,
          h_1, h_2, h_3
        )
      ))) > 1
    ) {
      cli::cli_alert_info(
        "Warning: The observed status (s_observed), tenure (g), and unemployment duration (h)
      vectors for all three periods should all have the same length."
      )
    }
    if (length(st1_observed) > 1) {
      cli::cli_alert_info(
        "Warning: inputs should be of length 1, otherwise use `get_ar1_tenure_individual_likelihood`"
      )
    }
  }

  #_____________________________________________________________________________
  # Calculation-----------------------------------------------------------------
  lik <- data.table::data.table(
    true_111 = get_ar1_tenure_joint_probability(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      st1_true     = 1,
      st2_true     = 1,
      st3_true     = 1,
      g_1          = g_1,
      g_2          = g_2,
      g_3          = g_3,
      h_1          = h_1,
      h_2          = h_2,
      h_3          = h_3,
      theta_1      = theta_1,
      theta_2      = theta_2,
      sigma        = sigma,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks
    ),
    true_110 = get_ar1_tenure_joint_probability(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      st1_true     = 1,
      st2_true     = 1,
      st3_true     = 0,
      g_1          = g_1,
      g_2          = g_2,
      g_3          = g_3,
      h_1          = h_1,
      h_2          = h_2,
      h_3          = h_3,
      theta_1      = theta_1,
      theta_2      = theta_2,
      sigma        = sigma,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks
    ),
    true_101 = get_ar1_tenure_joint_probability(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      st1_true     = 1,
      st2_true     = 0,
      st3_true     = 1,
      g_1          = g_1,
      g_2          = g_2,
      g_3          = g_3,
      h_1          = h_1,
      h_2          = h_2,
      h_3          = h_3,
      theta_1      = theta_1,
      theta_2      = theta_2,
      sigma        = sigma,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks
    ),
    true_100 = get_ar1_tenure_joint_probability(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      st1_true     = 1,
      st2_true     = 0,
      st3_true     = 0,
      g_1          = g_1,
      g_2          = g_2,
      g_3          = g_3,
      h_1          = h_1,
      h_2          = h_2,
      h_3          = h_3,
      theta_1      = theta_1,
      theta_2      = theta_2,
      sigma        = sigma,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks
    ),
    true_011 = get_ar1_tenure_joint_probability(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      st1_true     = 0,
      st2_true     = 1,
      st3_true     = 1,
      g_1          = g_1,
      g_2          = g_2,
      g_3          = g_3,
      h_1          = h_1,
      h_2          = h_2,
      h_3          = h_3,
      theta_1      = theta_1,
      theta_2      = theta_2,
      sigma        = sigma,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks
    ),
    true_010 = get_ar1_tenure_joint_probability(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      st1_true     = 0,
      st2_true     = 1,
      st3_true     = 0,
      g_1          = g_1,
      g_2          = g_2,
      g_3          = g_3,
      h_1          = h_1,
      h_2          = h_2,
      h_3          = h_3,
      theta_1      = theta_1,
      theta_2      = theta_2,
      sigma        = sigma,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks
    ),
    true_001 = get_ar1_tenure_joint_probability(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      st1_true     = 0,
      st2_true     = 0,
      st3_true     = 1,
      g_1          = g_1,
      g_2          = g_2,
      g_3          = g_3,
      h_1          = h_1,
      h_2          = h_2,
      h_3          = h_3,
      theta_1      = theta_1,
      theta_2      = theta_2,
      sigma        = sigma,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks
    ),
    true_000 = get_ar1_tenure_joint_probability(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      st1_true     = 0,
      st2_true     = 0,
      st3_true     = 0,
      g_1          = g_1,
      g_2          = g_2,
      g_3          = g_3,
      h_1          = h_1,
      h_2          = h_2,
      h_3          = h_3,
      theta_1      = theta_1,
      theta_2      = theta_2,
      sigma        = sigma,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks
    )
  )

  if (!by_true_status) {
    lik <- rowSums(lik)
  }

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
#'
#' @examples
#' get_ar1_tenure_joint_probability(
#' st1_observed = 1,
#' st2_observed = 1,
#' st3_observed = 1,
#' st1_true     = 1,
#' st2_true     = 1,
#' st3_true     = 1,
#' g_1          = 0.5,
#' g_2          = 0.75,
#' g_3          = 1.1,
#' h_1          = 0,
#' h_2          = 0,
#' h_3          = 0,
#' err          = 1.65,
#' mu           = 0.1,
#' theta_1      = 1.65,
#' theta_2      = -1.65,
#' sigma        = 0.1,
#' lambda_h     = 2,
#' lambda_g     = 2
#' )
get_ar1_tenure_joint_probability <- function(
    st1_observed,
    st2_observed,
    st3_observed,
    st1_true,
    st2_true,
    st3_true,
    g_1,
    g_2,
    g_3,
    h_1,
    h_2,
    h_3,
    theta_1,
    theta_2,
    sigma,
    lambda_g,
    lambda_h,
    err,
    #mu,
    skip_checks  = FALSE
  ) {

  #_____________________________________________________________________________
  # Arguments-------------------------------------------------------------------
  # p_err     <- pnorm(err)
  # p_theta_1 <- pnorm(theta_1)
  # p_theta_2 <- pnorm(theta_2)
   p_mu      <- (pnorm(theta_2))/(1 - pnorm(theta_1) + pnorm(theta_2))
  # p_sigma   <- exp(sigma)
  if (!skip_checks) {
    if (lambda_h < 0) {
      cli::cli_abort("Exponential parameter `lambda_h` must be non-negative")
    }
    if (lambda_g < 0) {
      cli::cli_abort("Exponential parameter `lambda_g` must be non-negative")
    }
    if (p_sigma < 0) {
      cli::cli_abort("Gaussian parameter `p_sigma` must be non-negative")
    }
  }

  #_____________________________________________________________________________
  # Misclassify Q1--------------------------------------------------------------
  p_misclass_1 <- pnorm(err)^(st1_observed * st1_true + (1 - st1_observed) * (1 - st1_true)) * # no error
    (1 - pnorm(err))^((1 - st1_observed) * st1_true + st1_observed * (1 - st1_true))           # error

  #_____________________________________________________________________________
  # Misclassify Q2--------------------------------------------------------------
  p_misclass_2 <- pnorm(err)^(st2_observed * st2_true + (1 - st2_observed) * (1 - st2_true)) * # no error
    (1 - pnorm(err))^((1 - st2_observed) * st2_true + st2_observed * (1 - st2_true))           # error

  #_____________________________________________________________________________
  # Misclassify Q3--------------------------------------------------------------
  p_misclass_3 <- pnorm(err)^(st3_observed * st3_true + (1 - st3_observed) * (1 - st3_true)) * # no error
    (1 - pnorm(err))^((1 - st3_observed) * st3_true + st3_observed * (1 - st3_true))           # error

  #_____________________________________________________________________________
  # S1* in Q1 ----------------------------------------------------------------
  p_employed_q1 <- p_mu^st1_true * # empl
    (1 - p_mu)^(1 - st1_true)      # unempl

  #_____________________________________________________________________________
  # S2*|S1*---------------------------------------------------------------------
  p_trans_q2 <- pnorm(theta_1)^(st1_true * st2_true) *      # S*2 = 1 | S*1 = 1
    (1 - pnorm(theta_1))^(st1_true * (1 - st2_true)) *      # S*2 = 0 | S*1 = 1
    pnorm(theta_2)^((1 - st1_true) * st2_true) *            # S*2 = 1 | S*1 = 0
    (1 - pnorm(theta_2))^((1 - st1_true) * (1 - st2_true))  # S*2 = 0 | S*1 = 0

  #_____________________________________________________________________________
  # S3*|S2*---------------------------------------------------------------------
  p_trans_q3 <- pnorm(theta_1)^(st2_true * st3_true) *      # S*3 = 1 | S*2 = 1
    (1 - pnorm(theta_1))^(st2_true * (1 - st3_true)) *      # S*3 = 0 | S*2 = 1
    pnorm(theta_2)^((1 - st2_true) * st3_true) *            # S*3 = 1 | S*2 = 0
    (1 - pnorm(theta_2))^((1 - st2_true) * (1 - st3_true))  # S*3 = 0 | S*2 = 0

  #_____________________________________________________________________________
  # tenure Q1-------------------------------------------------------------------
  p_tenure_q1 <- ex_gaussian_density(x      = g_1,
                                     sigma  = exp(sigma),
                                     lambda = lambda_g)^st1_observed           # s1 = 1

  #_____________________________________________________________________________
  # unemployment duration Q1----------------------------------------------------
  p_unempl_dur_q1 <- ex_gaussian_density(x      = h_1,
                                         sigma  = exp(sigma),
                                         lambda = lambda_h)^(1 - st1_observed) # s1 = 0

  #_____________________________________________________________________________
  # tenure Q2-------------------------------------------------------------------
  p_tenure_q2 <-
    dnorm(
      x      = g_2 - g_1 - 0.25,
      mean   = 0,
      sd     = sqrt(2 * exp(sigma)^2)
    )^(st2_observed*st1_observed*st1_true*st2_true) *  # S*2 = 1, S*1 = 1, S2 = 1, S1 = 1
    dnorm(
      x      = g_2 - 0.125,
      mean   = 0,
      sd     = exp(sigma)
    )^(st2_observed*st2_true*(1 - st1_true)) *         # S*2 = 1, S*1 = 0, S2 = 1, S1 = 0/1
    ex_gaussian_density(
      x      = g_2,
      sigma  = exp(sigma),
      lambda = lambda_g
    )^(st2_observed*(1 - st2_true)) *# S*2 = 0, S2 = 1
    ex_gaussian_density(
      x      = g_2 - 0.25,
      sigma  = exp(sigma),
      lambda = lambda_g
    )^(st2_observed*(st2_true*(1 - st1_observed)*st1_true))   # S*2 = 1, S*1 = 1, S2 = 1, S1 = 0
  #_____________________________________________________________________________
  # unemployment duration Q2----------------------------------------------------
  p_unempl_dur_q2 <-
    dnorm(
      x      = h_2 - h_1 - 0.25,
      mean   = 0,
      sd     = sqrt(2 * exp(sigma)^2)
    )^((1 - st2_observed)*(1 - st1_observed)*(1 - st1_true)*(1 - st2_true)) * # S*2=0, S*1=0, S2=0, S1=0
    dnorm(
      x      = h_2 - 0.125,
      mean   = 0,
      sd     = exp(sigma)
    )^( (1 - st2_observed)*(1 - st2_true)*st1_true ) *            #
    ex_gaussian_density(
      x      = h_2,
      sigma  = exp(sigma),
      lambda = lambda_h
    )^((1 - st2_observed)*(st2_true)) *            #
    ex_gaussian_density(
      x      = h_2 - 0.25,
      sigma  = exp(sigma),
      lambda = lambda_h
    )^((1 - st2_observed)*(1 - st2_true)*st1_observed*(1 - st1_true))

  #_____________________________________________________________________________
  # tenure Q3-------------------------------------------------------------------
  p_tenure_q3 <-
    dnorm(
      x    = g_3 - g_1 - 0.5,
      mean = 0,
      sd   = sqrt(2*exp(sigma)^2)
    )^(st3_observed*(1 - st2_observed)*st1_observed*st3_true*st2_true*st1_true) *
    dnorm(
      x    = g_3 - g_2 - 0.25,
      mean = 0,
      sd   = sqrt(2 * exp(sigma)^2)
    )^(st3_observed*st2_observed*st2_true*st3_true) *
    dnorm(
      x    = g_3 - 0.125,
      mean = 0,
      sd   = exp(sigma)
    )^(st3_observed*st3_true*(1 - st2_true)) *
    dnorm(
      x    = g_3 - 0.375,
      mean = 0,
      sd   = exp(sigma)
    )^( st3_observed*st3_true*(1 - st2_observed)*st2_true*(1 - st1_true) ) *
    ex_gaussian_density(
      x      = g_3,
      sigma  = exp(sigma),
      lambda = lambda_g
    )^(st3_observed*(1 - st3_true)) *
    ex_gaussian_density(
      x      = g_3 - 0.5,
      sigma  = exp(sigma),
      lambda = lambda_g
    )^(st3_observed*st3_true*(1 - st2_observed)*st2_true*(1 - st1_observed)*st1_true)

  #_____________________________________________________________________________
  # unemployment duration Q3----------------------------------------------------
  p_unempl_dur_q3 <-
    dnorm(
      x    = h_3 - h_1 - 0.5,
      mean = 0,
      sd   = sqrt(2 * exp(sigma)^2)
    )^((1 - st3_observed)*st2_observed*(1 - st1_observed)*(1 - st3_true)*(1 - st2_true)*(1 - st1_true)) *
    dnorm(
      x    = h_3 - h_2 - 0.25,
      mean = 0,
      sd   = sqrt(2 * exp(sigma)^2)
    )^((1 - st3_observed)*(1 - st2_observed)*(1 - st2_true)*(1 - st3_true)) *
    dnorm(
      x    = h_3 - 0.125,
      mean = 0,
      sd   = exp(sigma)
    )^((1 - st3_observed)*(1 - st3_true)*st2_true) *
    dnorm(
      x    = h_3 - 0.375,
      mean = 0,
      sd   = exp(sigma)
    )^( (1 - st3_observed)*(1 - st3_true)*st2_observed*(1 - st2_true)*st1_true ) *
    ex_gaussian_density(
      x      = h_3,
      sigma  = exp(sigma),
      lambda = lambda_h
    )^((1 - st3_observed)*st3_true) *
    ex_gaussian_density(
      x      = h_3 - 0.5,
      sigma  = exp(sigma),
      lambda = lambda_h
    )^((1 - st3_observed)*(1 - st3_true)*st2_observed*(1 - st2_true)*st1_observed*(1 - st1_true))


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
  lik <- pmax(lik, 1e-100)
  lik

}




#' Density of exGaussian distribution
#'
#' Sum of a gaussian and exponential random variable.
#' Former distributed N(0,sigma), latter distributed exp(lambda)
#'
#' @param x numeric value
#' @param sigma numeric: standard deviation for gaussian component
#' @param lambda numeric: parameter of exponential component - mean of exponential
#'
#' @return numeric atomic vector length 1
#' @keywords internal
ex_gaussian_density <- function(x, sigma, lambda) {

  #_____________________________________________________________________________
  # Arguments-------------------------------------------------------------------
   if (lambda <= 0) {
     cli::cli_alert_info("`lambda` parameter must be non-negative, as parameter of exponential distribution")
     lambda <- 1e-100
   }
   if (any(sigma <= 0)) {
     cli::cli_alert_info("`sigma` parameter must be postive, as parameter of gaussian distribution. It has been set to 1e-100")
     sigma <- 1e-100
     }
    #___________________________________________________________________________
    # Use package---------------------------------------------------------------
    result <- gamlss.dist::dexGAUS(
      x     = x,
      mu    = 0,
      sigma = sigma,
      nu    = lambda # mean = lambda
    )
    #___________________________________________________________________________
    # Return--------------------------------------------------------------------
    return(result)

}



















