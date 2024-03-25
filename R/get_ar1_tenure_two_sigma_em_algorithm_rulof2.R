
#' Log-likelihood of full observed sample - tenure/unemployment duration
#'
#' @inheritParams get_ar1_tenure_individual_likelihood
#' @inheritParams get_ar1_sym_full_likelihood
#'
#' @return numeric: log-likelihood
#' @export
get_ar1_tenure_full_likelihood_two_sigma <- function(
    st1_observed,
    st2_observed,
    st3_observed,
    g_1,
    g_2,
    g_3,
    h_1,
    h_2,
    h_3,
    posterior_probs,
    # theta_1,
    # theta_2,
    sigma_g,
    sigma_h,
    lambda_g,
    lambda_h,
    err,
    #mu,
    skip_checks = FALSE,
    return_posterior = FALSE,
    store_likelihoods = FALSE
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
  # print(lambda_g)
  # print(lambda_h)
# print(theta_1, " ", theta_2, " ", sigma_g, " ", sigma_h, " ", lambda_g, " ", lambda_h, " ", err)

  #_____________________________________________________________________________
  # Computations ---------------------------------------------------------------
  llik <- get_ar1_tenure_individual_likelihood_two_sigma(
    st1_observed   = st1_observed,
    st2_observed   = st2_observed,
    st3_observed   = st3_observed,
    g_1            = g_1,
    g_2            = g_2,
    g_3            = g_3,
    h_1            = h_1,
    h_2            = h_2,
    h_3            = h_3,
    posterior_probs = posterior_probs,
    # theta_1        = theta_1,
    # theta_2        = theta_2,
    sigma_g        = sigma_g,
    sigma_h        = sigma_h,
    lambda_g       = lambda_g,
    lambda_h       = lambda_h,
    err            = err,
    #mu             = mu,
    by_true_status = FALSE,
    skip_checks    = skip_checks,
    return_posterior = return_posterior,
    store_likelihoods = store_likelihoods
  )

  #_____________________________________________________________________________
  # Return --------------------------------------------------------------------
  llik

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
get_ar1_tenure_individual_likelihood_two_sigma <- function(
    st1_observed,
    st2_observed,
    st3_observed,
    g_1,
    g_2,
    g_3,
    h_1,
    h_2,
    h_3,
    posterior_probs,
    # theta_1,
    # theta_2,
    sigma_g,
    sigma_h,
    lambda_g,
    lambda_h,
    err,
    #mu,
    by_true_status = FALSE,
    skip_checks    = FALSE,
    return_posterior = return_posterior,
    store_likelihoods = store_likelihoods
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
    # if (sigma < 0) {
    #   cli::cli_alert_info("Gaussian parameter `sigma` must be non-negative - not incl exp transform")
    # }
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

  lambda_h  <- exp(lambda_h)
  lambda_g  <- exp(lambda_g)

  p_tenure_q2_case1 <- dnorm(
    x      = g_2 - g_1 - 0.25,
    mean   = 0,
    sd     = sqrt(2 * exp(sigma_g)^2)
  )

  p_tenure_q2_case2 <- dnorm(
    x      = g_2 - 0.125,
    mean   = 0,
    sd     = exp(sigma_g)
  )

  p_tenure_q2_case3 <- ex_gaussian_density(
    x      = g_2,
    sigma  = exp(sigma_g),
    lambda = lambda_g
  )

  p_tenure_q2_case4 <-  ex_gaussian_density(
    x      = g_2 - 0.25,
    sigma  = exp(sigma_g),
    lambda = lambda_g
  )

  p_unempl_dur_q2_case1 <- dnorm(
    x      = h_2 - h_1 - 0.25,
    mean   = 0,
    sd     = sqrt(2 * exp(sigma_h)^2)
  )

  p_unempl_dur_q2_case2 <- dnorm(
    x      = h_2 - 0.125,
    mean   = 0,
    sd     = exp(sigma_h)
  )

  p_unempl_dur_q2_case3 <- ex_gaussian_density(
    x      = h_2,
    sigma  = exp(sigma_h),
    lambda = lambda_h
  )

  p_unempl_dur_q2_case4 <- ex_gaussian_density(
    x      = h_2 - 0.25,
    sigma  = exp(sigma_h),
    lambda = lambda_h
  )

  p_tenure_q3_case1 <- dnorm(
    x    = g_3 - g_1 - 0.5,
    mean = 0,
    sd   = sqrt(2*exp(sigma_g)^2)
  )

  p_tenure_q3_case2 <- dnorm(
    x    = g_3 - g_2 - 0.25,
    mean = 0,
    sd   = sqrt(2 * exp(sigma_g)^2)
  )

  p_tenure_q3_case3 <- dnorm(
    x    = g_3 - 0.125,
    mean = 0,
    sd   = exp(sigma_g)
  )

  p_tenure_q3_case4 <- dnorm(
    x    = g_3 - 0.375,
    mean = 0,
    sd   = exp(sigma_g)
  )

  p_tenure_q3_case5 <- ex_gaussian_density(
    x      = g_3,
    sigma  = exp(sigma_g),
    lambda = lambda_g
  )

  p_tenure_q3_case6 <- ex_gaussian_density(
    x      = g_3 - 0.5,
    sigma  = exp(sigma_g),
    lambda = lambda_g
  )

  p_unempl_dur_q3_case1 <- dnorm(
    x    = h_3 - h_1 - 0.5,
    mean = 0,
    sd   = sqrt(2 * exp(sigma_h)^2)
  )

  p_unempl_dur_q3_case2 <- dnorm(
    x    = h_3 - h_2 - 0.25,
    mean = 0,
    sd   = sqrt(2 * exp(sigma_h)^2)
  )

  p_unempl_dur_q3_case3 <- dnorm(
    x    = h_3 - 0.125,
    mean = 0,
    sd   = exp(sigma_h)
  )

  p_unempl_dur_q3_case4 <- dnorm(
    x    = h_3 - 0.375,
    mean = 0,
    sd   = exp(sigma_h)
  )

  p_unempl_dur_q3_case5 <- ex_gaussian_density(
    x      = h_3,
    sigma  = exp(sigma_h),
    lambda = lambda_h
  )

  p_unempl_dur_q3_case6 <- ex_gaussian_density(
    x      = h_3 - 0.5,
    sigma  = exp(sigma_h),
    lambda = lambda_h
  )

    true_111 = get_ar1_tenure_joint_probability_two_sigma(
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
      # theta_1      = theta_1,
      # theta_2      = theta_2,
      sigma_g      = sigma_g,
      sigma_h      = sigma_h,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks,
      store_likelihoods = store_likelihoods
    )

    true_110 = get_ar1_tenure_joint_probability_two_sigma(
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
      # theta_1      = theta_1,
      # theta_2      = theta_2,
      sigma_g      = sigma_g,
      sigma_h      = sigma_h,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks,
      store_likelihoods = store_likelihoods
    )

    true_101 = get_ar1_tenure_joint_probability_two_sigma(
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
      # theta_1      = theta_1,
      # theta_2      = theta_2,
      sigma_g      = sigma_g,
      sigma_h      = sigma_h,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks,
      store_likelihoods = store_likelihoods
    )

    true_100 = get_ar1_tenure_joint_probability_two_sigma(
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
      # theta_1      = theta_1,
      # theta_2      = theta_2,
      sigma_g      = sigma_g,
      sigma_h      = sigma_h,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks,
      store_likelihoods = store_likelihoods
    )

    true_011 = get_ar1_tenure_joint_probability_two_sigma(
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
      # theta_1      = theta_1,
      # theta_2      = theta_2,
      sigma_g      = sigma_g,
      sigma_h      = sigma_h,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks,
      store_likelihoods = store_likelihoods
    )

    true_010 = get_ar1_tenure_joint_probability_two_sigma(
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
      # theta_1      = theta_1,
      # theta_2      = theta_2,
      sigma_g      = sigma_g,
      sigma_h      = sigma_h,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks,
      store_likelihoods = store_likelihoods
    )

    true_001 = get_ar1_tenure_joint_probability_two_sigma(
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
      # theta_1      = theta_1,
      # theta_2      = theta_2,
      sigma_g      = sigma_g,
      sigma_h      = sigma_h,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks,
      store_likelihoods = store_likelihoods
    )

    true_000 = get_ar1_tenure_joint_probability_two_sigma(
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
      # theta_1      = theta_1,
      # theta_2      = theta_2,
      sigma_g      = sigma_g,
      sigma_h      = sigma_h,
      lambda_g     = lambda_g,
      lambda_h     = lambda_h,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6,
      err          = err,
      #mu           = mu,
      skip_checks  = skip_checks,
      store_likelihoods = store_likelihoods
    )

if(store_likelihoods == TRUE) {

  llik <- data.table::data.table(
    true_111$llik,
    true_110$llik,
    true_101$llik,
    true_100$llik,
    true_011$llik,
    true_010$llik,
    true_001$llik,
    true_000$llik
  )

  empl_lik <- data.table::data.table(
    true_111$empl_lik,
    true_110$empl_lik,
    true_101$empl_lik,
    true_100$empl_lik,
    true_011$empl_lik,
    true_010$empl_lik,
    true_001$empl_lik,
    true_000$empl_lik
  )
# browser()
  duration_lik <- data.table::data.table(
    true_111$duration_lik,
    true_110$duration_lik,
    true_101$duration_lik,
    true_100$duration_lik,
    true_011$duration_lik,
    true_010$duration_lik,
    true_001$duration_lik,
    true_000$duration_lik
  )

  p_tenure_q1 <- data.table::data.table(
    true_111$p_tenure_q1,
    true_110$p_tenure_q1,
    true_101$p_tenure_q1,
    true_100$p_tenure_q1,
    true_011$p_tenure_q1,
    true_010$p_tenure_q1,
    true_001$p_tenure_q1,
    true_000$p_tenure_q1
  )

  p_tenure_q2 <- data.table::data.table(
    true_111$p_tenure_q2,
    true_110$p_tenure_q2,
    true_101$p_tenure_q2,
    true_100$p_tenure_q2,
    true_011$p_tenure_q2,
    true_010$p_tenure_q2,
    true_001$p_tenure_q2,
    true_000$p_tenure_q2
  )

  p_tenure_q3 <- data.table::data.table(
    true_111$p_tenure_q3,
    true_110$p_tenure_q3,
    true_101$p_tenure_q3,
    true_100$p_tenure_q3,
    true_011$p_tenure_q3,
    true_010$p_tenure_q3,
    true_001$p_tenure_q3,
    true_000$p_tenure_q3
  )

  p_unempl_dur_q1 <- data.table::data.table(
    true_111$p_unempl_dur_q1,
    true_110$p_unempl_dur_q1,
    true_101$p_unempl_dur_q1,
    true_100$p_unempl_dur_q1,
    true_011$p_unempl_dur_q1,
    true_010$p_unempl_dur_q1,
    true_001$p_unempl_dur_q1,
    true_000$p_unempl_dur_q1
  )

  p_unempl_dur_q2 <- data.table::data.table(
    true_111$p_unempl_dur_q2,
    true_110$p_unempl_dur_q2,
    true_101$p_unempl_dur_q2,
    true_100$p_unempl_dur_q2,
    true_011$p_unempl_dur_q2,
    true_010$p_unempl_dur_q2,
    true_001$p_unempl_dur_q2,
    true_000$p_unempl_dur_q2
  )

  p_unempl_dur_q3 <- data.table::data.table(
    true_111$p_unempl_dur_q3,
    true_110$p_unempl_dur_q3,
    true_101$p_unempl_dur_q3,
    true_100$p_unempl_dur_q3,
    true_011$p_unempl_dur_q3,
    true_010$p_unempl_dur_q3,
    true_001$p_unempl_dur_q3,
    true_000$p_unempl_dur_q3
  )

  duration_llik <- rowSums(duration_lik)
  empl_llik <- rowSums(empl_lik)
  lp_tenure_q1 <- rowSums(log(p_tenure_q1))
  lp_tenure_q2 <- rowSums(log(p_tenure_q2))
  lp_tenure_q3 <- rowSums(log(p_tenure_q3))
  lp_unempl_dur_q1 <- rowSums(log(p_unempl_dur_q1))
  lp_unempl_dur_q2 <- rowSums(log(p_unempl_dur_q2))
  lp_unempl_dur_q3 <- rowSums(log(p_unempl_dur_q3))
  # browser()
  print(paste0("Total: ", sum(llik)))
  print(paste0("Empl: ", sum(empl_llik)))
  print(paste0("Duration: ", sum(duration_llik)))
  print(paste0("Tenure1: ", sum(lp_tenure_q1)))
  print(paste0("Tenure2: ", sum(lp_tenure_q2)))
  print(paste0("Tenure3: ", sum(lp_tenure_q3)))
  print(paste0("Unempl_dur1: ", sum(lp_unempl_dur_q1)))
  print(paste0("Unempl_dur2: ", sum(lp_unempl_dur_q2)))
  print(paste0("Unempl_dur3: ", sum(lp_unempl_dur_q3)))

}

  #_____________________________________________________________________________
  # Return----------------------------------------------------------------------

  if(store_likelihoods == FALSE) {
    llik <- data.table::data.table(
      true_111,
      true_110,
      true_101,
      true_100,
      true_011,
      true_010,
      true_001,
      true_000
    )
  }

  if(return_posterior == T) {
    posterior_probs <- exp(llik) %>%
      mutate(tot = true_111 + true_110 + true_101 + true_100 + true_011 + true_010 + true_001 + true_000) %>%
      mutate(
        true_111 = true_111/tot,
        true_110 = true_110/tot,
        true_101 = true_101/tot,
        true_100 = true_100/tot,
        true_011 = true_011/tot,
        true_010 = true_010/tot,
        true_001 = true_001/tot,
        true_000 = true_000/tot
      ) %>%
      select(-tot)
    posterior_probs
  } else if(return_posterior == F) {

    weighted_llik = posterior_probs*llik

    if (!by_true_status) {
      llik <- rowSums(weighted_llik)
    }
    llik
  }

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
get_ar1_tenure_joint_probability_two_sigma <- function(
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
    # theta_1,
    # theta_2,
    sigma_g,
    sigma_h,
    lambda_g,
    lambda_h,
    p_tenure_q2_case1,
    p_tenure_q2_case2,
    p_tenure_q2_case3,
    p_tenure_q2_case4,
    p_unempl_dur_q2_case1,
    p_unempl_dur_q2_case2,
    p_unempl_dur_q2_case3,
    p_unempl_dur_q2_case4,
    p_tenure_q3_case1,
    p_tenure_q3_case2,
    p_tenure_q3_case3,
    p_tenure_q3_case4,
    p_tenure_q3_case5,
    p_tenure_q3_case6,
    p_unempl_dur_q3_case1,
    p_unempl_dur_q3_case2,
    p_unempl_dur_q3_case3,
    p_unempl_dur_q3_case4,
    p_unempl_dur_q3_case5,
    p_unempl_dur_q3_case6,
    err,
    #mu,
    skip_checks  = FALSE,
    store_likelihoods = FALSE

) {

  #_____________________________________________________________________________
  # Arguments-------------------------------------------------------------------
  # p_err     <- pnorm(err)
  # p_theta_1 <- pnorm(theta_1)
  # p_theta_2 <- pnorm(theta_2)
  # p_mu      <- (pnorm(theta_2))/(1 - pnorm(theta_1) + pnorm(theta_2))

  # p_sigma   <- exp(sigma)
  if (!skip_checks) {
    if (lambda_h < 0) {
      cli::cli_abort("Exponential parameter `lambda_h` must be non-negative")
    }
    if (lambda_g < 0) {
      cli::cli_abort("Exponential parameter `lambda_g` must be non-negative")
    }
    # if (p_sigma < 0) {
    #   cli::cli_abort("Gaussian parameter `p_sigma` must be non-negative")
    # }
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
  # p_employed_q1 <- p_mu^st1_true * # empl
  #   (1 - p_mu)^(1 - st1_true)      # unempl
  #
  # #_____________________________________________________________________________
  # # S2*|S1*---------------------------------------------------------------------
  # p_trans_q2 <- pnorm(theta_1)^(st1_true * st2_true) *      # S*2 = 1 | S*1 = 1
  #   (1 - pnorm(theta_1))^(st1_true * (1 - st2_true)) *      # S*2 = 0 | S*1 = 1
  #   pnorm(theta_2)^((1 - st1_true) * st2_true) *            # S*2 = 1 | S*1 = 0
  #   (1 - pnorm(theta_2))^((1 - st1_true) * (1 - st2_true))  # S*2 = 0 | S*1 = 0
  #
  # #_____________________________________________________________________________
  # # S3*|S2*---------------------------------------------------------------------
  # p_trans_q3 <- pnorm(theta_1)^(st2_true * st3_true) *      # S*3 = 1 | S*2 = 1
  #   (1 - pnorm(theta_1))^(st2_true * (1 - st3_true)) *      # S*3 = 0 | S*2 = 1
  #   pnorm(theta_2)^((1 - st2_true) * st3_true) *            # S*3 = 1 | S*2 = 0
  #   (1 - pnorm(theta_2))^((1 - st2_true) * (1 - st3_true))  # S*3 = 0 | S*2 = 0

  #_____________________________________________________________________________
  # tenure Q1-------------------------------------------------------------------
  p_tenure_q1 <- ex_gaussian_density(x      = g_1,
                                     sigma  = exp(sigma_g),
                                     lambda = lambda_g)^st1_observed           # s1 = 1
  # print(lambda_h)
  # print(exp(sigma_g))
  #_____________________________________________________________________________
  # unemployment duration Q1----------------------------------------------------
  p_unempl_dur_q1 <- ex_gaussian_density(x      = h_1,
                                         sigma  = exp(sigma_h),
                                         lambda = lambda_h)^(1 - st1_observed) # s1 = 0

  #_____________________________________________________________________________
  # tenure Q2-------------------------------------------------------------------

  # p_tenure_q2_case1 <- dnorm(
  #   x      = g_2 - g_1 - 0.25,
  #   mean   = 0,
  #   sd     = sqrt(2 * exp(sigma_g)^2)
  # )
  #
  # p_tenure_q2_case2 <- dnorm(
  #   x      = g_2 - 0.125,
  #   mean   = 0,
  #   sd     = exp(sigma_g)
  # )
  #
  # p_tenure_q2_case3 <- ex_gaussian_density(
  #   x      = g_2,
  #   sigma  = exp(sigma_g),
  #   lambda = lambda_g
  # )
  #
  # p_tenure_q2_case4 <-  ex_gaussian_density(
  #   x      = g_2 - 0.25,
  #   sigma  = exp(sigma_g),
  #   lambda = lambda_g
  # )

  p_tenure_q2 <-
    p_tenure_q2_case1^(st2_observed*st1_observed*st1_true*st2_true) *  # S*2 = 1, S*1 = 1, S2 = 1, S1 = 1
    p_tenure_q2_case2^(st2_observed*st2_true*(1 - st1_true)) *         # S*2 = 1, S*1 = 0, S2 = 1, S1 = 0/1
    p_tenure_q2_case3^(st2_observed*(1 - st2_true)) *# S*2 = 0, S2 = 1
    p_tenure_q2_case4^(st2_observed*(st2_true*(1 - st1_observed)*st1_true))   # S*2 = 1, S*1 = 1, S2 = 1, S1 = 0

  p_tenure_q2_case <-
    case_when(
      (st2_observed*st1_observed*st1_true*st2_true) == 1 ~ 1,
      (st2_observed*st2_true*(1 - st1_true)) == 1 ~ 2,
      (st2_observed*(1 - st2_true)) == 1 ~ 3,
      (st2_observed*(st2_true*(1 - st1_observed)*st1_true)) == 1 ~ 4,
      TRUE ~ 0
    )

  #_____________________________________________________________________________
  # unemployment duration Q2----------------------------------------------------

  # p_unempl_dur_q2_case1 <- dnorm(
  #   x      = h_2 - h_1 - 0.25,
  #   mean   = 0,
  #   sd     = sqrt(2 * exp(sigma_h)^2)
  # )
  #
  # p_unempl_dur_q2_case2 <- dnorm(
  #   x      = h_2 - 0.125,
  #   mean   = 0,
  #   sd     = exp(sigma_h)
  # )
  #
  # p_unempl_dur_q2_case3 <- ex_gaussian_density(
  #   x      = h_2,
  #   sigma  = exp(sigma_h),
  #   lambda = lambda_h
  # )
  #
  # p_unempl_dur_q2_case4 <- ex_gaussian_density(
  #   x      = h_2 - 0.25,
  #   sigma  = exp(sigma_h),
  #   lambda = lambda_h
  # )

  p_unempl_dur_q2 <-
    p_unempl_dur_q2_case1^((1 - st2_observed)*(1 - st1_observed)*(1 - st1_true)*(1 - st2_true)) * # S*2=0, S*1=0, S2=0, S1=0
    p_unempl_dur_q2_case2^((1 - st2_observed)*(1 - st2_true)*st1_true) *            #
    p_unempl_dur_q2_case3^((1 - st2_observed)*(st2_true)) *            #
    p_unempl_dur_q2_case4^((1 - st2_observed)*(1 - st2_true)*st1_observed*(1 - st1_true))

  p_unempl_dur_q2_case <-
    case_when(
      ((1 - st2_observed)*(1 - st1_observed)*(1 - st1_true)*(1 - st2_true)) == 1 ~ 1,
      ((1 - st2_observed)*(1 - st2_true)*st1_true) == 1 ~ 2,
      ((1 - st2_observed)*(st2_true)) == 1 ~ 3,
      ((1 - st2_observed)*(1 - st2_true)*st1_observed*(1 - st1_true)) == 1 ~ 4,
      TRUE ~ 0
    )

  #_____________________________________________________________________________
  # tenure Q3-------------------------------------------------------------------

  # p_tenure_q3_case1 <- dnorm(
  #   x    = g_3 - g_1 - 0.5,
  #   mean = 0,
  #   sd   = sqrt(2*exp(sigma_g)^2)
  # )
  #
  # p_tenure_q3_case2 <- dnorm(
  #   x    = g_3 - g_2 - 0.25,
  #   mean = 0,
  #   sd   = sqrt(2 * exp(sigma_g)^2)
  # )
  #
  # p_tenure_q3_case3 <- dnorm(
  #   x    = g_3 - 0.125,
  #   mean = 0,
  #   sd   = exp(sigma_g)
  # )
  #
  # p_tenure_q3_case4 <- dnorm(
  #   x    = g_3 - 0.375,
  #   mean = 0,
  #   sd   = exp(sigma_g)
  # )
  #
  # p_tenure_q3_case5 <- ex_gaussian_density(
  #   x      = g_3,
  #   sigma  = exp(sigma_g),
  #   lambda = lambda_g
  # )
  #
  # p_tenure_q3_case6 <- ex_gaussian_density(
  #   x      = g_3 - 0.5,
  #   sigma  = exp(sigma_g),
  #   lambda = lambda_g
  # )

  p_tenure_q3 <-
    p_tenure_q3_case1^(st3_observed*(1 - st2_observed)*st1_observed*st3_true*st2_true*st1_true) *
    p_tenure_q3_case2^(st3_observed*st2_observed*st2_true*st3_true) *
    p_tenure_q3_case3^(st3_observed*st3_true*(1 - st2_true)) *
    p_tenure_q3_case4^(st3_observed*st3_true*(1 - st2_observed)*st2_true*(1 - st1_true)) *
    p_tenure_q3_case5^(st3_observed*(1 - st3_true)) *
    p_tenure_q3_case6^(st3_observed*st3_true*(1 - st2_observed)*st2_true*(1 - st1_observed)*st1_true)

  p_tenure_q3_case <-
    case_when(
      (st3_observed*(1 - st2_observed)*st1_observed*st3_true*st2_true*st1_true) == 1 ~ 1,
      (st3_observed*st2_observed*st2_true*st3_true) == 1 ~ 2,
      (st3_observed*st3_true*(1 - st2_true)) == 1 ~ 3,
      (st3_observed*st3_true*(1 - st2_observed)*st2_true*(1 - st1_true)) == 1 ~ 4,
      (st3_observed*(1 - st3_true)) == 1 ~ 5,
      (st3_observed*st3_true*(1 - st2_observed)*st2_true*(1 - st1_observed)*st1_true) == 1 ~ 6,
      TRUE ~ 0
    )


  #_____________________________________________________________________________
  # unemployment duration Q3----------------------------------------------------

  # p_unempl_dur_q3_case1 <- dnorm(
  #   x    = h_3 - h_1 - 0.5,
  #   mean = 0,
  #   sd   = sqrt(2 * exp(sigma_h)^2)
  # )
  #
  # p_unempl_dur_q3_case2 <- dnorm(
  #   x    = h_3 - h_2 - 0.25,
  #   mean = 0,
  #   sd   = sqrt(2 * exp(sigma_h)^2)
  # )
  #
  # p_unempl_dur_q3_case3 <- dnorm(
  #   x    = h_3 - 0.125,
  #   mean = 0,
  #   sd   = exp(sigma_h)
  # )
  #
  # p_unempl_dur_q3_case4 <- dnorm(
  #   x    = h_3 - 0.375,
  #   mean = 0,
  #   sd   = exp(sigma_h)
  # )
  #
  # p_unempl_dur_q3_case5 <- ex_gaussian_density(
  #   x      = h_3,
  #   sigma  = exp(sigma_h),
  #   lambda = lambda_h
  # )
  #
  # p_unempl_dur_q3_case6 <- ex_gaussian_density(
  #   x      = h_3 - 0.5,
  #   sigma  = exp(sigma_h),
  #   lambda = lambda_h
  # )

  p_unempl_dur_q3 <-
    p_unempl_dur_q3_case1^((1 - st3_observed)*st2_observed*(1 - st1_observed)*(1 - st3_true)*(1 - st2_true)*(1 - st1_true)) *
    p_unempl_dur_q3_case2^((1 - st3_observed)*(1 - st2_observed)*(1 - st2_true)*(1 - st3_true)) *
    p_unempl_dur_q3_case3^((1 - st3_observed)*(1 - st3_true)*st2_true) *
    p_unempl_dur_q3_case4^((1 - st3_observed)*(1 - st3_true)*st2_observed*(1 - st2_true)*st1_true ) *
    p_unempl_dur_q3_case5^((1 - st3_observed)*st3_true) *
    p_unempl_dur_q3_case6^((1 - st3_observed)*(1 - st3_true)*st2_observed*(1 - st2_true)*st1_observed*(1 - st1_true))


  p_unempl_dur_q3_case <-
    case_when(
      ((1 - st3_observed)*st2_observed*(1 - st1_observed)*(1 - st3_true)*(1 - st2_true)*(1 - st1_true)) == 1 ~ 1,
      ((1 - st3_observed)*(1 - st2_observed)*(1 - st2_true)*(1 - st3_true)) == 1 ~ 2,
      ((1 - st3_observed)*(1 - st3_true)*st2_true) == 1 ~ 3,
      ((1 - st3_observed)*(1 - st3_true)*st2_observed*(1 - st2_true)*st1_true ) == 1 ~ 4,
      ((1 - st3_observed)*st3_true) == 1 ~ 5,
      ((1 - st3_observed)*(1 - st3_true)*st2_observed*(1 - st2_true)*st1_observed*(1 - st1_true)) == 1 ~ 6,
      TRUE ~ 0
    )

  p_misclass_1 <- ifelse(is.na(p_misclass_1) | p_misclass_1 < 1e-07, 1e-07, p_misclass_1)
  lp_misclass_1 <- log(p_misclass_1)


  p_misclass_2 <- ifelse(is.na(p_misclass_2) | p_misclass_2 < 1e-07, 1e-07, p_misclass_2)
  lp_misclass_2 <- log(p_misclass_2)

  p_misclass_3 <- ifelse(is.na(p_misclass_3) | p_misclass_3 < 1e-07, 1e-07, p_misclass_3)
  lp_misclass_3 <- log(p_misclass_3)

  # p_employed_q1 <- ifelse(is.na(p_employed_q1) | p_employed_q1 < 1e-07, 1e-07, p_employed_q1)
  # lp_employed_q1 <- log(p_employed_q1)
  #
  # p_trans_q2 <- ifelse(is.na(p_trans_q2) | p_trans_q2 < 1e-07, 1e-07, p_trans_q2)
  # lp_trans_q2 <- log(p_trans_q2)
  #
  # p_trans_q3 <- ifelse(is.na(p_trans_q3) | p_trans_q3 < 1e-07, 1e-07, p_trans_q3)
  # lp_trans_q3 <- log(p_trans_q3)

  p_tenure_q1 <- ifelse(is.na(p_tenure_q1) | p_tenure_q1 < 1e-07, 1e-07, p_tenure_q1)
  lp_tenure_q1 <- log(p_tenure_q1)

  p_tenure_q2 <- ifelse(is.na(p_tenure_q2) | p_tenure_q2 < 1e-07, 1e-07, p_tenure_q2)
  lp_tenure_q2 <- log(p_tenure_q2)

  p_tenure_q3 <- ifelse(is.na(p_tenure_q3) | p_tenure_q3 < 1e-07, 1e-07, p_tenure_q3)
  lp_tenure_q3 <- log(p_tenure_q3)

  p_unempl_dur_q1 <- ifelse(is.na(p_unempl_dur_q1) | p_unempl_dur_q1 < 1e-07, 1e-07, p_unempl_dur_q1)
  lp_unempl_dur_q1 <- log(p_unempl_dur_q1)

  p_unempl_dur_q2 <- ifelse(is.na(p_unempl_dur_q2) | p_unempl_dur_q2 < 1e-07, 1e-07, p_unempl_dur_q2)
  lp_unempl_dur_q2 <- log(p_unempl_dur_q2)

  p_unempl_dur_q3 <- ifelse(is.na(p_unempl_dur_q3) | p_unempl_dur_q3 < 1e-07, 1e-07, p_unempl_dur_q3)
  lp_unempl_dur_q3 <- log(p_unempl_dur_q3)

  #_____________________________________________________________________________
  # Likelihood------------------------------------------------------------------
  llik <- lp_misclass_1 +
    lp_misclass_2 +
    lp_misclass_3 +
    # lp_employed_q1 +
    # lp_trans_q2 +
    # lp_trans_q3 +
    lp_tenure_q1 +
    lp_unempl_dur_q1 +
    lp_tenure_q2 +
    lp_unempl_dur_q2 +
    lp_tenure_q3 +
    lp_unempl_dur_q3

  if(store_likelihoods == TRUE) {
    empl_lik <- lp_misclass_1 +
      lp_misclass_2 +
      lp_misclass_3
      # lp_employed_q1 +
      # lp_trans_q2 +
      # lp_trans_q3

    duration_lik <- lp_tenure_q1 +
      lp_unempl_dur_q1 +
      lp_tenure_q2 +
      lp_unempl_dur_q2 +
      lp_tenure_q3 +
      lp_unempl_dur_q3

    lik_all <- data.table::data.table(
      st1_observed = st1_observed,
      st2_observed = st2_observed,
      st3_observed = st3_observed,
      p_misclass_1 = p_misclass_1,
      p_misclass_2 = p_misclass_2,
      p_misclass_3 = p_misclass_3,
      # p_employed_q1 = p_employed_q1,
      # p_trans_q2 = p_trans_q2,
      # p_trans_q3 = p_trans_q3,
      p_tenure_q1 = p_tenure_q1,
      p_tenure_q2 = p_tenure_q2,
      p_tenure_q3 = p_tenure_q3,
      p_unempl_dur_q1 = p_unempl_dur_q1,
      p_unempl_dur_q2 = p_unempl_dur_q2,
      p_unempl_dur_q3 = p_unempl_dur_q3,
      p_tenure_q2_case = p_tenure_q2_case,
      p_tenure_q2_case1 = p_tenure_q2_case1,
      p_tenure_q2_case2 = p_tenure_q2_case2,
      p_tenure_q2_case3 = p_tenure_q2_case3,
      p_tenure_q2_case4 = p_tenure_q2_case4,
      p_unempl_dur_q2_case = p_unempl_dur_q2_case,
      p_unempl_dur_q2_case1 = p_unempl_dur_q2_case1,
      p_unempl_dur_q2_case2 = p_unempl_dur_q2_case2,
      p_unempl_dur_q2_case3 = p_unempl_dur_q2_case3,
      p_unempl_dur_q2_case4 = p_unempl_dur_q2_case4,
      p_tenure_q3_case = p_tenure_q3_case,
      p_tenure_q3_case1 = p_tenure_q3_case1,
      p_tenure_q3_case2 = p_tenure_q3_case2,
      p_tenure_q3_case3 = p_tenure_q3_case3,
      p_tenure_q3_case4 = p_tenure_q3_case4,
      p_tenure_q3_case5 = p_tenure_q3_case5,
      p_tenure_q3_case6 = p_tenure_q3_case6,
      p_unempl_dur_q3_case = p_unempl_dur_q3_case,
      p_unempl_dur_q3_case1 = p_unempl_dur_q3_case1,
      p_unempl_dur_q3_case2 = p_unempl_dur_q3_case2,
      p_unempl_dur_q3_case3 = p_unempl_dur_q3_case3,
      p_unempl_dur_q3_case4 = p_unempl_dur_q3_case4,
      p_unempl_dur_q3_case5 = p_unempl_dur_q3_case5,
      p_unempl_dur_q3_case6 = p_unempl_dur_q3_case6
    )

    saveRDS(lik_all, paste0("data/df_lik_all_empl_obs", st1_true, st2_true, st3_true, ".rds"))
    print(paste0("empl_obs", st1_true, st2_true, st3_true))
    list_lik <- list(llik = llik, empl_lik = empl_lik, duration_lik = duration_lik,
                     p_tenure_q1 = p_tenure_q1,
                     p_tenure_q2 = p_tenure_q2,
                     p_tenure_q3 = p_tenure_q3,
                     p_unempl_dur_q1 = p_unempl_dur_q1,
                     p_unempl_dur_q2 = p_unempl_dur_q2,
                     p_unempl_dur_q3 = p_unempl_dur_q3
    )
    list_lik

  } else if(store_likelihoods == FALSE) {
    llik
  }


  #_____________________________________________________________________________
  # Return------------------------------------------------------------------

  # browser()

  # print(paste("True employment Case: ", st1_true, st2_true, st3_true))
  #
  # nums <- which(lik == 0)
  #
  # if(sum(p_tenure_q1[nums] == 0) > 0) {
  #   print(paste("# obs with P(p_tenure_q1) = 0: ", sum(p_tenure_q1[nums] == 0)))
  #   nums_tenure_q1 <- which(p_tenure_q1 == 0)
  #   print("Problematic observed employment sequences:")
  #   print(unique(df_sim_ten[nums_tenure_q1, c(7,8,9)]))
  #   print("Problematic true employment sequences:")
  #   print(unique(df_sim_ten[nums_tenure_q1, c(1,2,3)]))
  # }
  #
  #
  # if(sum(p_tenure_q2[nums] == 0) > 0) {
  #   print(paste("# obs with P(p_tenure_q2) = 0: ", sum(p_tenure_q2[nums] == 0)))
  #   nums_tenure_q2 <- which(p_tenure_q2 == 0)
  #   print("Problematic observed employment sequences:")
  #   print(unique(df_sim_ten[nums_tenure_q2, c(7,8,9)]))
  #   print("Problematic true employment sequences:")
  #   print(unique(df_sim_ten[nums_tenure_q2, c(1,2,3)]))
  # }
  #
  # if(sum(p_tenure_q3[nums] == 0) > 0) {
  #   print(paste("# obs with P(p_tenure_q3) = 0: ", sum(p_tenure_q3[nums] == 0)))
  #   nums_tenure_q3 <- which(p_tenure_q3 == 0)
  #   print("Problematic observed employment sequences:")
  #   print(unique(df_sim_ten[nums_tenure_q3, c(7,8,9)]))
  #   print("Problematic true employment sequences:")
  #   print(unique(df_sim_ten[nums_tenure_q3, c(1,2,3)]))
  # }
  #
  #
  # if(sum(p_unempl_dur_q1[nums] == 0) > 0) {
  #   print(paste("# obs with P(p_unempl_dur_q1) = 0: ", sum(p_unempl_dur_q1[nums] == 0)))
  #   nums_unempl_dur_q1 <- which(p_unempl_dur_q1 == 0)
  #   print("Problematic observed employment sequences:")
  #   print(unique(df_sim_ten[nums_unempl_dur_q1, c(7,8,9)]))
  #   print("Problematic true employment sequences:")
  #   print(unique(df_sim_ten[nums_unempl_dur_q1, c(1,2,3)]))
  # }
  #
  #
  # if(sum(p_unempl_dur_q2[nums] == 0) > 0) {
  #   print(paste("# obs with P(p_unempl_dur_q2) = 0: ", sum(p_unempl_dur_q2[nums] == 0)))
  #   nums_unempl_dur_q2 <- which(p_unempl_dur_q2 == 0)
  #   print("Problematic observed employment sequences:")
  #   print(unique(df_sim_ten[nums_unempl_dur_q2, c(7,8,9)]))
  #   print("Problematic true employment sequences:")
  #   print(unique(df_sim_ten[nums_unempl_dur_q2, c(1,2,3)]))
  # }
  #
  # if(sum(p_unempl_dur_q3[nums] == 0) > 0) {
  #   print(paste("# obs with P(p_unempl_dur_q3) = 0: ", sum(p_unempl_dur_q3[nums] == 0)))
  #   nums_unempl_dur_q3 <- which(p_unempl_dur_q3 == 0)
  #   print("Problematic observed employment sequences:")
  #   print(unique(df_sim_ten[nums_unempl_dur_q3, c(7,8,9)]))
  #   print("Problematic true employment sequences:")
  #   print(unique(df_sim_ten[nums_unempl_dur_q3, c(1,2,3)]))
  # }

  # browser()
  # print(sum(p_tenure_q1[nums] == 0))
  # print(sum(p_tenure_q2[nums] == 0))
  # print(sum(p_tenure_q3[nums] == 0))
  #
  # print(sum(p_unempl_dur_q1[nums] == 0))
  # print(sum(p_unempl_dur_q2[nums] == 0))
  # print(sum(p_unempl_dur_q3[nums] == 0))


  # print(" ")

  # browser()
  # llik <- ifelse(is.na(llik)|is.infinite(llik), log(0.0000001), llik)
  # empl_lik <- ifelse(is.na(empl_lik), 0, empl_lik)
  # duration_lik <- ifelse(is.na(duration_lik), 0, duration_lik)
  # print(mean(p_misclass_1*p_misclass_2*p_misclass_3*p_employed_q1*p_trans_q2*p_trans_q3))
  # print(mean(p_tenure_q1*p_unempl_dur_q1*p_tenure_q2*p_unempl_dur_q2*p_tenure_q3*p_unempl_dur_q3))
  # print(sum(log(lik)))


  #


}




















