
#' Estimate AR1 tenure classification error estimator
#'
#' @inheritParams estimate_ar1_sym
#' @inheritParams get_ar1_tenure_full_likelihood
#' @param globa_n number of samples design - global estimation.
#' See `mlrMBO:::generateDesign`
#'
#' @return list
#' @export
estimate_ar1_tenure_two_sigma <- function(
    st3_observed         = NULL,
    st2_observed         = NULL,
    st1_observed         = NULL,
    gt3                  = NULL,
    gt2                  = NULL,
    gt1                  = NULL,
    ht3                  = NULL,
    ht2                  = NULL,
    ht1                  = NULL,
    posterior_probs      = NULL,
    weights              = NULL,
    model_name           = "SME_tenure",
    init_params          = NULL,
    global_est           = FALSE,
    global_specification = NULL,
    verbose              = FALSE,
    global_n             = 100,
    global_iterations    = 50
){

  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------
  if (is.null(weights)) {
    weights <- rep(1, length(st1_observed))
  }

  #_____________________________________________________________________________
  # Obj Function----------------------------------------------------------------
  fn_ll <- function(par_vec){

    # Specify Parameters
    # theta_1  <- par_vec[1]
    # theta_2  <- par_vec[2]
    sigma_g  <- par_vec[1]
    sigma_h  <- par_vec[2]
    lambda_g <- par_vec[3]
    lambda_h <- par_vec[4]
    err      <- par_vec[5]
    #mu       <- par_vec[7]

    #print(par_vec)
    log_lik <- get_ar1_tenure_full_likelihood_two_sigma(
      st1_observed   = st1_observed,
      st2_observed   = st2_observed,
      st3_observed   = st3_observed,
      g_1            = gt1,
      g_2            = gt2,
      g_3            = gt3,
      h_1            = ht1,
      h_2            = ht2,
      h_3            = ht3,
      posterior_probs = posterior_probs,
      # theta_1        = theta_1,
      # theta_2        = theta_2,
      sigma_g        = sigma_g,
      sigma_h        = sigma_h,
      lambda_g       = lambda_g,
      lambda_h       = lambda_h,
      err            = err,
      skip_checks    = TRUE
    )

    log_lik <- (log_lik*weights) |>
      collapse::fsum()

    # Return
    return(log_lik)

  }

  #_____________________________________________________________________________

  #_____________________________________________________________________________
  # Local Estimation -----------------------------------------------------------

  sme_estimation <- maxLik::maxBFGSR(
    fn_ll,
    start = init_params,
    print.level = 2,
    control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
  )

  if (verbose) {
    print(sme_estimation |> summary())
    print(sme_estimation$estimate |> names())
  }

# browser()
  #_____________________________________________________________________________
  # Results object -------------------------------------------------------------
  # sme_results <- list(
  #   # "model_type"         = model_name,
  #   "estimated_model"    = sme_estimation,
  #   # "mbo_model"          = global_model,
  #   "loglik"             = stats::logLik(sme_estimation),
  #   "model_summary"      = sme_estimation |> summary(),
  #   # "poverty_entry_rate" = SME_JobEntry_rate,
  #   # "poverty_exit_rate"  = SME_JobExit_rate,
  #   # #"poverty_rate"       = SME_Empl_rate,
  #   # "misclass_rate"      = SME_Misclass_rate,
  #   # "poverty_entry_se"   = SME_JobEntry_se,
  #   # "poverty_exit_se"    = SME_JobExit_se,
  #   # #"poverty_se"         = SME_Empl_se,
  #   # "misclass_se"        = SME_Misclass_se,
  #   # "implied_CI_5perc"   = SME_CI_5,
  #   # "implied_CI_1perc"   = SME_CI_1,
  #   # "implied_CI_01perc"  = SME_CI_01
  # )


  #_____________________________________________________________________________
  # Return -----------------------------------------------------------
  return(sme_estimation)


}

