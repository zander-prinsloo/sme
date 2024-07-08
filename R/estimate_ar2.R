#' Estimate the AR(2) misclassification error model
#'
#' @param st4_observed vector of observed binary status in the 4th wave - either logical or binary 0-1
#' @param st3_observed vector of observed binary status in the 3rd wave - either logical or binary 0-1
#' @param st2_observed vector of observed binary status in the 2nd wave - either logical or binary 0-1
#' @param st1_observed vector of observed binary status in the 1st wave - either logical or binary 0-1
#' @param weights vector of weights summing to population. NULL (default) receives equal weights
#' @param model_name name of the model
#' @param init_params initial parameters for the optimization
#' @param global_est logical indicating whether to estimate the model globally first, then
#' use results as starting values for the local optimization
#' @param global_specification list of bounds for the global optimization
#' @param verbose logical indicating whether to print information as the estimation goes
#'
#' @return list of estimated parameters and results
#' @export
estimate_ar2 <- function(
    st4_observed         = NULL,
    st3_observed         = NULL,
    st2_observed         = NULL,
    st1_observed         = NULL,
    weights              = NULL,
    model_name           = "SOME",
    init_params          = NULL,
    global_est           = FALSE,
    global_specification = NULL,
    verbose              = FALSE
){

  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------
  if (is.null(weights)) {
    weights <- rep(1, length(st1_observed))
  }

  check_names <- c("theta_11",
                   "theta_12",
                   "theta_21",
                   "theta_22",
                   "err",
                   "mu")

  if (is.null(init_params)) {
    init_params <- c(
      "theta_11" = 0.059,
      "theta_12" = 0.029,
      "theta_21" = 0.029,
      "theta_22" = 0.029,
      "err"      = 1.96,
      "mu"       = 0
    )
  } else if (all(!names(init_params) %chin% check_names)) {
    names(init_params) <- check_names
  }
  #_____________________________________________________________________________
  # Obj Function----------------------------------------------------------------
  fn_ll <- function(par_vec){

    # Specify Parameters
    theta_11 <- par_vec[1]
    theta_12 <- par_vec[2]
    theta_21 <- par_vec[3]
    theta_22 <- par_vec[4]
    err      <- par_vec[5]
    err      <- pnorm(err)
    mu       <- par_vec[6]
    mu       <- pnorm(mu)

    log_lik <- get_ar2_expected_log_likelihood(
      s4 = st4_observed,
      s3 = st3_observed,
      s2 = st2_observed,
      s1 = st1_observed,
      theta11 = theta11,
      theta12 = theta12,
      theta21 = theta21, # 0.713 + 0.058
      theta22 = theta22,         # 0.058
      mu      = mu,
      err     = err,
      weights = weights
    )

    log_lik <- (log_lik*weights)

    # Return
    return(log_lik)

  }

  #_____________________________________________________________________________
  # Global Estimation ----------------------------------------------------------
  if (isTRUE(global_est)) {

    if (is.null(global_specification)) {

      theta_11_lower <- 0
      theta_11_upper <- 1
      theta_12_lower <- 0
      theta_12_upper <- 1
      theta_21_lower <- 0
      theta_21_upper <- 1
      theta_22_lower <- 0
      theta_22_upper <- 1
      err_lower      <- -2.5
      err_upper      <- 2.5
      mu_lower       <- 0.5
      mu_upper       <- -0.5

    } else if (
      !is.list(global_specification) |
      !names(global_specification) %chin% c(
        "theta_11_lower",
        "theta_11_upper",
        "theta_12_lower",
        "theta_12_upper",
        "theta_21_lower",
        "theta_21_upper",
        "theta_22_lower",
        "theta_22_upper",
        "err_lower",
        "err_upper",
        "mu_lower",
        "mu_upper",
      )
    ) {
      cli::cli_abort(
        "If global specification is non NULL, then should be a list giving upper
        and lower bounds for each parameter"
      )
    } else {
      theta_11_lower <- global_specification$theta_11_lower
      theta_11_upper <- global_specification$theta_11_upper
      theta_12_lower <- global_specification$theta_12_lower
      theta_12_upper <- global_specification$theta_12_upper
      theta_21_lower <- global_specification$theta_21_lower
      theta_21_upper <- global_specification$theta_21_upper
      theta_22_lower <- global_specification$theta_22_lower
      theta_22_upper <- global_specification$theta_22_upper
      err_lower      <- global_specification$err_lower
      err_upper      <- global_specification$err_upper
      mu_lower       <- global_specification$mu_lower
      mu_upper       <- global_specification$mu_upper
    }

    # Create Parameter Space
    par_space <- makeParamSet(
      makeNumericParam(
        "theta_11",
        lower = theta_11_lower,
        upper = theta_11_upper
      ),
      makeNumericParam(
        "theta_12",
        lower = theta_12_lower,
        upper = theta_12_upper
      ),
      makeNumericParam(
        "theta_21",
        lower = theta_21_lower,
        upper = theta_21_upper
      ),
      makeNumericParam(
        "theta_22",
        lower = theta_22_lower,
        upper = theta_22_upper
      ),
      makeNumericParam(
        "err",
        lower = err_lower,
        upper = err_upper
      ),
      makeNumericParam(
        "mu",
        lower = mu_lower,
        upper = mu_upper
      )
    )

    obj_func <- smoof::makeSingleObjectiveFunction(
      fn       = fn_ll,
      minimize = FALSE,
      noisy    = TRUE,
      par.set  = par_space
    )
    sur_learner <- mlr::makeLearner(
      cl           = "regr.km",
      predict.type = "se",
      covtype      = "matern3_2",
      control      = list(trace = FALSE
      )
    )
    control_object <- mlrMBO::makeMBOControl() |>
      mlrMBO::setMBOControlTermination(iters = 50) |>
      mlrMBO::setMBOControlInfill(crit = mlrMBO::makeMBOInfillCritEI())

    # Define initial search ----
    set.seed(1234)
    df_initial_search <- generateDesign(
      n       = 100,
      par.set = par_space,
      fun     = lhs::randomLHS
    )

    df_initial_search <- df_initial_search |>
      rbind(
        c("theta_11" = 0.06,
          "theta_12" = 0.03,
          "theta_21" = 0.95,
          "theta_22" = 0.035,
          "err"      = 1.96,
          "mu"       = 0),
        c("theta_11" = 0.138,
          "theta_12" = 0.165,
          "theta_21" = 0.771,
          "theta_22" = 0.058,
          "err"      = 2.45,
          "mu"       = -0.05),
        init_params
      )
    df_initial_search$y <- apply(
      df_initial_search,
      1,
      obj_func
    )
    print(df_initial_search)

    # Do Global estimation ----
    global_model <- mlrMBO::mbo(
      fun       = obj_func,
      learner   = sur_learner,
      control   = control_object,
      design    = df_initial_search,
      show.info = verbose
    )

    init_params <- global_model$x |>
      unlist()

  } else {
    global_model <- NULL
  }

  #_____________________________________________________________________________
  # Local Estimation -----------------------------------------------------------

  sme_estimation <- maxLik(
    fn_ll,
    start = init_params,
    method = "NM"
  )
  if (verbose) {
    print(sme_estimation |> summary())
  }
  if (maxLik::returnCode(sme_estimation) == 1) {
    if (verbose) {
      cli::cli_alert_info("Re-estimating local by refreshing initial params")
    }
    sme_estimation <- maxLik(
      fn_ll,
      start = sme_estimation$estimate,
      method = "NM"
    )
  }
  if (maxLik::returnCode(sme_estimation) == 1) {
    if (verbose) {
      cli::cli_alert_info("Re-estimating local by refreshing initial params")
    }
    sme_estimation <- maxLik(
      fn_ll,
      start = sme_estimation$estimate,
      method = "NM"
    )
  }
  if (maxLik::returnCode(sme_estimation) == 1) {
    if (verbose) {
      cli::cli_alert_info("Re-estimating local by refreshing initial params")
    }
    sme_estimation <- maxLik(
      fn_ll,
      start = sme_estimation$estimate,
      method = "NM"
    )

    if (verbose) {
      if (maxLik::returnCode(sme_estimation) == 1) {
        cli::cli_alert_info(
          "User should use output params of local estimation
           as the initial params and rerun the estimation, but do not
          include global estimation"
        )
      }
    }

  }

  #_____________________________________________________________________________
  # Implied Probs --------------------------------------------------------------

  # Implied Probs ----
  true_111 <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "theta_12 + theta_21 + theta-11 - 2*theta_22",
    parameterNames = check_names
  )
  true_111_rate <- true_111[1][[1]]
  true_111_se   <- true_111[2][[1]]

  true_110 <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "theta_12",
    parameterNames = check_names
  )
  true_110_rate <- true_110[1][[1]]
  true_110_se   <- true_110[2][[1]]

  true_101 <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "theta_21",
    parameterNames = check_names
  )
  true_101_rate <- true_101[1][[1]]
  true_101_se   <- true_101[2][[1]]

  true_100 <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "theta_22",
    parameterNames = check_names
  )
  true_100_rate <- true_100[1][[1]]
  true_100_se   <- true_100[2][[1]]

  true_level <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "pnorm(mu)",
    parameterNames = check_names
  )
  true_level_rate <- true_level[1][[1]]
  true_level_se   <- true_level[2][[1]]

  misclass <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "1 - pnorm(err)",
    parameterNames = check_names
  )
  misclass_rate <- misclass[1][[1]]
  misclass_se   <- misclass[2][[1]]


  # Implied CIs ----
  coef5  <- qnorm(1 - 0.05/2)
  coef1  <- qnorm(1 - 0.01/2)
  coef01 <- qnorm(1 - 0.001/2)
  SME_CI_5 <- data.table(
    "Prob" = c("T_given_TT", "T_given_TF", "", "T_given_FT", "T_given_FF", "level", "misclass"),
    "LowerBound" = c(true_111_rate - coef5*true_111_se,
                     true_110_rate - coef5*true_110_se,
                     true_101_rate - coef5*true_101_se,
                     true_100_rate - coef5*true_100_se,
                     level_rate    - coef5*level_se,
                     misclass_rate + coef5*misclass_se),
    "UpperBound" = c(true_111_rate + coef5*true_111_se,
                     true_110_rate + coef5*true_110_se,
                     true_101_rate + coef5*true_101_se,
                     true_100_rate + coef5*true_100_se,
                     level_rate    + coef5*level_se,
                     misclass_rate + coef5*misclass_se)
  )
  SME_CI_1 <- data.table(
    "Prob" = c("T_given_TT", "T_given_TF", "", "T_given_FT", "T_given_FF", "level", "misclass"),
    "LowerBound" = c(true_111_rate - coef1*true_111_se,
                     true_110_rate - coef1*true_110_se,
                     true_101_rate - coef1*true_101_se,
                     true_100_rate - coef1*true_100_se,
                     level_rate    - coef1*level_se,
                     misclass_rate + coef1*misclass_se),
    "UpperBound" = c(true_111_rate + coef1*true_111_se,
                     true_110_rate + coef1*true_110_se,
                     true_101_rate + coef1*true_101_se,
                     true_100_rate + coef1*true_100_se,
                     level_rate    + coef1*level_se,
                     misclass_rate + coef1*misclass_se)
  )
  SME_CI_01 <- data.table(
    "Prob" = c("T_given_TT", "T_given_TF", "", "T_given_FT", "T_given_FF", "level", "misclass"),
    "LowerBound" = c(true_111_rate - coef01*true_111_se,
                     true_110_rate - coef01*true_110_se,
                     true_101_rate - coef01*true_101_se,
                     true_100_rate - coef01*true_100_se,
                     level_rate    - coef01*level_se,
                     misclass_rate + coef01*misclass_se),
    "UpperBound" = c(true_111_rate + coef01*true_111_se,
                     true_110_rate + coef01*true_110_se,
                     true_101_rate + coef01*true_101_se,
                     true_100_rate + coef01*true_100_se,
                     level_rate    + coef01*level_se,
                     misclass_rate + coef01*misclass_se)
  )


  #_____________________________________________________________________________
  # Results object -------------------------------------------------------------
  sme_results <- list(
    "lags"               = "AR2",
    "model_type"         = model_name,
    "estimated_model"    = sme_estimation,
    "mbo_model"          = global_model,
    "loglik"             = logLik(sme_estimation),
    "model_summary"      = sme_estimation %>% summary(),
    "true_111"           = true_111,
    "true_110"           = true_110,
    "true_101"           = true_101,
    "true_100"           = true_100,
    "level"              = level,
    "misclass"           = misclass,
    "implied_CI_5perc"   = SME_CI_5,
    "implied_CI_1perc"   = SME_CI_1,
    "implied_CI_01perc"  = SME_CI_01
  )


  #_____________________________________________________________________________
  # Return -----------------------------------------------------------
  return(sme_results)


}
