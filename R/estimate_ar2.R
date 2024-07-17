#' Estimate the AR(2) misclassification error model
#'
#' @param st4_observed vector of observed binary status in the 4th wave - either logical or binary 0-1
#' @param st3_observed vector of observed binary status in the 3rd wave - either logical or binary 0-1
#' @param st2_observed vector of observed binary status in the 2nd wave - either logical or binary 0-1
#' @param st1_observed vector of observed binary status in the 1st wave - either logical or binary 0-1
#' @param weights vector of weights summing to population. NULL (default) receives equal weights
#' @param model_name name of the model
#' @param mu_end logical: should `mu` parameter be endogenous i.e. determined as function of thetas. Default is TRUE
#' @param init_params named double vector: initial parameters for the optimization
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
    mu_end               = TRUE,
    init_params          = if (isTRUE(mu_end)) {
      c("theta11" = 0.05,
        "theta12" = 0.95,
        "theta21" = 0.04,
        "theta22" = 0.03,
        "err"      = 1.96)
     } else{c(
      "theta11" = 0.059,
      "theta12" = 0.029,
      "theta21" = 0.029,
      "theta22" = 0.029,
      "err"     = 1.96,
      "mu"      = 0
    )},
    global_est           = FALSE,
    global_specification = NULL,
    verbose              = FALSE
){

  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------
  if (is.null(weights)) {
    weights <- rep(1, length(st1_observed))
  }

  if (isFALSE(mu_end)) {
    check_names <- c("theta11",
                     "theta12",
                     "theta21",
                     "theta22",
                     "err",
                     "mu")
  } else {
    check_names <- c("theta11",
                     "theta12",
                     "theta21",
                     "theta22",
                     "err")
  }


  if (is.null(init_params) & isFALSE(mu_end)) {
    init_params <- c(
      "theta11" = 0.059,
      "theta12" = 0.029,
      "theta21" = 0.029,
      "theta22" = 0.029,
      "err"      = 1.96,
      "mu"       = 0
    )
  } else if (all(!names(init_params) %chin% check_names)) {
    names(init_params) <- check_names
  } else if (is.null(init_params) & isTRUE(mu_end)) {
    init_params <- c(
      "theta11" = 0.059,
      "theta12" = 0.029,
      "theta21" = 0.029,
      "theta22" = 0.029,
      "err"      = 1.96
    )
  }
  #_____________________________________________________________________________
  # Obj Function----------------------------------------------------------------
  fn_ll <- function(par_vec){

    # Specify Parameters
    theta11 <- par_vec[1]
    theta12 <- par_vec[2]
    theta21 <- par_vec[3]
    theta22 <- par_vec[4]
    err     <- par_vec[5]
    err     <- 1 - pnorm(err)
    if (isFALSE(mu_end)) {
      mu      <- par_vec[6]
      mu      <- pnorm(mu)
    } else {
      mu <- NULL
    }


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
    print(log_lik)

    # Return
    return(log_lik)

  }

  #_____________________________________________________________________________
  # Global Estimation ----------------------------------------------------------
  if (isTRUE(global_est)) {

    if (is.null(global_specification)) {

      theta11_lower <- 0
      theta11_upper <- 1
      theta12_lower <- 0
      theta12_upper <- 1
      theta21_lower <- 0
      theta21_upper <- 1
      theta22_lower <- 0
      theta22_upper <- 1
      err_lower      <- -2.5
      err_upper      <- 2.5
      if (isFALSE(mu_end)) {
        mu_lower       <- -1
        mu_upper       <- 1
      }

    } else if (
      !is.list(global_specification) |
      !names(global_specification) %chin% c(
        "theta11_lower",
        "theta11_upper",
        "theta12_lower",
        "theta12_upper",
        "theta21_lower",
        "theta21_upper",
        "theta22_lower",
        "theta22_upper",
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
      theta11_lower <- global_specification$theta11_lower
      theta11_upper <- global_specification$theta11_upper
      theta12_lower <- global_specification$theta12_lower
      theta12_upper <- global_specification$theta12_upper
      theta21_lower <- global_specification$theta21_lower
      theta21_upper <- global_specification$theta21_upper
      theta22_lower <- global_specification$theta22_lower
      theta22_upper <- global_specification$theta22_upper
      err_lower     <- global_specification$err_lower
      err_upper     <- global_specification$err_upper
      if (isFALSE(mu_end)) {
        mu_lower      <- global_specification$mu_lower
        mu_upper      <- global_specification$mu_upper
      }
    }

    # Create Parameter Space
    if (isFALSE(mu_end)) {
      par_space <- makeParamSet(
        makeNumericParam(
          "theta11",
          lower = theta11_lower,
          upper = theta11_upper
        ),
        makeNumericParam(
          "theta12",
          lower = theta12_lower,
          upper = theta12_upper
        ),
        makeNumericParam(
          "theta21",
          lower = theta21_lower,
          upper = theta21_upper
        ),
        makeNumericParam(
          "theta22",
          lower = theta22_lower,
          upper = theta22_upper
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
    } else {
      par_space <- makeParamSet(
        makeNumericParam(
          "theta11",
          lower = theta11_lower,
          upper = theta11_upper
        ),
        makeNumericParam(
          "theta12",
          lower = theta12_lower,
          upper = theta12_upper
        ),
        makeNumericParam(
          "theta21",
          lower = theta21_lower,
          upper = theta21_upper
        ),
        makeNumericParam(
          "theta22",
          lower = theta22_lower,
          upper = theta22_upper
        ),
        makeNumericParam(
          "err",
          lower = err_lower,
          upper = err_upper
        )
      )
    }
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

    if (isFALSE(mu_end)) {
      df_initial_search <- df_initial_search |>
        rbind(
          c("theta11" = 0.06,
            "theta12" = 0.03,
            "theta21" = 0.95,
            "theta22" = 0.035,
            "err"      = 1.96,
            "mu"       = 0),
          c("theta11" = 0.138,
            "theta12" = 0.165,
            "theta21" = 0.771,
            "theta22" = 0.058,
            "err"      = 2.45,
            "mu"       = -0.05),
          init_params
        )
    } else {
      df_initial_search <- df_initial_search |>
        rbind(
          c("theta11" = 0.06,
            "theta12" = 0.03,
            "theta21" = 0.95,
            "theta22" = 0.035,
            "err"      = 1.96),
          c("theta11" = 0.138,
            "theta12" = 0.165,
            "theta21" = 0.771,
            "theta22" = 0.058,
            "err"      = 2.45),
          init_params
        )
    }

    print("hello2")
    print(list(df_initial_search, obj_func))
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

    if (verbose) {
      print(sme_estimation |> summary())
    }
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
      print(sme_estimation |> summary())
    }
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
    if (verbose) {
      print(sme_estimation |> summary())
    }

  }


  #_____________________________________________________________________________
  # Implied Probs --------------------------------------------------------------

  # Implied Probs ----
  true_111 <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "theta12 + theta21 + theta11 - 2*theta22",
    parameterNames = check_names
  )
  true_111_rate <- true_111[1][[1]]
  true_111_se   <- true_111[2][[1]]

  true_110 <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "theta12",
    parameterNames = check_names
  )
  true_110_rate <- true_110[1][[1]]
  true_110_se   <- true_110[2][[1]]

  true_101 <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "theta21",
    parameterNames = check_names
  )
  true_101_rate <- true_101[1][[1]]
  true_101_se   <- true_101[2][[1]]

  true_100 <- car::deltaMethod(
    object         = sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "theta22",
    parameterNames = check_names
  )
  true_100_rate <- true_100[1][[1]]
  true_100_se   <- true_100[2][[1]]

  if (isFALSE(mu_end)) {
    level <- car::deltaMethod(
      object         = sme_estimation,
      vcov.          = vcov(sme_estimation),
      g              = "pnorm(mu)",
      parameterNames = check_names
    )
    level_rate <- level[1][[1]]
    level_se   <- level[2][[1]]
  }
  # else {
  #   level <- car::deltaMethod(
  #     object         = sme_estimation,
  #     vcov.          = vcov(sme_estimation),
  #     g              = "mu_function(theta11, theta12, theta21, theta22)",
  #     parameterNames = check_names
  #   )
  # }



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
  if (isFALSE(mu_end)) {
    SME_CI_5 <- data.table(
      "Prob" = c("T_given_TT", "T_given_TF", "", "T_given_FT", "T_given_FF", "misclass"),
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
      "Prob" = c("T_given_TT", "T_given_TF", "", "T_given_FT", "T_given_FF", "misclass"),
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
      "Prob" = c("T_given_TT", "T_given_TF", "", "T_given_FT", "T_given_FF", "misclass"),
      "LowerBound" = c(true_111_rate - coef01*true_111_se,
                       true_110_rate - coef01*true_110_se,
                       true_101_rate - coef01*true_101_se,
                       true_100_rate - coef01*true_100_se,
                       misclass_rate + coef01*misclass_se),
      "UpperBound" = c(true_111_rate + coef01*true_111_se,
                       true_110_rate + coef01*true_110_se,
                       true_101_rate + coef01*true_101_se,
                       true_100_rate + coef01*true_100_se,
                       misclass_rate + coef01*misclass_se)
    )

  } else {
    SME_CI_5 <- data.table(
      "Prob" = c("T_given_TT", "T_given_TF", "T_given_FT", "T_given_FF", "misclass"),
      "LowerBound" = c(true_111_rate - coef5*true_111_se,
                       true_110_rate - coef5*true_110_se,
                       true_101_rate - coef5*true_101_se,
                       true_100_rate - coef5*true_100_se,
                       misclass_rate + coef5*misclass_se),
      "UpperBound" = c(true_111_rate + coef5*true_111_se,
                       true_110_rate + coef5*true_110_se,
                       true_101_rate + coef5*true_101_se,
                       true_100_rate + coef5*true_100_se,
                       misclass_rate + coef5*misclass_se)
    )
    SME_CI_1 <- data.table(
      "Prob" = c("T_given_TT", "T_given_TF", "T_given_FT", "T_given_FF", "misclass"),
      "LowerBound" = c(true_111_rate - coef1*true_111_se,
                       true_110_rate - coef1*true_110_se,
                       true_101_rate - coef1*true_101_se,
                       true_100_rate - coef1*true_100_se,
                       misclass_rate + coef1*misclass_se),
      "UpperBound" = c(true_111_rate + coef1*true_111_se,
                       true_110_rate + coef1*true_110_se,
                       true_101_rate + coef1*true_101_se,
                       true_100_rate + coef1*true_100_se,
                       misclass_rate + coef1*misclass_se)
    )
    SME_CI_01 <- data.table(
      "Prob" = c("T_given_TT", "T_given_TF", "T_given_FT", "T_given_FF", "misclass"),
      "LowerBound" = c(true_111_rate - coef01*true_111_se,
                       true_110_rate - coef01*true_110_se,
                       true_101_rate - coef01*true_101_se,
                       true_100_rate - coef01*true_100_se,
                       misclass_rate + coef01*misclass_se),
      "UpperBound" = c(true_111_rate + coef01*true_111_se,
                       true_110_rate + coef01*true_110_se,
                       true_101_rate + coef01*true_101_se,
                       true_100_rate + coef01*true_100_se,
                       misclass_rate + coef01*misclass_se)
    )

  }


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
    "misclass"           = misclass,
    "implied_CI_5perc"   = SME_CI_5,
    "implied_CI_1perc"   = SME_CI_1,
    "implied_CI_01perc"  = SME_CI_01
  )


  #_____________________________________________________________________________
  # Return -----------------------------------------------------------
  return(sme_results)


}
