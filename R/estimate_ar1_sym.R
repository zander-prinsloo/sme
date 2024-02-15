
#' Estimate AR1 symmetric classification error estimator
#'
#' @param st3_observed numeric vector: observed status in t3
#' @param st2_observed numeric vector: observed status in t2
#' @param st1_observed numeric vector: observed status in t1
#' @param weights numeric: survey weights. Default is NULL, giving
#' equal weight to each observation.
#' @param model_name character: used to identify different estimated results
#' @param init_params named atomic vector: initial values for parameters
#' @param global_est logical: should a global estimator be used to determine
#' the starting values for local estimation
#' @param global_specification list: parameter bounds for global estimation
#' @param verbose logical: should estimation info be printed to console
#'
#' @return list
#' @export
estimate_ar1_sym <- function(
    st3_observed         = NULL,
    st2_observed         = NULL,
    st1_observed         = NULL,
    weights              = NULL,
    model_name           = "SME",
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
  if (is.null(init_params)) {
    init_params <- c(
      "theta_01" = 1.5,
      "theta_02" = -1.5,
      "err"      = 1.96,
      "mu"       = 0.5
    )
  } else if (all(!names(init_params) %chin% c("theta_01", "theta_02", "err", "mu"))) {
    names(init_params) <- c("theta_01", "theta_02", "err", "mu")
  }
  #_____________________________________________________________________________
  # Obj Function----------------------------------------------------------------
  fn_ll <- function(par_vec){

    # Specify Parameters
    theta_01 <- par_vec[1]
    theta_02 <- par_vec[2]
    err      <- par_vec[3]
    mu       <- par_vec[4]

    log_lik <- get_ar1_sym_full_likelihood(
      st3_observed = st3_observed,
      st2_observed = st2_observed,
      st1_observed = st1_observed,
      par_vec      = c(
        theta_01,
        theta_02,
        err,
        mu
      )
      # theta_01     = theta_01,
      # theta_02     = theta_02,
      # err          = err,
      # mu           = mu
    ) |>
      log()

    log_lik <- (log_lik*weights) |>
      fsum()

    # Return
    return(log_lik)

  }

  #_____________________________________________________________________________
  # Global Estimation ----------------------------------------------------------
  if (isTRUE(global_est)) {

    if (is.null(global_specification)) {

      theta_01_lower <- -3
      theta_01_upper <-  3
      theta_02_lower <- -3
      theta_02_upper <-  3
      err_lower      <- -3
      err_upper      <-  3
      mu_lower       <- -3
      mu_upper       <-  3

    } else if (
      !is.list(global_specification) |
      !names(global_specification) %chin% c(
        "theta_01_lower",
        "theta_01_upper",
        "theta_02_lower",
        "theta_02_upper",
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
      theta_01_lower <- global_specification$theta_01_lower
      theta_01_upper <- global_specification$theta_01_upper
      theta_02_lower <- global_specification$theta_02_lower
      theta_02_upper <- global_specification$theta_02_upper
      err_lower      <- global_specification$err_lower
      err_upper      <- global_specification$err_upper
      mu_lower       <- global_specification$mu_lower
      mu_upper       <- global_specification$mu_upper
    }

    # Create Parameter Space
    par_space <- makeParamSet(
      makeNumericParam(
        "theta_01",
        lower = theta_01_lower,
        upper = theta_01_upper
      ),
      makeNumericParam(
        "theta_02",
        lower = theta_02_lower,
        upper = theta_02_upper
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
        c("theta_01" = 0.957,
          "theta_02" = -1.465,
          "err"       = 1.929,
          "mu"       = -0.0038),
        c(
          "theta_01" = 1.199927414,
          "theta_02" = -2,
          "err"       = 2.9,
          "mu"       = -0.018828444
        ),
        c(
          "theta_01" = 1.98,
          "theta_02" = -1.89,
          "err"       = 1.86,
          "mu"       = -0.02
        ),
        c(
          "theta_01" = 1.389,
          "theta_02" = -1.364,
          "err"       = 3,
          "mu"       = -0.02
        ),
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
  # Implied Probs ===-----------------------------------------------------------

  # Implied Probs ----
  SME_JobExit <- car::deltaMethod(
    object = sme_estimation,
    vcov. = vcov(sme_estimation),
    g = "1 - pnorm(theta_01)",
    parameterNames = c("theta_01", "theta_02", "err", "mu")
  )
  SME_JobExit_rate <- SME_JobExit[1][[1]]
  SME_JobEntry <- car::deltaMethod(
    object = sme_estimation,
    vcov. = vcov(sme_estimation),
    g = "pnorm(theta_02)",
    parameterNames = c("theta_01", "theta_02", "err", "mu")
  )
  SME_JobEntry_rate <- SME_JobEntry[1][[1]]
  SME_Empl <- car::deltaMethod(
    object = sme_estimation,
    vcov. = vcov(sme_estimation),
    g = "pnorm(mu)",
    parameterNames = c("theta_01", "theta_02", "err", "mu")
  )
  SME_Empl_rate <- SME_Empl[1][[1]]
  SME_Misclass <- car::deltaMethod(
    object = sme_estimation,
    vcov. = vcov(sme_estimation),
    g = "1 - pnorm(err)",
    parameterNames = c("theta_01", "theta_02", "err", "mu")
  )
  SME_Misclass_rate <- SME_Misclass[1][[1]]

  # Implied SEs ----
  SME_JobExit_se <- SME_JobExit[2][[1]]
  SME_JobEntry_se <- SME_JobEntry[2][[1]]
  SME_Empl_se <- SME_Empl[2][[1]]
  SME_Misclass_se <- SME_Misclass[2][[1]]


  # Implied CIs ----
  coef5 <- qnorm(1 - 0.05/2)
  coef1 <- qnorm(1 - 0.01/2)
  coef01 <- qnorm(1 - 0.001/2)
  SME_CI_5 <- data.table(
    "Prob" = c("poverty_exit", "poverty_entry", "level", "misclass"),
    "LowerBound" = c(SME_JobExit_rate  - coef5*SME_JobExit_se,
                     SME_JobEntry_rate - coef5*SME_JobEntry_se,
                     SME_Empl_rate     - coef5*SME_Empl_se,
                     SME_Misclass_rate - coef5*SME_Misclass_se),
    "UpperBound" = c(SME_JobExit_rate  + coef5*SME_JobExit_se,
                     SME_JobEntry_rate + coef5*SME_JobEntry_se,
                     SME_Empl_rate     + coef5*SME_Empl_se,
                     SME_Misclass_rate + coef5*SME_Misclass_se)
  )
  SME_CI_1 <- data.table(
    "Prob" = c("poverty_exit", "poverty_entry", "level", "misclass"),
    "LowerBound" = c(SME_JobExit_rate  - coef1*SME_JobExit_se,
                     SME_JobEntry_rate - coef1*SME_JobEntry_se,
                     SME_Empl_rate     - coef1*SME_Empl_se,
                     SME_Misclass_rate - coef1*SME_Misclass_se),
    "UpperBound" = c(SME_JobExit_rate  + coef1*SME_JobExit_se,
                     SME_JobEntry_rate + coef1*SME_JobEntry_se,
                     SME_Empl_rate     + coef1*SME_Empl_se,
                     SME_Misclass_rate + coef1*SME_Misclass_se)
  )
  SME_CI_01 <- data.table(
    "Prob" = c("poverty_exit", "poverty_entry", "level", "misclass"),
    "LowerBound" = c(SME_JobExit_rate  - coef01*SME_JobExit_se,
                     SME_JobEntry_rate - coef01*SME_JobEntry_se,
                     SME_Empl_rate     - coef01*SME_Empl_se,
                     SME_Misclass_rate - coef01*SME_Misclass_se),
    "UpperBound" = c(SME_JobExit_rate  + coef01*SME_JobExit_se,
                     SME_JobEntry_rate + coef01*SME_JobEntry_se,
                     SME_Empl_rate     + coef01*SME_Empl_se,
                     SME_Misclass_rate + coef01*SME_Misclass_se)
  )


  #_____________________________________________________________________________
  # Results object -------------------------------------------------------------
  sme_results <- list(
    "model_type"         = model_name,
    "estimated_model"    = sme_estimation,
    "mbo_model"          = global_model,
    "loglik"             = logLik(sme_estimation),
    "model_summary"      = sme_estimation %>% summary(),
    "poverty_entry_rate" = SME_JobEntry_rate,
    "poverty_exit_rate"  = SME_JobExit_rate,
    "poverty_rate"       = SME_Empl_rate,
    "misclass_rate"      = SME_Misclass_rate,
    "poverty_entry_se"   = SME_JobEntry_se,
    "poverty_exit_se"    = SME_JobExit_se,
    "poverty_se"         = SME_Empl_se,
    "misclass_se"        = SME_Misclass_se,
    "implied_CI_5perc"   = SME_CI_5,
    "implied_CI_1perc"   = SME_CI_1,
    "implied_CI_01perc"  = SME_CI_01
  )


  #_____________________________________________________________________________
  # Return -----------------------------------------------------------
  return(sme_results)


}



