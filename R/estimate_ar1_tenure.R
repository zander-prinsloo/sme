
#' Estimate AR1 tenure classification error estimator
#'
#' @inheritParams estimate_ar1_sym
#' @inheritParams get_ar1_tenure_full_likelihood
#' @param globa_n number of samples design - global estimation.
#' See `mlrMBO:::generateDesign`
#'
#' @return list
#' @export
estimate_ar1_tenure <- function(
    st3_observed         = NULL,
    st2_observed         = NULL,
    st1_observed         = NULL,
    gt3                  = NULL,
    gt2                  = NULL,
    gt1                  = NULL,
    ht3                  = NULL,
    ht2                  = NULL,
    ht1                  = NULL,
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
  # if (is.null(init_params)) {
  #   init_params <- c(
  #     "theta_01" = 1.5,
  #     "theta_02" = -1.5,
  #     "err"      = 1.96,
  #     "mu"       = 0.5
  #   )
  # } else if (!names(init_params) %chin% c("theta_01", "theta_02", "err", "mu")) {
  #   names(init_params) <- c("theta_01", "theta_02", "err", "mu")
  # }

  #_____________________________________________________________________________
  # Obj Function----------------------------------------------------------------
  fn_ll <- function(par_vec){

    # Specify Parameters
    theta_1  <- par_vec[1]
    theta_2  <- par_vec[2]
    sigma    <- par_vec[3]
    lambda_g <- par_vec[4]
    lambda_h <- par_vec[5]
    err      <- par_vec[6]
    mu       <- par_vec[7]

    log_lik <- get_ar1_tenure_full_likelihood(
      st1_observed   = st1_observed,
      st2_observed   = st2_observed,
      st3_observed   = st3_observed,
      g_1            = gt1,
      g_2            = gt2,
      g_3            = gt3,
      h_1            = ht1,
      h_2            = ht2,
      h_3            = ht3,
      theta_1        = theta_1,
      theta_2        = theta_2,
      sigma          = sigma,
      lambda_g       = lambda_g,
      lambda_h       = lambda_h,
      err            = err,
      mu             = mu,
      skip_checks    = TRUE
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

      theta_1_lower  <-  1
      theta_1_upper  <-  2.5

      theta_2_lower  <- -2.5
      theta_2_upper  <- -1

      sigma_lower    <-  0.00001
      sigma_upper    <-  1

      lambda_g_lower <-  5
      lambda_g_upper <-  8

      lambda_h_lower <-  1
      lambda_h_upper <-  4

      err_lower      <-  0.5
      err_upper      <-  3

      mu_lower       <- -0.5
      mu_upper       <-  0.5

    } #else if (
      #!is.list(global_specification) |
      # !names(global_specification) %chin% c(
      #   "theta_1_lower",
      #   "theta_1_upper",
      #   "theta_2_lower",
      #   "theta_2_upper",
      #   "sigma_upper",
      #   "sigma_lower",
      #   "lambda_g_upper",
      #   "lambda_g_lower",
      #   "lambda_h_upper",
      #   "lambda_h_lower",
      #   "err_lower",
      #   "err_upper",
      #   "mu_lower",
      #   "mu_upper"
      # )
    # ) {
    #   cli::cli_abort(
    #     "If global specification is non NULL, then should be a list giving upper
    #     and lower bounds for each parameter"
    #   )
    #}
    else {
      theta_1_lower  <- global_specification$theta_1_lower
      theta_1_upper  <- global_specification$theta_1_upper
      theta_2_lower  <- global_specification$theta_2_lower
      theta_2_upper  <- global_specification$theta_2_upper
      sigma_upper    <- global_specification$sigma_upper
      sigma_lower    <- global_specification$sigma_lower
      lambda_g_upper <- global_specification$lambda_g_upper
      lambda_g_lower <- global_specification$lambda_g_lower
      lambda_h_upper <- global_specification$lambda_h_upper
      lambda_h_lower <- global_specification$lambda_h_lower
      err_lower      <- global_specification$err_lower
      err_upper      <- global_specification$err_upper
      mu_lower       <- global_specification$mu_lower
      mu_upper       <- global_specification$mu_upper
    }

    # Create Parameter Space
    par_space <- makeParamSet(
      makeNumericParam(
        "theta_1",
        lower = theta_1_lower,
        upper = theta_1_upper
      ),
      makeNumericParam(
        "theta_2",
        lower = theta_2_lower,
        upper = theta_2_upper
      ),
      makeNumericParam(
        "sigma",
        lower = sigma_lower,
        upper = sigma_upper
      ),
      makeNumericParam(
        "lambda_g",
        lower = lambda_g_lower,
        upper = lambda_g_upper
      ),
      makeNumericParam(
        "lambda_h",
        lower = lambda_h_lower,
        upper = lambda_h_upper
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
      noisy    = verbose,
      par.set  = par_space
    )
    sur_learner <- mlr::makeLearner(
      cl           = "regr.km",
      predict.type = "se",
      covtype      = "matern3_2",
      control      = list(trace = FALSE)
    )
    control_object <- mlrMBO::makeMBOControl() |>
      mlrMBO::setMBOControlTermination(
                                       iters = global_iterations) |>
      mlrMBO::setMBOControlInfill(
                                  crit = mlrMBO::makeMBOInfillCritEI())

    # Define initial search ----
    set.seed(123)
    df_initial_search <- generateDesign(
      n       = global_n,
      par.set = par_space,
      fun     = lhs::randomLHS
    )

    df_initial_search <- df_initial_search |>
      rbind(
        init_params,
        # c(
        #   "theta_1"   = 1.6,
        #   "theta_2"   = -1.6,
        #   "sigma"     = 0.15,
        #   "lambda_g"  = 0.3,
        #   "lambda_h"  = 0.3,
        #   "err"       = 1.6,
        #   "mu"        = 0.1
        # ),
        # c(
        #   "theta_1"   = 1.7,
        #   "theta_2"   = -1.7,
        #   "sigma"     = 0.05,
        #   "lambda_g"  = 0.4,
        #   "lambda_h"  = 0.4,
        #   "err"       = 1.5,
        #   "mu"        = -0.1
        # ),
        # c(
        #   "theta_1"   = 1.32,
        #   "theta_2"   = -1.87,
        #   "sigma"     = 2.47e-05,
        #   "lambda_g"  = 4.85,
        #   "lambda_h"  = 2.52,
        #   "err"       = 0.00324,
        #   "mu"        = -0.292
        #),
        # c(
        #   "theta_1"   = 1.32,
        #   "theta_2"   = -1.87,
        #   "sigma"     = 2.47e-05,
        #   "lambda_g"  = 4.85,
        #   "lambda_h"  = 2.52,
        #   "err"       = 0.00324,
        #   "mu"        = -0.292
        # ),
        # c(
        #   "theta_1"   = 1.2950126839  ,
        #   "theta_2"   = -2.2073950339  ,
        #   "sigma"     = 0.0001000000 ,
        #   "lambda_g"  = 5.9594590423  ,
        #   "lambda_h"  = 1.5448323110,
        #   "err"       = 0.0001171644,
        #   "mu"        = 0.2465914012
        # ),
        c(
          "theta_1"   = 1.8,
          "theta_2"   = -1.8,
          "sigma"     = 0.2 ,
          "lambda_g"  = 6.5  ,
          "lambda_h"  = 2.7,
          "err"       = 1.98,
          "mu"        = 0.2465914012
        ),
        c(
          "theta_1"   = 1.96,
          "theta_2"   = -1.89,
          "sigma"     = 0.3 ,
          "lambda_g"  = 6.5  ,
          "lambda_h"  = 2.7,
          "err"       = 1.2,
          "mu"        = -0.2
        )

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
    names(init_params) <- names(global_model$x)

    if (verbose == T) {
      print(init_params)
    }
  } else {
    global_model <- NULL
  }

  #_____________________________________________________________________________
  # Local Estimation -----------------------------------------------------------
  # A_mat <- matrix(
  #   c(0, 0, 1, 0, 0, 0, 0,  # Constraint for sigma
  #     0, 0, 0, 1, 0, 0, 0,  # Constraint for lambda_g
  #     0, 0, 0, 0, 1, 0, 0,  # Constraint for lambda_h
  #     0, 0, 0, 0, 0, 1, 0), # Constraint for err
  #   ncol  = 7,
  #   byrow = TRUE
  # )
  # B_mat <- c(-0.000001, -0.000001, -0.000001, -0.000001)

  sme_estimation <- maxLik::maxLik(
    fn_ll,
    start = init_params,
    method = "NM"#,
    # constraints = list(
    #   ineqA = A_mat,
    #   ineqB = B_mat
    # )

  )
  if (verbose) {
    print(sme_estimation |> summary())
    print(sme_estimation$estimate |> names())
  }
  if (maxLik::returnCode(sme_estimation) == 1) {
    if (verbose) {
      cli::cli_alert_info("Re-estimating local by refreshing initial params")
    }
    sme_estimation <- maxLik::maxLik(
      fn_ll,
      start = sme_estimation$estimate,
      method = "NM"#,
      # constraints = list(
      #   ineqA = A_mat,
      #   ineqB = B_mat
      # )
    )
    if (verbose) {
      print(sme_estimation |> summary())
    }
  }
  if (maxLik::returnCode(sme_estimation) == 1) {
    if (verbose) {
      cli::cli_alert_info("Re-estimating local by refreshing initial params")
    }
    sme_estimation <- maxLik::maxLik(
      fn_ll,
      start = sme_estimation$estimate,
      method = "NM"#,
      # constraints = list(
      #   ineqA = A_mat,
      #   ineqB = B_mat
      # )
    )
    if (verbose) {
      print(sme_estimation |> summary())
    }
  }
  if (maxLik::returnCode(sme_estimation) == 1) {
    if (verbose) {
      cli::cli_alert_info("Re-estimating local by refreshing initial params")
    }
    sme_estimation <- maxLik::maxLik(
      fn_ll,
      start = sme_estimation$estimate,
      method = "NM"#,
      # constraints = list(
      #   ineqA = A_mat,
      #   ineqB = B_mat
      # )
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
  SME_JobExit <- car::deltaMethod(
    sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "1 - pnorm(theta_1)",
    parameterNames = c("theta_1",
                       "theta_2",
                       "sigma",
                       "lambda_g",
                       "lambda_h",
                       "err",
                       "mu")
  )
  SME_JobExit_rate <- SME_JobExit[1][[1]]
  SME_JobEntry <- car::deltaMethod(
    sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "pnorm(theta_2)",
    parameterNames = c("theta_1",
                       "theta_2",
                       "sigma",
                       "lambda_g",
                       "lambda_h",
                       "err",
                       "mu")
  )
  SME_JobEntry_rate <- SME_JobEntry[1][[1]]
  SME_Empl <- car::deltaMethod(
    sme_estimation,
    vcov.          = vcov(sme_estimation),
    g              = "pnorm(mu)",
    parameterNames = c("theta_1",
                       "theta_2",
                       "sigma",
                       "lambda_g",
                       "lambda_h",
                       "err",
                       "mu")
  )
  SME_Empl_rate <- SME_Empl[1][[1]]
  SME_Misclass <- car::deltaMethod(
    sme_estimation,
    vcov. = vcov(sme_estimation),
    g = "1 - pnorm(err)",
    parameterNames = c("theta_1",
                       "theta_2",
                       "sigma",
                       "lambda_g",
                       "lambda_h",
                       "err",
                       "mu")
  )
  SME_Misclass_rate <- SME_Misclass[1][[1]]

  # Implied SEs ----
  SME_JobExit_se  <- SME_JobExit[2][[1]]
  SME_JobEntry_se <- SME_JobEntry[2][[1]]
  SME_Empl_se     <- SME_Empl[2][[1]]
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
    "model_summary"      = sme_estimation |> summary(),
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

