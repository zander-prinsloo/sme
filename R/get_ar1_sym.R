
#' Log likelihood of observed data
#'
#' The log likelihood for vectors giving the observed status
#' for three periods
#'
#' @param st3_observed numeric vector: observed status in t3
#' @param st2_observed numeric vector: observed status in t2
#' @param st1_observed numeric vector: observed status in t1
#' @param par_vec parameter vector
#'
#' @return numeric
#' @export
get_ar1_sym_full_likelihood <- function(
    st3_observed,
    st2_observed,
    st1_observed,
    par_vec
){
  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------


    # Specify Parameters
    theta_01 <- par_vec[1]
    theta_02 <- par_vec[2]
    err      <- par_vec[3]
    mu       <- par_vec[4]

    log_lik <- get_ar1_sym_individual_likelihood(
      st3_observed = st3_observed,
      st2_observed = st2_observed,
      st1_observed = st1_observed,
      theta_01     = theta_01,
      theta_02     = theta_02,
      err          = err,
      mu           = mu
    ) |>
      log() |>
      sum()

    # Return
    return(log_lik)

  }







#' Individual likelihood function
#'
#' Likelihood of observed statuses in t1, t2, t3 for an
#' individual.
#'
#' @inheritParams get_ar1_sym_transition_probability
#'
#' @return numeric
#' @export
get_ar1_sym_individual_likelihood <- function(
    st3_observed,
    st2_observed,
    st1_observed,
    theta_01,
    theta_02,
    err,
    mu
  ){
  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------




  #_____________________________________________________________________________
  # Transition prob ------------------------------------------------------------
  ind_lik_a <- get_ar1_sym_transition_probability(
    st3_observed = st3_observed,
    st2_observed = st2_observed,
    st1_observed = st1_observed,
    theta_01     = theta_01,
    theta_02     = theta_02,
    err          = err,
    mu           = mu
  )

  #_____________________________________________________________________________
  # Joint prob -----------------------------------------------------------------
  ind_lik_b <- get_ar1_sym_joint_probability(
    st2_observed = st2_observed,
    st2_true     = c(
      rep(1, n),
      rep(0, n),
      rep(1, n),
      rep(0, n)
    ),
    st1_observed = st1_observed,
    st1_true     = c(
      rep(1, n),
      rep(1, n),
      rep(0, n),
      rep(0, n)
    ),
    theta_01     = theta_01,
    theta_02     = theta_02,
    err          = err,
    mu           = mu
  ) |>
    sum()

  # Final Product
  outcome <- ind_lik_a*ind_lik_b

  # Return
  return(outcome)

}








#' Observed ransition probability, which is conditional
#' distribution of observed period 3 (t3) status given
#' observed statuses in t2 and t1 - P(st3 | st2, st1)
#'
#' @param st3 numeric: 0 or 1 - observed status in period t2
#' @inheritParams get_ar1_sym_joint_probability
#'
#' @return numeric
#' @keywords internal
get_ar1_sym_transition_probability <- function(
    st3_observed,
    st2_observed,
    st1_observed,
    theta_01,
    theta_02,
    err,
    mu
){

  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------




  #_____________________________________________________________________________
  # Constant -------------------------------------------------------------------
  constant_1 <- (2*pnorm(err) - 1)*(
    pnorm(theta_02) + (
      pnorm(theta_01) -
        pnorm(theta_02)
    )*pnorm(theta_02)
  )

  constant_2 <- 1 - pnorm(err)

  constant_3 <- -(1 - pnorm(err))*(
    pnorm(theta_01) -
      pnorm(theta_02)
  )*(
    pnorm(theta_01) -
      pnorm(theta_02)
  )

  constant <- constant_1 + constant_2 + constant_3

  #_____________________________________________________________________________
  # Term 1 -----------------------------------------------------------------
  coef_1 <- (
    pnorm(theta_01) -
      pnorm(theta_02)
  )*(
    pnorm(theta_01) -
      pnorm(theta_02)
  )
  term_1 <- coef_1*st1

  #_____________________________________________________________________________
  # Term 2 -----------------------------------------------------------------
  term_2 <- 0

  #_____________________________________________________________________________
  # Term 3 - expectation epsilon t1 --------------------------------------------
  coef_3 <- (
    pnorm(theta_01) -
      pnorm(theta_02)
  )*(
    pnorm(theta_01) -
      pnorm(theta_02)
  )
  term_3 <- -coef_3*get_ar1_sym_expected_value_epsilon_t1(
    st2_observed = st2_observed,
    st1_observed = st1_observed,
    theta_01     = theta_01,
    theta_02     = theta_02,
    err          = err,
    mu           = mu
  )

  #_____________________________________________________________________________
  # Term 4 -----------------------------------------------------------------
  term_4 <- 0

  #_____________________________________________________________________________
  # Term 5 ---------------------------------------------------------------------
  coef_5 <- (
    pnorm(theta_01) -
      pnorm(theta_02)
  )
  coef_5 <- coef_5*(2*pnorm(err) - 1)
  term_5 <- coef_5*get_ar1_sym_expected_value_e_t2(
    st2_observed = st2_observed,
    st1_observed = st1_observed,
    theta_01     = theta_01,
    theta_02     = theta_02,
    err          = err,
    mu           = mu
  )

  # Prob transition to 1
  prob_1 <- constant + term_1 + term_2 + term_3 + term_4 + term_5

  # Transition in and out of employment
  output <- ifelse(
    st3 == 1,
    prob_1,
    1 - prob_1
  )

  # Return
  return(output)
  #list(
  #    constant, constant_1, constant_2, constant_3, term_1, term_2, term_3, term_4, term_5, output
  #)

}





#' Conditional expectation of $\epsilon_{t-2}$
#'
#' Note that $\epsilon_{t-2}$ is the notation in the paper. In the code
#' it is referred to as period 1.
#'
#' @inheritParams get_ar1_sym_joint_probability
#'
#' @return numeric: between 1 and -1
#' @keywords internal
get_ar1_sym_expected_value_epsilon_t1 <- function(
    st2_observed,
    st1_observed,
    theta_01,
    theta_02,
    err,
    mu
){
  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------

  n <- length(st2)
  # First component - first two terms
  a <- st1 - (1 - pnorm(err))

  # Final Term - numerator
  b_numer <- get_ar1_sym_joint_probability(
    st2_observed = st2_observed,
    st2_true     = c(
      rep(1, n),
      rep(0, n)
    ),
    st1_observed = st1_observed,
    st1_true     = c(
      rep(1, n),
      rep(1, n)
    ),
    theta_01     = theta_01,
    theta_02     = theta_02,
    err          = err,
    mu           = mu
  ) |>
    sum()


  # Final Term - denominator
  b_denom <- get_ar1_sym_joint_probability(
    st2_observed = st2_observed,
    st2_true     = c(
      rep(1, n),
      rep(0, n),
      rep(1, n),
      rep(0, n)
    ),
    st1_observed = st1_observed,
    st1_true     = c(
      rep(1, n),
      rep(1, n),
      rep(0, n),
      rep(0, n)
    ),
    theta_01     = theta_01,
    theta_02     = theta_02,
    err          = err,
    mu           = mu
  ) |>
    sum()

  b <- b_numer/b_denom
  b <- -(2*pnorm(err) - 1)*b

  # Final
  outcome <- a + b

  # Return
  return(outcome)

}






#' Conditional expectation of $e_{t-1}$
#'
#' Note that $e_{t-1}$ is the notation in the paper. In the code
#' it is referred to as period 2.
#'
#' @inheritParams get_ar1_sym_joint_probability
#'
#' @return numeric: between 1 and -1
#' @keywords internal
get_ar1_sym_expected_value_e_t2 <- function(
  st2_observed,
  st1_observed,
  theta_01,
  theta_02,
  err,
  mu
){
  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------

  n <- length(st2)
  #_____________________________________________________________________________
  # Joint Prob st1 & st2 -------------------------------------------------------
  joint_prob <- get_ar1_sym_joint_probability(
    st2_observed = st2_observed,
    st2_true     = c(
      rep(0, n),
      rep(0, n),
      rep(1, n),
      rep(1, n)
    ),
    st1_observed = st1_observed,
    st1_true     = c(
      rep(0, n),
      rep(1, n),
      rep(0, n),
      rep(1, n)
    ),
    theta_01     = theta_01,
    theta_02     = theta_02,
    err          = err,
    mu           = mu
  ) |>
    sum()

  #_____________________________________________________________________________
  # Calculations ---------------------------------------------------------------
  term_1a <- get_ar1_sym_joint_probability(
    st2_observed = st2_observed,
    st2_true     = c(
      rep(1, n),
      rep(1, n)
    ),
    st1_observed = st1_observed,
    st1_true     = c(
      rep(1, n),
      rep(0, n)
    ),
    theta_01     = theta_01,
    theta_02     = theta_02,
    err          = err,
    mu           = mu
  ) |>
    sum()
  term_1 <- term_1a/joint_prob

  # Term 2
  term_2 <- -pnorm(theta_02)

  # Term 3
  term_3 <- -(pnorm(theta_01) - pnorm(theta_02))
  term_3 <- term_3*(
    get_ar1_sym_joint_probability(
      st2_observed = st2_observed,
      st2_true     = c(
        rep(1, n),
        rep(0, n)
      ),
      st1_observed = st1_observed,
      st1_true     = c(
        rep(1, n),
        rep(1, n)
      ),
      theta_01     = theta_01,
      theta_02     = theta_02,
      err          = err,
      mu           = mu
    ) |>
      sum()
  )
  term_3 <- term_3/joint_prob

  # Combine
  output <- term_1 + term_2 + term_3

  # Return
  return(output)
}





#' Joint probability of t1 and t2 observed and true statuses
#'
#' This function gives the joint probability that some specified combination
#' of true and observed statuses (S) are experienced over to periods. The
#' This function requires specifying parameter values.
#'
#' @param st2_observed numeric: 0 or 1 - observed status in period t2
#' @param st2_true numeric: 0 or 1 - true status in period t2
#' @param st1_observed numeric: 0 or 1 - observed status in period t1
#' @param st1_true numeric: 0 or 1 - true status in period t1
#' @param theta_01 numeric: parameter of true transition process specifying
#' the probability (after applying probit) exiting status 1
#' @param theta_02 numeric: parameter of true transition process specifying
#' the probability (after applying probit) of entering status 1
#' @param err numeric: parameter of noise process. The probability of
#' misclassification is equal to 1-Phi(err), where Phi() is the probit
#' transformation
#' @param mu numeric: parameter of true process. Specifies
#' unconditional distribution in t1, where P(st1_true = 1) = mu
#'
#' @return numeric double: joint probability between 0 and 1
#' @keywords internal
get_ar1_sym_joint_probability <- function(
    st2_observed,
    st2_true,
    st1_observed,
    st1_true,
    theta_01,
    theta_02,
    err,
    mu
){

  #_____________________________________________________________________________
  # Arg Checks -----------------------------------------------------------------


  #_____________________________________________________________________________
  # Misclass Prob Q2 -----------------------------------------------------------
  a <- st2_observed*(
    pnorm(err)*st2_true +
      (1 - pnorm(err))*(1 - st2_true)
  ) +
    (1 - st2_observed)*(
      pnorm(err)*(1 - st2_true) +
        (1 - pnorm(err))*st2_true
    )

  #_____________________________________________________________________________
  # True Transition Prob Q1->Q2 ------------------------------------------------
  b_1 <- st2_true*(
    pnorm(theta_02) +
      (
        pnorm(theta_01) -
          pnorm(theta_02)
      )*st1_true
  )
  b_2 <- (1 - st2_true)*(
    1 - pnorm(theta_02) -
      (
        pnorm(theta_01) -
          pnorm(theta_02)
      )*st1_true
  )
  b <- b_1 + b_2

  #_____________________________________________________________________________
  # Misclass Prob Q1 -----------------------------------------------------------
  c <- st1_observed*(
    pnorm(err)*st1_true +
      (1 - pnorm(err))*(1 - st1_true)
  ) +
    (1 - st1_observed)*(
      pnorm(err)*(1 - st1_true) +
        (1 - pnorm(err))*st1_true
    )

  #_____________________________________________________________________________
  # Unconditional Q1 -----------------------------------------------------------
  d <- st1_true*pnorm(mu) +
    (1 - st1_true)*(1 - pnorm(mu))

  # Multiply them together
  final_product <- a*b*c*d

  # return
  return(final_product)

}
