# df_qlfs <- readRDS("C:/Users/rulof/Gitlab/sme/data/df_qlfs-8Feb2024.RDS")

ex_gaussian_density <- function(x, sigma, lambda) {
  # print(lambda)
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


library(tidyverse)
source("R/get_ar1_tenure_two_sigma_em_algorithm_rulof_trial.R", echo=TRUE)
source("R/estimate_ar1_tenure_two_sigma_em_algorithm_rulof.R", echo=TRUE)
# library(ParamHelpers)
# library(collapse)
# library(data.table)


fn_QLFS_Multiple_Wave <- function(){

  # Import
  df_1 <- haven::read_dta(
    "data/qlfs_raw_4waves.dta"
  )

  # Keep only age, educ, status for 3 waves
  df_1 <- df_1 %>%
    select(contains(c("age", "educ", "status", "tenure", "timegap", "formal"))) %>%  # keep age, educ, empl status
    select(-contains(c("4"))) # remove columns on wage, and 4th wave
  #select(-contains(c("wage", "4"))) # remove columns on wage, and 4th wave

  # Make Employment Binary
  df_1 <- df_1 %>%
    mutate(
      Status_Q1 = ifelse(
        status1 == 1, 1, 0
      ),
      Status_Q2 = ifelse(
        status2 == 1, 1, 0
      ),
      Status_Q3 = ifelse(
        status3 == 1, 1, 0
      ),
      tenure1 = if_else(Status_Q1 == 0, 0, tenure1),
      tenure2 = if_else(Status_Q2 == 0, 0, tenure2),
      tenure3 = if_else(Status_Q3 == 0, 0, tenure3),
      timegap1 = if_else(Status_Q1 == 1, 0, timegap1),
      timegap2 = if_else(Status_Q2 == 1, 0, timegap2),
      timegap3 = if_else(Status_Q3 == 1, 0, timegap3)
      # tenure1 = ifelse(is.na(tenure1), 0, tenure1),
      # tenure2 = ifelse(is.na(tenure2), 0, tenure2),
      # tenure3 = ifelse(is.na(tenure3), 0, tenure3),
      # timegap1 = ifelse(is.na(timegap1), 0,
      #                   ifelse(timegap1>7, 0, timegap1)),
      # timegap2 = ifelse(is.na(timegap2), 0,
      #                   ifelse(timegap2>7, 0, timegap2)),
      # timegap3 = ifelse(is.na(timegap3), 0,
      #                   ifelse(timegap2>7, 0, timegap3))

    ) %>%
    select(-c(status1, status2, status3)) # remove other statuses

  # Age between 18 and 55
  df_1 <- df_1 %>%
    filter(age1 > 17 & age1 < 56)

  # Return
  return(df_1)

}

# Use above function
df_qlfs <- fn_QLFS_Multiple_Wave()

# Create timegap variables
df_qlfs <- df_qlfs %>%
  mutate(
    TG_Q1 = case_when(
      timegap1 == 0 ~ 0,
      timegap1 == 1 ~ 1.5,
      timegap1 == 2 ~ 4.5,
      timegap1 == 3 ~ 7.5,
      timegap1 == 4 ~ 10.5,
      timegap1 == 5 ~ 24,
      timegap1 == 6 ~ 48,
      timegap1 == 7 ~ 90,
      timegap1 == 8 ~ NA_real_,
      timegap1 == 99 ~ NA_real_
    )
  )%>%
  mutate(
    TG_Q2 = case_when(
      timegap2 == 0 ~ 0,
      timegap2 == 1 ~ 1.5,
      timegap2 == 2 ~ 4.5,
      timegap2 == 3 ~ 7.5,
      timegap2 == 4 ~ 10.5,
      timegap2 == 5 ~ 24,
      timegap2 == 6 ~ 48,
      timegap2 == 7 ~ 90,
      timegap2 == 8 ~ NA_real_,
      timegap2 == 99 ~ NA_real_
    )
  )%>%
  mutate(
    TG_Q3 = case_when(
      timegap3 == 0 ~ 0,
      timegap3 == 1 ~ 1.5,
      timegap3 == 2 ~ 4.5,
      timegap3 == 3 ~ 7.5,
      timegap3 == 4 ~ 10.5,
      timegap3 == 5 ~ 24,
      timegap3 == 6 ~ 48,
      timegap3 == 7 ~ 90,
      timegap3 == 8 ~ NA_real_,
      timegap3 == 99 ~ NA_real_
    )
  )

# Age Inconsistency
df_qlfs <- df_qlfs %>%
  mutate(
    Age_Inc_Q1 = ifelse(
      (age1>age2 | abs(age1-age2)>1) & !(age2>age3 | abs(age2-age3)>1), 1, 0
    )) %>%
  mutate(
    Age_Inc_Q2 = ifelse(
      (age1>age2 | abs(age1-age2)>1) & (age2>age3 | abs(age2-age3)>1), 1, 0
    )) %>%
  mutate(
    Age_Inc_Q3 = ifelse(
      !(age1>age2 | abs(age1-age2)>1) & (age2>age3 | abs(age2-age3)>1), 1, 0
    ))

df_qlfs <- df_qlfs %>%
  filter(!is.na(Status_Q1)) %>%
  filter(!is.na(Status_Q2)) %>%
  filter(!is.na(Status_Q3)) %>%
  filter(!is.na(tenure1)) %>%
  filter(!is.na(tenure2)) %>%
  filter(!is.na(tenure3)) %>%
  filter(!is.na(TG_Q1)) %>%
  filter(!is.na(TG_Q2)) %>%
  filter(!is.na(TG_Q3))

# df_qlfs <- df_qlfs[1:10,]

df_posterior_probs0 <- get_ar1_tenure_full_likelihood_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  g_1                  = df_qlfs$tenure1/12, # must be annual
  g_2                  = df_qlfs$tenure2/12, # must be annual
  g_3                  = df_qlfs$tenure3/12, # must be annual
  h_1                  = df_qlfs$TG_Q1/12, # must be annual
  h_2                  = df_qlfs$TG_Q2/12, # must be annual
  h_3                  = df_qlfs$TG_Q3/12, # must be annual
  posterior_probs = NULL,
  theta_1   = 1.950167,
  theta_2   = -1.479017,
  sigma_g   = log(0.05),
  sigma_h   = log(2),
  lambda_g  = log(6.625821),
  lambda_h  = log(2.727932),
  err       = 1.906552,
  return_posterior = TRUE
)

rowSums(df_posterior_probs0)


em_estimates1 <- estimate_ar1_tenure_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  gt1                  = df_qlfs$tenure1/12, # must be annual
  gt2                  = df_qlfs$tenure2/12, # must be annual
  gt3                  = df_qlfs$tenure3/12, # must be annual
  ht1                  = df_qlfs$TG_Q1/12, # must be annual
  ht2                  = df_qlfs$TG_Q2/12, # must be annual
  ht3                  = df_qlfs$TG_Q3/12, # must be annual
  weights              = rep(1,
                             length(df_qlfs$Status_Q1)),
  model_name           = "sim_tenure_check", # not necessary
  global_est           = FALSE, # if TRUE, then first does mlrMBO global
  init_params = c(
    theta_1   = 1.950167,
    theta_2   = -1.479017,
    sigma_g   = log(0.05),
    sigma_h   = log(2),
    lambda_g  = log(6.625821),
    lambda_h  = log(2.727932),
    err       = 1.906552
  ),
  posterior_probs = df_posterior_probs0,
  verbose              = TRUE # gives info as estimation goes
  # global_n             = 50,   # size of initial grid for mlrMBO, only if global
  # global_iterations    = 50    # number of global est. iterations, only if global
)

df_posterior_probs1 <- get_ar1_tenure_full_likelihood_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  g_1                  = df_qlfs$tenure1/12, # must be annual
  g_2                  = df_qlfs$tenure2/12, # must be annual
  g_3                  = df_qlfs$tenure3/12, # must be annual
  h_1                  = df_qlfs$TG_Q1/12, # must be annual
  h_2                  = df_qlfs$TG_Q2/12, # must be annual
  h_3                  = df_qlfs$TG_Q3/12, # must be annual
  posterior_probs = NULL,
  theta_1   = 1.1746179,
  theta_2   = -0.9699788,
  sigma_g   = -3.0546923,
  sigma_h   = 0.1812088,
  lambda_g  = 1.8999608,
  lambda_h  = 1.2825365,
  err       = 0.9408933,
  return_posterior = TRUE
)

rowSums(df_posterior_probs1)


em_estimates2 <- estimate_ar1_tenure_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  gt1                  = df_qlfs$tenure1/12, # must be annual
  gt2                  = df_qlfs$tenure2/12, # must be annual
  gt3                  = df_qlfs$tenure3/12, # must be annual
  ht1                  = df_qlfs$TG_Q1/12, # must be annual
  ht2                  = df_qlfs$TG_Q2/12, # must be annual
  ht3                  = df_qlfs$TG_Q3/12, # must be annual
  weights              = rep(1,
                             length(df_qlfs$Status_Q1)),
  model_name           = "sim_tenure_check", # not necessary
  global_est           = FALSE, # if TRUE, then first does mlrMBO global
  init_params = c(
    theta_1   = 1.1746179,
    theta_2   = -0.9699788,
    sigma_g   = -3.0546923,
    sigma_h   = 0.1812088,
    lambda_g  = 1.8999608,
    lambda_h  = 1.2825365,
    err       = 0.9408933
  ),
  posterior_probs = df_posterior_probs1,
  verbose              = TRUE # gives info as estimation goes
  # global_n             = 50,   # size of initial grid for mlrMBO, only if global
  # global_iterations    = 50    # number of global est. iterations, only if global
)

df_posterior_probs2 <- get_ar1_tenure_full_likelihood_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  g_1                  = df_qlfs$tenure1/12, # must be annual
  g_2                  = df_qlfs$tenure2/12, # must be annual
  g_3                  = df_qlfs$tenure3/12, # must be annual
  h_1                  = df_qlfs$TG_Q1/12, # must be annual
  h_2                  = df_qlfs$TG_Q2/12, # must be annual
  h_3                  = df_qlfs$TG_Q3/12, # must be annual
  posterior_probs = NULL,
  theta_1   = 1.1495291,
  theta_2   = -0.9951772,
  sigma_g   = -3.2801397,
  sigma_h   = -0.3654636,
  lambda_g  = 1.9070932,
  lambda_h  = 1.2812205,
  err       = 0.7474981,
  return_posterior = TRUE
)

rowSums(df_posterior_probs2)

em_estimates3 <- estimate_ar1_tenure_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  gt1                  = df_qlfs$tenure1/12, # must be annual
  gt2                  = df_qlfs$tenure2/12, # must be annual
  gt3                  = df_qlfs$tenure3/12, # must be annual
  ht1                  = df_qlfs$TG_Q1/12, # must be annual
  ht2                  = df_qlfs$TG_Q2/12, # must be annual
  ht3                  = df_qlfs$TG_Q3/12, # must be annual
  weights              = rep(1,
                             length(df_qlfs$Status_Q1)),
  model_name           = "sim_tenure_check", # not necessary
  global_est           = FALSE, # if TRUE, then first does mlrMBO global
  init_params = c(
    theta_1   = 1.1495291,
    theta_2   = -0.9951772,
    sigma_g   = -3.2801397,
    sigma_h   = -0.3654636,
    lambda_g  = 1.9070932,
    lambda_h  = 1.2812205,
    err       = 0.7474981
  ),
  posterior_probs = df_posterior_probs2,
  verbose              = TRUE # gives info as estimation goes
  # global_n             = 50,   # size of initial grid for mlrMBO, only if global
  # global_iterations    = 50    # number of global est. iterations, only if global
)


df_posterior_probs3 <- get_ar1_tenure_full_likelihood_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  g_1                  = df_qlfs$tenure1/12, # must be annual
  g_2                  = df_qlfs$tenure2/12, # must be annual
  g_3                  = df_qlfs$tenure3/12, # must be annual
  h_1                  = df_qlfs$TG_Q1/12, # must be annual
  h_2                  = df_qlfs$TG_Q2/12, # must be annual
  h_3                  = df_qlfs$TG_Q3/12, # must be annual
  posterior_probs = NULL,
  theta_1   = 1.1785836,
  theta_2   = -1.0430183,
  sigma_g   = -3.5449949,
  sigma_h   = -0.8723082,
  lambda_g  = 1.8992190,
  lambda_h  = 1.3024427,
  err       = 0.6553965,
  return_posterior = TRUE
)

rowSums(df_posterior_probs3)

em_estimate3 <- estimate_ar1_tenure_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  gt1                  = df_qlfs$tenure1/12, # must be annual
  gt2                  = df_qlfs$tenure2/12, # must be annual
  gt3                  = df_qlfs$tenure3/12, # must be annual
  ht1                  = df_qlfs$TG_Q1/12, # must be annual
  ht2                  = df_qlfs$TG_Q2/12, # must be annual
  ht3                  = df_qlfs$TG_Q3/12, # must be annual
  weights              = rep(1,
                             length(df_qlfs$Status_Q1)),
  model_name           = "sim_tenure_check", # not necessary
  global_est           = FALSE, # if TRUE, then first does mlrMBO global
  init_params = c(
    theta_1   = 1.1785836,
    theta_2   = -1.0430183,
    sigma_g   = -3.5449949,
    sigma_h   = -0.8723082,
    lambda_g  = 1.8992190,
    lambda_h  = 1.3024427,
    err       = 0.6553965
  ),
  posterior_probs = df_posterior_probs3,
  verbose              = TRUE # gives info as estimation goes
  # global_n             = 50,   # size of initial grid for mlrMBO, only if global
  # global_iterations    = 50    # number of global est. iterations, only if global
)


df_posterior_probs4 <- get_ar1_tenure_full_likelihood_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  g_1                  = df_qlfs$tenure1/12, # must be annual
  g_2                  = df_qlfs$tenure2/12, # must be annual
  g_3                  = df_qlfs$tenure3/12, # must be annual
  h_1                  = df_qlfs$TG_Q1/12, # must be annual
  h_2                  = df_qlfs$TG_Q2/12, # must be annual
  h_3                  = df_qlfs$TG_Q3/12, # must be annual
  posterior_probs = NULL,
  theta_1   = 1.1868826,
  theta_2   = -1.0652721,
  sigma_g   = -3.8774172,
  sigma_h   = -1.2713265,
  lambda_g  = 1.8898253,
  lambda_h  = 1.3010004,
  err       = 0.5816923,
  return_posterior = TRUE
)

rowSums(df_posterior_probs4)

em_estimate4 <- estimate_ar1_tenure_two_sigma(
  st3_observed         = df_qlfs$Status_Q3,
  st2_observed         = df_qlfs$Status_Q2,
  st1_observed         = df_qlfs$Status_Q1,
  gt1                  = df_qlfs$tenure1/12, # must be annual
  gt2                  = df_qlfs$tenure2/12, # must be annual
  gt3                  = df_qlfs$tenure3/12, # must be annual
  ht1                  = df_qlfs$TG_Q1/12, # must be annual
  ht2                  = df_qlfs$TG_Q2/12, # must be annual
  ht3                  = df_qlfs$TG_Q3/12, # must be annual
  weights              = rep(1,
                             length(df_qlfs$Status_Q1)),
  model_name           = "sim_tenure_check", # not necessary
  global_est           = FALSE, # if TRUE, then first does mlrMBO global
  init_params = c(
    theta_1   = 1.1868826,
    theta_2   = -1.0652721,
    sigma_g   = -3.8774172,
    sigma_h   = -1.2713265,
    lambda_g  = 1.8898253,
    lambda_h  = 1.3010004,
    err       = 0.5816923
  ),
  posterior_probs = df_posterior_probs4,
  verbose              = TRUE # gives info as estimation goes
  # global_n             = 50,   # size of initial grid for mlrMBO, only if global
  # global_iterations    = 50    # number of global est. iterations, only if global
)


# li <- get_ar1_tenure_full_likelihood_two_sigma(
#   st3_observed         = df_qlfs$Status_Q3,
#   st2_observed         = df_qlfs$Status_Q2,
#   st1_observed         = df_qlfs$Status_Q1,
#   g_1                  = df_qlfs$tenure1/12, # must be annual
#   g_2                  = df_qlfs$tenure2/12, # must be annual
#   g_3                  = df_qlfs$tenure3/12, # must be annual
#   h_1                  = df_qlfs$TG_Q1/12, # must be annual
#   h_2                  = df_qlfs$TG_Q1/12, # must be annual
#   h_3                  = df_qlfs$TG_Q3/12, # must be annual
#   posterior_probs = df_posterior_probs0,
#   theta_1   = 1.950167,
#   theta_2   = -1.479017,
#   sigma_g   = log(0.05),
#   sigma_h   = log(2),
#   lambda_g  = log(6.625821),
#   lambda_h  = log(2.727932),
#   err       = 1.906552,
#   store_likelihoods = T
# )
# sum(li)

# results_tenure_sim <- estimate_ar1_tenure_two_sigma(
#   st3_observed         = df_qlfs$Status_Q3,
#   st2_observed         = df_qlfs$Status_Q2,
#   st1_observed         = df_qlfs$Status_Q1,
#   gt3                  = df_qlfs$tenure3/12 , # must be annual
#   gt2                  = df_qlfs$tenure2/12 , # must be annual
#   gt1                  = df_qlfs$tenure1/12, # must be annual
#   ht3                  = df_qlfs$TG_Q3/12, # must be annual
#   ht2                  = df_qlfs$TG_Q2/12, # must be annual
#   ht1                  = df_qlfs$TG_Q1/12, # must be annual
#   posterior_probs = df_posterior_probs0,
#   weights              = rep(1,
#                              length(df_qlfs$Status_Q1)),
#   global_est           = FALSE, # if TRUE, then first does mlrMBO global
#   verbose              = TRUE, # gives info as estimation goes
#   global_n             = 50,   # size of initial grid for mlrMBO, only if global
#   global_iterations    = 50    # number of global est. iterations, only if global
# )


# li <- get_ar1_tenure_full_likelihood_two_sigma(
#   st3_observed         = df_qlfs$Status_Q3,
#   st2_observed         = df_qlfs$Status_Q2,
#   st1_observed         = df_qlfs$Status_Q1,
#   g_1                  = df_qlfs$tenure1/12, # must be annual
#   g_2                  = df_qlfs$tenure2/12, # must be annual
#   g_3                  = df_qlfs$tenure3/12, # must be annual
#   h_1                  = df_qlfs$TG_Q1/12, # must be annual
#   h_2                  = df_qlfs$TG_Q1/12, # must be annual
#   h_3                  = df_qlfs$TG_Q3/12, # must be annual
#   posterior_probs = df_posterior_probs0,
#   theta_1   = 1.1711239,
#   theta_2   = -0.9626205,
#   sigma_g   = -3.0528036,
#   sigma_h   = 0.2180329,
#   lambda_g  = 1.9007228,
#   lambda_h  = 0.9083551,
#   err       = 0.9447185,
#   store_likelihoods = T
# )
# sum(li)
#
#
# df_posterior_probs1 <- get_ar1_tenure_full_likelihood_two_sigma(
#   st3_observed         = df_qlfs$Status_Q3,
#   st2_observed         = df_qlfs$Status_Q2,
#   st1_observed         = df_qlfs$Status_Q1,
#   g_1                  = df_qlfs$tenure1/12, # must be annual
#   g_2                  = df_qlfs$tenure2/12, # must be annual
#   g_3                  = df_qlfs$tenure3/12, # must be annual
#   h_1                  = df_qlfs$TG_Q1/12, # must be annual
#   h_2                  = df_qlfs$TG_Q1/12, # must be annual
#   h_3                  = df_qlfs$TG_Q3/12, # must be annual
#   posterior_probs = NULL,
#   theta_1   = 1.1711239,
#   theta_2   = -0.9626205,
#   sigma_g   = -3.0528036,
#   sigma_h   = 0.2180329,
#   lambda_g  = 1.9007228,
#   lambda_h  = 0.9083551,
#   err       = 0.9447185,
#   return_posterior = TRUE
# )

# rowSums(df_posterior_probs1)
#
# results_tenure_sim <- estimate_ar1_tenure_two_sigma(
#   st3_observed         = df_qlfs$Status_Q3,
#   st2_observed         = df_qlfs$Status_Q2,
#   st1_observed         = df_qlfs$Status_Q1,
#   gt3                  = df_qlfs$tenure3/12 , # must be annual
#   gt2                  = df_qlfs$tenure2/12 , # must be annual
#   gt1                  = df_qlfs$tenure1/12, # must be annual
#   ht3                  = df_qlfs$TG_Q3/12, # must be annual
#   ht2                  = df_qlfs$TG_Q2/12, # must be annual
#   ht1                  = df_qlfs$TG_Q1/12, # must be annual
#   posterior_probs = df_posterior_probs1,
#   weights              = rep(1,
#                              length(df_qlfs$Status_Q1)),
#   global_est           = TRUE, # if TRUE, then first does mlrMBO global
#   verbose              = TRUE, # gives info as estimation goes
#   global_n             = 50,   # size of initial grid for mlrMBO, only if global
#   global_iterations    = 50    # number of global est. iterations, only if global
# )
#
#
#
#
# start_time <- Sys.time()
# li <- get_ar1_tenure_full_likelihood_two_sigma(
#   st3_observed         = df_qlfs$Status_Q3,
#   st2_observed         = df_qlfs$Status_Q2,
#   st1_observed         = df_qlfs$Status_Q1,
#   g_1                  = df_qlfs$tenure1/12, # must be annual
#   g_2                  = df_qlfs$tenure2/12, # must be annual
#   g_3                  = df_qlfs$tenure3/12, # must be annual
#   h_1                  = df_qlfs$TG_Q1/12, # must be annual
#   h_2                  = df_qlfs$TG_Q1/12, # must be annual
#   h_3                  = df_qlfs$TG_Q3/12, # must be annual
#   posterior_probs = df_posterior_probs,
#   theta_1   = 1.950167,
#   theta_2   = -1.479017,
#   sigma_g   = log(0.05),
#   sigma_h   = log(2),
#   lambda_g  = log(6.625821),
#   lambda_h  = log(2.727932),
#   err       = 1.906552
# )
# end_time <- Sys.time()
# end_time - start_time
# sum(li)
#
# df_posterior_probs <- li$posterior_probs %>% select(-tot)
#
# li <- get_ar1_tenure_full_likelihood_two_sigma(
#   st3_observed         = df_qlfs$Status_Q3,
#   st2_observed         = df_qlfs$Status_Q2,
#   st1_observed         = df_qlfs$Status_Q1,
#   g_1                  = df_qlfs$tenure1/12, # must be annual
#   g_2                  = df_qlfs$tenure2/12, # must be annual
#   g_3                  = df_qlfs$tenure3/12, # must be annual
#   h_1                  = df_qlfs$TG_Q1/12, # must be annual
#   h_2                  = df_qlfs$TG_Q1/12, # must be annual
#   h_3                  = df_qlfs$TG_Q3/12, # must be annual
#   theta_1   = 1.960000,
#   theta_2   = -1.89000000,
#   sigma_g   = -2,
#   sigma_h   = -2,
#   lambda_g  = 6.5,
#   lambda_h  = 2.7,
#   err       = 1.2
# )
# sum(li)
#
#
#
# li <- get_ar1_tenure_full_likelihood_two_sigma(
#   st3_observed         = df_qlfs$Status_Q3,
#   st2_observed         = df_qlfs$Status_Q2,
#   st1_observed         = df_qlfs$Status_Q1,
#   g_1                  = df_qlfs$tenure1/12, # must be annual
#   g_2                  = df_qlfs$tenure2/12, # must be annual
#   g_3                  = df_qlfs$tenure3/12, # must be annual
#   h_1                  = df_qlfs$TG_Q1/12, # must be annual
#   h_2                  = df_qlfs$TG_Q1/12, # must be annual
#   h_3                  = df_qlfs$TG_Q3/12, # must be annual
#   theta_1   =  1.526654,
#   theta_2   = -1.612638,
#   sigma_g   = -37.662788,
#   sigma_h   = -0.178858,
#   lambda_g  = 1.867037,
#   lambda_h  = 0.93250,
#   err      = 0.377669
# )
# sum(li)






