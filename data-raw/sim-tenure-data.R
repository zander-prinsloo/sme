## code to prepare `sim-tenure-data` dataset goes here

sim_tenure_ar1 <- function(n        = 100000,
                           theta1   = 1.645,
                           theta2   = -1.645,
                           sigma    = 0.1, # 0.3-0.6
                           lambda_g = 3, # mean of 3
                           lambda_h = 3, # mean of 3
                           mu       = 0,
                           err      = 1.96,
                           seed     = 1234){
  set.seed(seed)
  # time 1-----------
  #__________________
  s1true <- rbinom(n    = n,
                   size = 1,
                   prob = pnorm(mu))
  # time 2-----------
  #__________________
  s2true <- sapply(
    s1true,
    FUN = function(x){
      if (x == 1) {
        r <- rbinom(n    = 1,
                    size = 1,
                    prob = pnorm(theta1))
      } else {
        r <- rbinom(n    = 1,
                    size = 1,
                    prob = pnorm(theta2))
      }
    }
  )
  # time 3-----------
  #__________________
  s3true <- sapply(
    s2true,
    FUN = function(x){
      if (x == 1) {
        r <- rbinom(n    = 1,
                    size = 1,
                    prob = pnorm(theta1))
      } else {
        r <- rbinom(n    = 1,
                    size = 1,
                    prob = pnorm(theta2))
      }
    }
  )


  # Misclass----------------------------------
  #___________________________________________
  err1 <- rbinom(n    = n,
                 size = 1,
                 prob = pnorm(err))
  err2 <- rbinom(n    = n,
                 size = 1,
                 prob = pnorm(err))
  err3 <- rbinom(n    = n,
                 size = 1,
                 prob = pnorm(err))
  # note: err == 1 => correctly classified

  # data frame-----------
  #______________________
  df <- data.frame(
    s1true = s1true,
    s2true = s2true,
    s3true = s3true,
    err1   = err1,
    err2   = err2,
    err3   = err3
  )

  # Obs Status----------------------------------
  #___________________________________________
  df <- df |>
    ftransform(
      s1obs = ifelse(err1 == 1, s1true, abs(1 - s1true)),
      s2obs = ifelse(err2 == 1, s2true, abs(1 - s2true)),
      s3obs = ifelse(err3 == 1, s3true, abs(1 - s3true))
    )

  s1obs <- df$s1obs
  s2obs <- df$s2obs
  s3obs <- df$s3obs
  # Tenure------------------------------------
  #___________________________________________
  g1true              <- rexp(n = n, rate = 1/lambda_g) # mean = 1/lambda
  g1true[s1true == 0] <- 0
  g2true              <- s2true*(g1true + 0.25) + (1 - s2true)*0
  g3true              <- s3true*(g2true + 0.25) + (1 - s3true)*0
  k1                  <- gamlss.dist::rexGAUS(n = n, mu = 0, sigma = sigma, nu = lambda_g) # mean = 1/lambda
  k2                  <- gamlss.dist::rexGAUS(n = n, mu = 0, sigma = sigma, nu = lambda_g)
  k3                  <- gamlss.dist::rexGAUS(n = n, mu = 0, sigma = sigma, nu = lambda_g)
  u_g1                <- rnorm(n = n, mean = 0, sd = sigma)
  u_g2                <- rnorm(n = n, mean = 0, sd = sigma)
  u_g3                <- rnorm(n = n, mean = 0, sd = sigma)
  g1obs               <- err1*(s1obs*(g1true + u_g1)) + (1 - err1)*(s1obs*k1)
  g2obs               <- err2*(s2obs*(g2true + u_g2)) + (1 - err2)*(s2obs*k2)
  g3obs               <- err3*(s3obs*(g3true + u_g3)) + (1 - err3)*(s3obs*k3)

  # Unempl Duration---------------------------
  #___________________________________________
  h1true              <- rexp(n = n, rate = 1/lambda_h) # mean = 1/lambda
  h1true[s1true == 1] <- 0
  h2true              <- (1 - s2true)*(h1true + 0.25) + s2true*0
  h3true              <- (1 - s3true)*(h2true + 0.25) + s3true*0
  l1                  <- gamlss.dist::rexGAUS(n = n, mu = 0, sigma = sigma, nu = lambda_h) # mean = 1/lambda
  l2                  <- gamlss.dist::rexGAUS(n = n, mu = 0, sigma = sigma, nu = lambda_h)
  l3                  <- gamlss.dist::rexGAUS(n = n, mu = 0, sigma = sigma, nu = lambda_h)
  u_h1                <- rnorm(n = n, mean = 0, sd = sigma)
  u_h2                <- rnorm(n = n, mean = 0, sd = sigma)
  u_h3                <- rnorm(n = n, mean = 0, sd = sigma)
  h1obs               <- err1*((1 - s1obs)*(h1true + u_h1)) + (1 - err1)*((1 - s1obs)*l1)
  h2obs               <- err2*((1 - s2obs)*(h2true + u_h2)) + (1 - err2)*((1 - s2obs)*l2)
  h3obs               <- err3*((1 - s3obs)*(h3true + u_h3)) + (1 - err3)*((1 - s3obs)*l3)

  # Final data--------------------------------
  #___________________________________________
  df <- df |>
    ftransform(
      g1true = g1true,
      g2true = g2true,
      g3true = g3true,
      g1obs  = g1obs,
      g2obs  = g2obs,
      g3obs  = g3obs,
      h1true = h1true,
      h2true = h2true,
      h3true = h3true,
      h1obs  = h1obs,
      h2obs  = h2obs,
      h3obs  = h3obs,
      k1     = k1,
      k2     = k2,
      k3     = k3,
      l1     = l1,
      l2     = l2,
      l3     = l3
    )
  df

}
# df_sim |> View()
# df_sim <- sim_tenure_ar1(
#   n = 10000
# )
#
# sim_true_ar1(n = 100)
#
#
# d <- sim_true_ar1(n = 100000)
#
# d$s1true |> mean()
# d$s2 |> mean()
# d$s3 |> mean()
#
# d |>
#   fgroup_by(s3obs, s2obs) |>
#   fnobs() |>
#   fungroup() |>
#   #fgroup_by(s1) |>
#   ftransform(job_rate = s1obs/(50000))

# ______________________________________________________________________________
# create data set in /data------------------------------------------------------
# usethis::use_data(df_sim,
#                   overwrite = TRUE)
#
#
#
#
#
#
#
# # TRUE
# #_________________________________________________________
# (s2true[s1true == 0] |> sum())/(n - (s1true |> sum()))
# 1 - (s2true[s1true == 1] |> sum())/(s1true |> sum())
#
# (s3true[s2true == 0] |> sum())/(n - (s2true |> sum()))
# 1 - (s3true[s2true == 1] |> sum())/(s2true |> sum())
#
# (s3true[s1true == 0] |> sum())/(n - (s1true |> sum()))
# 1 - (s3true[s1true == 1] |> sum())/(s1true |> sum())
#
#
# # OBSERVED
# #_________________________________________________________
# (s2obs[s1obs == 0] |> sum())/(n - (s1obs |> sum()))
# 1 - (s2obs[s1obs == 1] |> sum())/(s1obs |> sum())
#
# (s3obs[s2obs == 0] |> sum())/(n - (s2obs |> sum()))
# 1 - (s3obs[s2obs == 1] |> sum())/(s2obs |> sum())
#
# (s3obs[s1obs == 0] |> sum())/(n - (s1obs |> sum()))
# 1 - (s3obs[s1obs == 1] |> sum())/(s1obs |> sum())
#
#


# MISCLASSIFIED
#_________________________________________________________






