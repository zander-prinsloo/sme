# df_sim <- sim_tenure_ar1(
#   n        = 10000,
#   theta1   = 1.65,
#   theta2   = -1.65,
#   sigma    = 0.1,
#   lambda_g = 3,
#   lambda_h = 3,
#   mu       = 0,
#   err      = 1.65
# )
# table(df_sim$s1obs, df_sim$s2obs)
# table(df_sim$s2obs, df_sim$s3obs)
# table(df_sim$s1obs, df_sim$s3obs)
#
# table(df_sim$s1true, df_sim$s2true)
# table(df_sim$s2true, df_sim$s3true)
# table(df_sim$s1true, df_sim$s3true)
#
# get_ar1_sym_joint_probability(
#   st2_observed = 0,
#   st1_observed = 0,
#   st2_true     = c(0, 1, 0, 1),
#   st1_true     = c(0, 0, 1, 1),
#   theta_1      = 1.65,
#   theta_2      = -1.65,
#   err          = 1.65,
#   mu           = 0
# ) |>
#   sum()
#
#
# get_ar1_sym_individual_likelihood(
#   st3_observed = 1,
#   st2_observed = 1,
#   st1_observed = 1,
#   theta_1      = 1.65,
#   theta_2      = -1.65,
#   err          = 1.65,
#   mu           = 0
# )*10000
# df_sim |>
#   fsubset(
#     s1obs == 1 & s2obs == 1 & s3obs == 1
#   ) |>
#   nrow()
#
# get_ar1_sym_individual_likelihood(
#   st3_observed = 0,
#   st2_observed = 1,
#   st1_observed = 1,
#   theta_1      = 1.65,
#   theta_2      = -1.65,
#   err          = 1.65,
#   mu           = 0
# )*10000
# df_sim |>
#   fsubset(
#     s1obs == 1 & s2obs == 1 & s3obs == 0
#   ) |>
#   nrow()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# results_ar1_sym_simulation <- maxLik::maxLik(
#   logLik = function(params) {
#     theta_1 <- params[1]
#     theta_2 <- params[2]
#     err     <- params[3]
#     mu      <- params[4]
#
#     lik <- get_ar1_sym_full_likelihood(
#       st3_observed = df_sim$s3obs,
#       st2_observed = df_sim$s2obs,
#       st1_observed = df_sim$s1obs,
#       par_vec = c(
#         1.6,
#         -1,6,
#         1.6,
#         0
#       )
#       # theta_1      = theta_1,
#       # theta_2      = theta_2,
#       # err          = err,
#       # mu           = mu
#     ) |>
#       log() |>
#       sum()
#
#     lik
#   },
#   method = "NM",
#   start = c(
#     "theta_1"      = 1.6,
#     "theta_2"      = -1.6,
#     "err"          = 1.6,
#     "mu"           = 0
#   ),
#   reltol = 1e-12
# )
#
# get_ar1_sym_individual_likelihood(
#   st3_observed = df_sim$s3obs,
#   st2_observed = df_sim$s2obs,
#   st1_observed = df_sim$s1obs,
#   theta_1      = 1.65,
#   theta_2      = -1.65,
#   err          = 1.65,
#   mu           = 0
# ) |>
#   log() |>
#   sum()
#
# get_ar1_sym_individual_likelihood(
#   st3_observed = 1,
#   st2_observed = 1,
#   st1_observed = 1,
#   theta_1      = 1.65,
#   theta_2      = -1.65,
#   err          = 1.65,
#   mu           = 0
# )
#
#

#
#
#
# results_ar1_sym_simulation_FUNC <- estimate_ar1_sym(
#   st3_observed     = df_sim$s3obs,
#   st2_observed     = df_sim$s2obs,
#   st1_observed     = df_sim$s1obs,
#   weights          = rep(1, length(df_sim$s1obs)),
#   init_params      = c(
#     "theta_1"      = 1.6,
#     "theta_2"      = -1.6,
#     "err"          = 1.6,
#     "mu"           = 0
#   ),
#   global_est       = TRUE,
#   verbose          = TRUE
# )
#
#










