#  df_sim_ten <- sim_tenure_ar1(
#    n        = 10000,
#    theta1   = 1.65,
#    theta2   = -1.65,
#    sigma    = 0.05,
#    lambda_g = 4,
#    lambda_h = 4,
#    mu       = 0,
#    err      = 1.65
# )
#
# results_tenure_sim <- estimate_ar1_tenure(
#    st3_observed         = df_sim_ten$s3obs,
#    st2_observed         = df_sim_ten$s2obs,
#    st1_observed         = df_sim_ten$s1obs,
#    gt3                  = df_sim_ten$g3obs, # must be annual
#    gt2                  = df_sim_ten$g2obs, # must be annual
#    gt1                  = df_sim_ten$g1obs, # must be annual
#    ht3                  = df_sim_ten$h3obs, # must be annual
#    ht2                  = df_sim_ten$h2obs, # must be annual
#    ht1                  = df_sim_ten$h1obs, # must be annual
#    weights              = rep(1,
#                               length(df_sim_ten$g3obs)),
#    model_name           = "sim_tenure_check", # not necessary
#    init_params          = c(                  # if not global, then
#                                               # starting params for maxLik
#                                               # if global, then added to
#                                               # initial grid
#      "theta_1"   = 1.6,
#      "theta_2"   = -1.6,
#      "sigma"     = -2.3,
#      "lambda_g"  = 3.1,
#      "lambda_h"  = 2.9,
#      "err"       = 1.6,
#      "mu"        = -0.1
#    ),
#    global_est           = TRUE, # if TRUE, then first does mlrMBO global
#    verbose              = TRUE, # gives info as estimation goes
#    global_n             = 50,   # size of initial grid for mlrMBO, only if global
#    global_iterations    = 50    # number of global est. iterations, only if global
# )


# # Specify your own bounds for global specification
# results_tenure_sim <- estimate_ar1_tenure(
#    st3_observed         = df_sim_ten$s3obs,
#    st2_observed         = df_sim_ten$s2obs,
#    st1_observed         = df_sim_ten$s1obs,
#    gt3                  = df_sim_ten$g3obs,
#    gt2                  = df_sim_ten$g2obs,
#    gt1                  = df_sim_ten$g1obs,
#    ht3                  = df_sim_ten$h3obs,
#    ht2                  = df_sim_ten$h2obs,
#    ht1                  = df_sim_ten$h1obs,
#    weights              = rep(1, length(df_sim_ten$g3obs)),
#    model_name           = "sim_tenure_check",
#    init_params          = c(
#                    "theta_1"   = 1.6,
#                    "theta_2"   = -1.6,
#                    "sigma"     = -2.3,
#                    "lambda_g"  = 3.1,
#                    "lambda_h"  = 2.9,
#                    "err"       = 1.6,
#                    "mu"        = -0.1
#    ),
#    global_specification = list(         # customize global bounds for mlrMBO
#                theta_1_lower   =  1,
#                theta_1_upper   =  2.5,
#                theta_2_lower   = -2.5,
#                theta_2_upper   = -1,
#                sigma_lower     =  -5,
#                sigma_upper     =  1,
#                lambda_g_lower  =  1,
#                lambda_g_upper  =  8,
#                lambda_h_lower  =  1,
#                lambda_h_upper  =  8,
#                err_lower       =  0.5,
#                err_upper       =  3,
#                mu_lower        = -0.5,
#                mu_upper        =  0.5
#    ),
#    global_est           = TRUE,
#    verbose              = TRUE,
#    global_n             = 50,
#    global_iterations    = 50
# )


