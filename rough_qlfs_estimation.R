dr <- "P:/02.personal/wb612474/research/4wave_empl"
pacman::p_load(collapse,
               tidyverse,
               maxLik,
               mlrMBO,
               data.table)

dta <- haven::read_dta(fs::path(dr,
                                "qlfs_raw_4waves.dta"))
dta <- dta |>
  ftransform(
    status1 = fcase(status1 == 0, 0,
                    status1 == 1, 1,
                    status1 == 2, 0,
                    status1 == 3, 0,
                    status1 == 4, 0,
                    default = NA),
    status2 = fcase(status1 == 0, 0,
                    status1 == 1, 1,
                    status1 == 2, 0,
                    status1 == 3, 0,
                    status1 == 4, 0,
                    default = NA),
    status3 = fcase(status1 == 0, 0,
                    status1 == 1, 1,
                    status1 == 2, 0,
                    status1 == 3, 0,
                    status1 == 4, 0,
                    default = NA),

    status4 = fcase(status1 == 0, 0,
                    status1 == 1, 1,
                    status1 == 2, 0,
                    status1 == 3, 0,
                    status1 == 4, 0,
                    default = NA)


  )


start_param_vals <- c(
  "theta11" = 0.00628550512,
  "theta12" = 0.95134742945,
  "theta21" = 0.04358475319,
  "theta22" = 0.00157493384,
  "err"     = 3.19535003589,
  "mu"      = 0.65967553971
)
# estimate
est_ar2_local_mu <- estimate_ar2(
  st4_observed = dta$status4,
  st3_observed = dta$status3,
  st2_observed = dta$status2,
  st1_observed = dta$status1,
  weights      = NULL,
  verbose      = TRUE,
  global_est   = FALSE,
  init_params  = start_param_vals,
  mu_end       = FALSE
)


#
# Estimate    Std. error
# theta11 0.00628550512 0.00002007432
# theta12 0.95134742945 0.00001650614
# theta21 0.04358475319 0.00001117982
# theta22 0.00157493384 0.00000071753
# err     3.19535003589 0.00002684253
# mu      0.65967553971 0.00001507750
# t value               Pr(> t)
# theta11    313.11 < 0.00000000000000022 ***
#   theta12  57635.99 < 0.00000000000000022 ***
#   theta21   3898.52 < 0.00000000000000022 ***
#   theta22   2194.94 < 0.00000000000000022 ***
#   err     119040.58 < 0.00000000000000022 ***
#   mu       43752.33 < 0.00000000000000022 ***
