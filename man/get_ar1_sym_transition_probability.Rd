% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ar1_sym.R
\name{get_ar1_sym_transition_probability}
\alias{get_ar1_sym_transition_probability}
\title{Observed ransition probability, which is conditional
distribution of observed period 3 (t3) status given
observed statuses in t2 and t1 - P(st3 | st2, st1)}
\usage{
get_ar1_sym_transition_probability(
  st3_observed,
  st2_observed,
  st1_observed,
  theta_01,
  theta_02,
  err,
  mu
)
}
\arguments{
\item{st2_observed}{numeric: 0 or 1 - observed status in period t2}

\item{st1_observed}{numeric: 0 or 1 - observed status in period t1}

\item{theta_01}{numeric: parameter of true transition process specifying
the probability (after applying probit) exiting status 1}

\item{theta_02}{numeric: parameter of true transition process specifying
the probability (after applying probit) of entering status 1}

\item{err}{numeric: parameter of noise process. The probability of
misclassification is equal to 1-Phi(err), where Phi() is the probit
transformation}

\item{mu}{numeric: parameter of true process. Specifies
unconditional distribution in t1, where P(st1_true = 1) = mu}

\item{st3}{numeric: 0 or 1 - observed status in period t2}
}
\value{
numeric
}
\description{
Observed ransition probability, which is conditional
distribution of observed period 3 (t3) status given
observed statuses in t2 and t1 - P(st3 | st2, st1)
}
\keyword{internal}
