# trajmsm 0.1.7

- Fix weight handling in `trajmsm_ipw()` for binary and Gaussian outcomes.
  These models use the supplied or computed weights and cluster by the
  specified ID column, regardless of how the caller names its data. Results
  may change if a previous fit silently ignored weights or used a different
  clustering column.
- Retain the 0.1.6 Cox inference correction described below. Function
  arguments, return values, and `treshold` behavior are unchanged.

# trajmsm 0.1.6

- Fix standard errors, Wald p-values, and log-hazard confidence intervals in
  `trajmsm_ipw(family = "survival")`. Standard errors now come from the
  fitted Cox model covariance rather than the exponentiated coefficient.
  The fitted model, its coefficients, and its weighting and clustering
  behavior are unchanged from CRAN 0.1.5.
- Add a regression test comparing ordinary and robust inference to direct
  Cox fits, including subject clustering for repeated start-stop rows.

This update does not change the history-restricted estimators or the
`treshold` weight handling.
