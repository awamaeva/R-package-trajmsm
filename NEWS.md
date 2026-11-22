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
