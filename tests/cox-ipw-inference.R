library(trajmsm)
library(survival)

# Fixed, synthetic data with distinct event times and non-integer weights.
d <- data.frame(
  id = 1:10, x = c(0, 1, 1, 0, 1, 0, 0, 1, 0, 1),
  z = c(2, 0, -1, 1, -2, -1, 0, 3, 1, -3), start = 0, stop = 1:10,
  event = c(1, 1, 0, 1, 1, 0, 1, 1, 0, 0),
  w = c(.45, .60, 1.1, 1.2, 2.1, .8, 3.5, .9, 7.5, 1.4)
)
repeated <- rbind(transform(d, stop = stop / 2, event = 0),
                  transform(d, start = stop / 2))

check_against_coxph <- function(data, formula, weights, se_column, threshold = .99) {
  direct <- do.call(coxph, list(formula = formula, data = data, weights = weights))
  sm <- coef(summary(direct))
  expected_se <- sm[, se_column]
  estimate <- coef(direct)
  expected <- cbind(estimate, expected_se, sm[, "Pr(>|z|)"],
                    estimate - 1.96 * expected_se, estimate + 1.96 * expected_se)
  colnames(expected) <- c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")

  got <- trajmsm_ipw(formula2 = formula, family = "survival",
                     identifier = "id", obsdata = data, weights = weights,
                     treshold = threshold)
  stopifnot(isTRUE(all.equal(got, expected, tolerance = 1e-10)),
            any(abs(got[, "Std.Error"] - sm[, "exp(coef)"]) > 1e-5))
  got
}

f <- Surv(start, stop, event) ~ x + z
invisible(check_against_coxph(d, f, rep(1, nrow(d)), "se(coef)"))
weighted <- check_against_coxph(d, f, d$w, "robust se")
stopifnot(identical(weighted, check_against_coxph(d, f, d$w, "robust se", 1)))

# A caller-specified subject cluster must use the cluster-robust variance.
clustered <- check_against_coxph(repeated, update(f, . ~ . + cluster(id)),
                                 repeated$w, "robust se")
unclustered <- check_against_coxph(repeated, f, repeated$w, "robust se")
stopifnot(isTRUE(all.equal(clustered[, "Estimate"], unclustered[, "Estimate"])),
          any(abs(clustered[, "Std.Error"] - unclustered[, "Std.Error"]) > 1e-5))
