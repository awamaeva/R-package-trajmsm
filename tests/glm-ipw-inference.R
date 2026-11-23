library(trajmsm)

# Binary group log odds differ by log(2) when weighted and zero otherwise.
# The Gaussian weighted group difference is 4/5 - 2/3 instead of zero.
wide <- data.frame(id = rep(101:104, 2L), y = rep(c(0, 0, 1, 1), 2L),
                   group = rep(0:1, each = 4L))
w <- c(1, 1, 1, 3, 1, 1, 1, 7)

glm_reference <- function(data, formula, family, weight) {
  fit_data <- transform(data, reference_weight = weight)
  direct <- glm(formula, family = family, data = fit_data, weights = reference_weight)
  retained <- match(rownames(model.frame(direct)), rownames(fit_data))
  estimate <- coef(direct)
  se <- sqrt(diag(sandwich::vcovCL(direct, cluster = fit_data$id[retained])))
  expected <- cbind(estimate, se, 2 * pnorm(-abs(estimate) / se),
                    estimate - 1.96 * se, estimate + 1.96 * se)
  colnames(expected) <- c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")
  expected
}

check_fit <- function(data, weight, family, shadow = FALSE, complete = TRUE) {
  # A formula can call functions defined by the user's own environment.
  user_transform <- function(x) x
  f <- y ~ user_transform(group)
  if (shadow) {
    # The unrelated caller data have different weights and ID values.
    obsdata <- transform(data, weights = 1, id = seq_len(nrow(data)))
  }
  got <- trajmsm_ipw(formula1 = f, family = family, identifier = "id",
                    obsdata = data, weights = weight)
  expected <- glm_reference(data, f, family, weight)
  stopifnot(isTRUE(all.equal(got, expected, tolerance = 1e-8)))
  if (complete) {
    group_effect <- if (family == "binomial") log(2) else 2/15
    stopifnot(abs(got[2, "Estimate"] - group_effect) < 1e-8)
  }
}

for (family in c("binomial", "gaussian")) {
  check_fit(wide, w, family)
  check_fit(wide, w, family, shadow = TRUE)

  # Missing model values and weights must omit matching rows from clustering.
  incomplete <- wide
  incomplete$y[1] <- NA
  incomplete_w <- replace(w, 5L, NA)
  for (na_action in c("na.omit", "na.exclude")) {
    prior <- options(na.action = na_action)
    check_fit(incomplete, incomplete_w, family, complete = FALSE)
    options(prior)
  }
}

# Do not let a missing column silently switch to per-row robust inference.
bad_id <- tryCatch(
  trajmsm_ipw(formula1 = y ~ group, family = "gaussian", identifier = "i",
              obsdata = wide, weights = w),
  error = conditionMessage
)
stopifnot(identical(bad_id, "'identifier' must name a column in 'obsdata'."))

# Check computed weights against a direct fit using the current flooring rule.
sim <- gendata(n = 800, start_year = 2011, total_followup = 3,
               format = "wide", seed = 845)
treatment <- paste0("statins", 2011:2013)
covariates <- lapply(2011:2013, function(year) paste0(c("hyper", "bmi"), year))
baseline <- c("age", "sex")
threshold <- 0.75
computed <- inverse_probability_weighting(
  identifier = "id", baseline = baseline, covariates = covariates,
  treatment = treatment, obsdata = sim, numerator = "stabilized"
)[[1]][, length(treatment)]
floored <- pmax(computed, quantile(computed, threshold, na.rm = TRUE))
stopifnot(any(computed < floored), any(computed == floored))
actual <- trajmsm_ipw(formula1 = y ~ statins2013, family = "gaussian",
                      identifier = "id", baseline = baseline, covariates = covariates,
                      treatment = treatment, obsdata = sim, treshold = threshold)
expected <- glm_reference(sim, y ~ statins2013, "gaussian", floored)
stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-8)))
