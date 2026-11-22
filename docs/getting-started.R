# A single-cohort LCGA-MSM example with stabilized treatment weights.
# Run from the repository root: Rscript docs/getting-started.R
library(trajmsm)

years <- 2011:2014
long <- gendata(n = 800, start_year = years[1], total_followup = length(years),
                format = "long", seed = 845)
stopifnot(!anyDuplicated(long[c("id", "time")]))

trajectory <- build_traj(
  obsdata = long, formula = cbind(statins, 1 - statins) ~ I(time - 2011),
  number_traj = 2, identifier = "id", seed = 945
)
classes <- trajectory$data_post
stopifnot(!anyDuplicated(classes$id), nlevels(droplevels(classes$class)) == 2L)
print(table(classes$class))
print(BIC(trajectory$traj_model))

wide <- reshape(long[c("id", "time", "age", "sex", "y", "statins", "hyper", "bmi")],
                direction = "wide", idvar = "id", timevar = "time",
                v.names = c("statins", "hyper", "bmi"), sep = "")
class_row <- match(wide$id, classes$id)
stopifnot(!anyDuplicated(wide$id), !anyNA(class_row))
wide$trajectory <- droplevels(classes$class[class_row])

# The numeric class labels are arbitrary. Use the group with lower observed
# average adherence as the reference in this illustration.
observed_adherence <- rowMeans(wide[paste0("statins", years)])
adherence_by_class <- tapply(observed_adherence, wide$trajectory, mean)
wide$trajectory <- relevel(wide$trajectory, ref = names(which.min(adherence_by_class)))

# Optional plot: observed averages within allocated classes, not model predictions.
if (interactive() || "--plot" %in% commandArgs(trailingOnly = TRUE)) {
  write_pdf <- !interactive()
  if (write_pdf) pdf("observed-trajectories.pdf", width = 7, height = 4.5)
  group_by_visit <- with(long, tapply(statins, list(time, wide$trajectory[match(id, wide$id)]), mean))
  matplot(years, group_by_visit, type = "b", lty = 1:2, pch = c(1, 19),
          col = "black", xaxt = "n", ylim = c(0, 1), xlab = "Year",
          ylab = "Observed proportion treated")
  axis(1, at = years)
  legend("bottomright", legend = paste("Class", colnames(group_by_visit)),
         lty = 1:2, pch = c(1, 19), col = "black", bty = "n")
  if (write_pdf) dev.off()
}

treatment <- paste0("statins", years)
covariates <- lapply(years, function(year) paste0(c("hyper", "bmi"), year))
weights_by_time <- inverse_probability_weighting(
  numerator = "stabilized", identifier = "id", baseline = c("age", "sex"),
  treatment = treatment, covariates = covariates, include_censor = FALSE,
  obsdata = wide
)[[1]]
stopifnot(nrow(weights_by_time) == nrow(wide), ncol(weights_by_time) == length(years))
w <- weights_by_time[, length(years)]
stopifnot(all(is.finite(w)), all(w > 0))

ess <- function(weight) sum(weight)^2 / sum(weight^2)
print(adherence_by_class)
print(quantile(w, c(0, .01, .50, .95, .99, 1)))
print(c(n = nrow(wide), effective_n = ess(w)))
print(tapply(w, wide$trajectory, ess))

# Save the weights against unique IDs so an analysis can reindex them if needed.
weights_by_id <- setNames(w, wide$id)
wide$iptw <- unname(weights_by_id[as.character(wide$id)])
stopifnot(identical(wide$iptw, w))

# Use the diagnosed weights unchanged for the primary illustration.
fit <- glm(y ~ trajectory, family = quasibinomial(), data = wide, weights = iptw)
coefficient_table <- function(model) {
  estimate <- coef(model)
  se <- sqrt(diag(sandwich::vcovCL(model, cluster = ~ id)))
  critical <- qnorm(.975)
  cbind(Estimate = estimate, SE = se,
        `Lower CI` = estimate - critical * se,
        `Upper CI` = estimate + critical * se)
}
estimates <- coefficient_table(fit)
stopifnot(nrow(estimates) == 2L, all(is.finite(estimates)))
print(estimates)

# Coefficients and confidence limits are log odds. Do not exponentiate SEs.
print(exp(estimates[grepl("^trajectory", rownames(estimates)),
              c("Estimate", "Lower CI", "Upper CI"), drop = FALSE]))

# Optional illustrative sensitivity fit: cap only the upper tail.
# Choose and justify any cap for a real study before interpreting its estimate.
cap <- unname(quantile(w, .99))
w_capped <- pmin(w, cap)
stopifnot(all(w_capped <= w), all(w_capped > 0))
print(c(n_capped = sum(w > cap), effective_n_capped = ess(w_capped)))
wide$iptw_capped <- pmin(wide$iptw, cap)
fit_capped <- glm(y ~ trajectory, family = quasibinomial(), data = wide, weights = iptw_capped)
estimates_capped <- coefficient_table(fit_capped)
stopifnot(identical(rownames(estimates), rownames(estimates_capped)), all(is.finite(estimates_capped)))
print(estimates_capped)

# Versions used to produce the example.
print(sessionInfo())
