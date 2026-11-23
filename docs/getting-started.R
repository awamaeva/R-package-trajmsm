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

# Optional figure: observed group proportions, not fitted probabilities.
if (interactive() || "--plot" %in% commandArgs(trailingOnly = TRUE)) {
  # Join by ID; use the same relevelled group labels as the outcome analysis.
  plot_data <- transform(long, group = wide$trajectory[match(id, wide$id)])
  stopifnot(!anyNA(plot_data$group))
  # ggtraj plots supplied values, so compute the proportions first.
  visit_means <- aggregate(statins ~ time + group, data = plot_data, FUN = mean)
  p <- trajmsm::ggtraj(
    traj_data = visit_means, treatment = "statins", time = "time",
    identifier = "id", class = "group"
  ) +
    ggplot2::scale_color_manual(values = c("#193b68", "#984818")) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, .25),
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::labs(title = NULL, x = "Visit year", y = "Proportion treated",
                  color = "Group", shape = "Group", linetype = "Group") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (interactive()) {
    print(p)
  } else {
    grDevices::svg("observed-trajectories.svg", width = 5.2, height = 4.1,
                   bg = "white")
    print(p)
    grDevices::dev.off()
    ggplot2::ggsave("observed-trajectories.pdf", p, width = 5.2, height = 4.1)
  }
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

# Compare the package wrapper with the same explicit fit when the lookup fix
# is installed. Older versions still run the direct fit above.
if (packageVersion("trajmsm") >= "0.1.7") {
  wrapped <- trajmsm_ipw(
    formula1 = y ~ trajectory, family = "quasibinomial", identifier = "id",
    obsdata = wide, weights = wide$iptw
  )
  direct_se <- sqrt(diag(sandwich::vcovCL(fit, cluster = ~ id)))
  stopifnot(isTRUE(all.equal(wrapped[, "Estimate"], coef(fit))),
            isTRUE(all.equal(wrapped[, "Std.Error"], direct_se)))
  message("trajmsm_ipw() agrees with the direct weighted fit.")
}

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
