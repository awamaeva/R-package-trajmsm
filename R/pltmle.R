#' @title Counterfactual means for a Pooled LTMLE
#' @description function to estimate counterfactual means for a pooled LTMLE.
#' @name sub_pltmle
#' @param formula specification of the model for the outcome to be fitted.
#' @param identifier  name of the column for unique identifiant.
#' @param covariates covariates.
#' @param treatment time-varying treatment.
#' @param time name of the time variable.
#' @param total_followup number of measuring times per interval.
#' @param ntimes_interval length of a time-interval (s).
#' @param number_traj an integer to choose the number of trajectory groups.
#' @param obsdata observed data in wide format.
#' @returns  \item{list_pltmle_countermeans}{Counterfactual means and influence functions with the pooled ltmle.}
#' \item{D}{Influence functions}
#' @import e1071
#' @examples
#' Obsdata_long = gendata_trajmsm(n = 1000, format = "long", seed = 945)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' time = "Time"
#' time_values <- c(2011,2012,2013)
#' formulaA = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = Obsdata_long, number_traj = 3, formula = formulaA, identifier = "id")
#' Datapost = restraj$data_post
#' trajmsm_long <- merge(Obsdata_long, Datapost, by = "id")
#'     AggFormula <- as.formula(paste("statins", "~", "time", "+", "class"))
#'     AggTrajData <- aggregate(AggFormula, data = trajmsm_long, FUN = mean)
#'     AggTrajData
#' trajmsm_long[ , "traj_group"] <- trajmsm_long[ , "class"]
#' trajmsm_wide = reshape(trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper"), timevar = "time", sep ="")
#' formulaY =  as.formula(" y ~ statins2011 + statins2012 + statins2013 + hyper2011 + bmi2011 + hyper2012 + bmi2012 +
#'                                     hyper2013 + bmi2013 + age + sex ")
#' class = factor(predict_traj(identifier = "id", total_followup = 3,
#'         treatment = "statins", time = "time", time_values = time_values,
#'         trajmodel = restraj$traj_model)$post_class);
#' traj_indic=t(sapply(1:nregimes,function(x)sapply(1:number_traj,function(i) ifelse(class[x]==i,1,0))))
#' traj_indic[,1]=1
#' pltmle(formula = formulaY, outcome = outcome,treatment = treatment_var,
#'                   covariates = covar, baseline = baseline_var, ntimes_interval = 3, number_traj = 3,
#'                  time =  "Time",time_values = time_values,identifier = "id",obsdata = trajmsm_wide,traj=traj_indic)
#' @author Awa Diop, Denis Talbot


pltmle <- function(formula, outcome, treatment, covariates, baseline, ntimes_interval, number_traj,
                   time, time_values, identifier, obsdata, traj, total_followup) {
  # Initialize variables
  D = NULL
  obsdata0.all = list()
  obsdata.all = list()
  nregimes = 2^ntimes_interval  # Number of treatment regimes

  # Generate treatment regimes
  dat_combn = bincombinations(ntimes_interval);
  list.regimes=lapply(1:nregimes,function(x){;
    regime = dat_combn[x,];
    return(regime);
  })

  # Loop over each regime
  for (regime in list.regimes) {
    obsdata0 = obsdata
    obsdata0[, treatment[1:ntimes_interval]] = sapply(regime, function(r) {
      return(rep(r, nrow(obsdata0)))
    })

    obsdata0.all[[length(obsdata0.all) + 1]] = obsdata0
    obsdata.all[[length(obsdata.all) + 1]] = obsdata
  }

  obsdata.all2 = do.call(rbind, obsdata.all)
  obsdata0.all2 = do.call(rbind, obsdata0.all)

  # Fit the model
  modQs = glm(formula, family = binomial(), data = obsdata)

  # Predict the outcome for all different regimes of treatment
  Qs = lapply(1:nregimes, function(i) predict(modQs, newdata = obsdata0.all[[i]], type = "response"))

  # Compute the weights for all different regimes of treatment
  Weights = inverse_probability_weighting(identifier = identifier, covariates = covariates,
                                          treatment = treatment, baseline = baseline,
                                          total_follow_up = total_followup, numerator = "unstabilized",
                                          include_censor = include_censor, censor = censor,obsdata = obsdata)[[1]];

  # Compute Hs for each regime
  Hs.all = lapply(1:nregimes, function(x) {
    regime = list.regimes[[x]]
    Hs = as.matrix(rowSums(sapply(1:ntimes_interval, function(i) obsdata[, treatment[i]] == regime[i])) == ntimes_interval) * Weights[, ntimes_interval]
    return(Hs %*% t(traj[x, 1:number_traj]))
  })

  Hs = do.call(rbind, Hs.all)

  # Update the risk for each regime of treatment
  modEs = glm(unlist(Qstar) ~-1+offset(qlogis(unlist(Qs))) + Hs, family = binomial, data = obsdata.all2);

  coef_Es = ifelse(is.na(coef(modEs)), 0, coef(modEs))

  Qstar = lapply(1:nregimes,function(x)plogis(qlogis(Qs[[x]]) +
                                                as.numeric(t(t(as.matrix( coef_Es))%*%t(Hs.all[[x]])))));

  # Influence curve for each regime of treatment
  D = lapply(1:nregimes, function(x) {
    sapply(1:number_traj, function(y) Hs.all[[x]][, y] * (obsdata[, outcome] - Qstar[[x]]))
  })

  # Loop to update the risk for each time interval
  for (i in (ntimes_interval - 1):1) {
    current_treatments = paste(treatment[1:i], collapse = "+")
    current_covariates = paste(unlist(covariates[1:i]), collapse = "+")
    baseline_covariates = paste(baseline, collapse = "+")

    modQs = lapply(1:nregimes, function(x) {
      Qstar_current = Qstar[[x]]
      formula_update = as.formula(paste("Qstar_current ~", current_treatments, "+", current_covariates, "+", baseline_covariates))
      glm(formula_update, family = binomial(), data = obsdata)
    })

    Qs = lapply(1:nregimes, function(x) predict(modQs[[x]], newdata = obsdata0.all[[x]], type = "response"))

    # Update Hs for each regime
    Hs.all = lapply(1:nregimes, function(x) {
      regime = list.regimes[[x]]
      Hs = as.matrix(rowSums(sapply(1:i, function(i) obsdata[, treatment[i]] == regime[i])) == i) * Weights[, i]
      return(Hs %*% t(traj[x, 1:number_traj]))
    })

    Hs = do.call(rbind, Hs.all)

    modEs = glm(unlist(Qstar) ~ -1 + offset(qlogis(unlist(Qs))) + Hs, family = binomial(), data = obsdata.all2)
    coef_Es = ifelse(is.na(coef(modEs)), 0, coef(modEs))
    Qstar = lapply(1:nregimes,function(x)plogis(qlogis(Qs[[x]]) +
                                                  as.numeric(t(t(as.matrix( coef_Es))%*%t(Hs.all[[x]])))));
  }

  # Aggregate the results
  Q0 = unlist(Qstar)
  obsdata0.all2$outcome = Q0
  obsdata0.all2$atraj = apply(obsdata0.all2[, treatment], 1, paste0, collapse = "")
  obsdataT2 = aggregate(outcome ~ atraj, data = obsdata0.all2, FUN = mean)
  obsdataT2$outcome = as.numeric(as.character(obsdataT2$outcome))

  list_pltmle_countermeans = list(counter.means = obsdataT2$outcome, D = D)
  return(list_pltmle_countermeans)
}


