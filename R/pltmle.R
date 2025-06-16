#' @title Counterfactual means for a Pooled LTMLE
#' @description Function to estimate counterfactual means for a pooled LTMLE.
#' @name pltmle
#' @param formula Specification of the model for the outcome to be fitted.
#' @param identifier  Name of the column of the unique identifier.
#' @param baseline Name of baseline covariates.
#' @param covariates Covariates.
#' @param treatment Time-varying treatment.
#' @param outcome Name of the outcome variable.
#' @param time Name of the time variable.
#' @param time_values Measuring times.
#' @param total_followup Number of measuring times per interval.
#' @param ntimes_interval Length of a time-interval (s).
#' @param traj Matrix of indicators for the trajectory groups.
#' @param number_traj An integer to choose the number of trajectory groups.
#' @param treshold For weight truncation.
#' @param class_var Name of the trajectory group variable.
#' @param class_pred Vector of predicted trajectory groups.
#' @param obsdata Observed data in wide format.
#' @returns  \item{list_pltmle_countermeans}{ Counterfactual means and influence functions with the pooled ltmle.}
#' \item{D}{Influence functions}
#' @import e1071
#' @export
#' @examples
#' \donttest{
#' obsdata_long = gendata(n = 2000, format = "long",total_followup = 3, seed = 945)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),
#' c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' time_values <- c(2011,2012,2013)
#' formulaA = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = obsdata_long, number_traj = 3,
#' formula = formulaA, identifier = "id")
#' datapost = restraj$data_post
#' trajmsm_long <- merge(obsdata_long, datapost, by = "id")
#'     AggFormula <- as.formula(paste("statins", "~", "time", "+", "class"))
#'     AggTrajData <- aggregate(AggFormula, data = trajmsm_long, FUN = mean)
#'     AggTrajData
#' trajmsm_long[ , "traj_group"] <- trajmsm_long[ , "class"]
#' obsdata= reshape(trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper"), timevar = "time", sep ="")
#' formula =  as.formula(" y ~ statins2011 + statins2012 + statins2013 +
#' hyper2011 + bmi2011 + hyper2012 + bmi2012 +
#'  hyper2013 + bmi2013 + age + sex ")
#' class = factor(predict_traj(identifier = "id", total_followup = 3,
#'         treatment = "statins", time = "time", time_values = time_values,
#'         trajmodel = restraj$traj_model)$post_class);
#' traj=t(sapply(1:8,function(x)sapply(1:3,function(i)ifelse(class[x]==i,1,0))))
#' traj[,1]=1
#' res_pltmle = pltmle(formula = formula, outcome = "y",treatment = treatment_var,
#' covariates = covariates, baseline = baseline_var, ntimes_interval = 3, number_traj = 3,
#'  time =  "time",time_values = time_values,identifier = "id",obsdata = obsdata,
#' traj=traj, treshold = 0.99, class_pred= class, class_var = "class")
#' res_pltmle$counter_means
#' }
#' @author Awa Diop, Denis Talbot

pltmle <- function (formula, outcome, treatment, covariates, baseline,
                    ntimes_interval, number_traj, time, time_values, identifier,
                    obsdata, traj, total_followup, treshold = treshold, class_var, class_pred)
{

  D = NULL
  D_list = list()
  obsdata0.all = list()
  obsdata.all = list()
  nregimes = 2^ntimes_interval
  dat_combn = bincombinations(ntimes_interval)
  list.regimes = lapply(1:nregimes, function(x) {
    regime = dat_combn[x, ]
    return(regime)
  })

  for (regime in list.regimes) {
    obsdata0 = obsdata
    obsdata0[, treatment[1:ntimes_interval]] = sapply(regime,
                                                      function(r) {
                                                        return(rep(r, nrow(obsdata0)))
                                                      })
    obsdata0.all[[length(obsdata0.all) + 1]] = obsdata0
    obsdata.all[[length(obsdata.all) + 1]] = obsdata
  }

  obsdata.all2 = do.call(rbind, obsdata.all)
  obsdata0.all2 = do.call(rbind, obsdata0.all)

  obsdata.all.new <- do.call(rbind,replicate(number_traj,obsdata.all2, simplify = FALSE))
  modQs = glm(formula, family = binomial, data = obsdata)

  Qs = lapply(1:nregimes, function(i) predict(modQs, newdata = obsdata0.all[[i]],
                                              type = "response"))

  Weights = inverse_probability_weighting(identifier = identifier,
                                          covariates = covariates, treatment = treatment, baseline = baseline,
                                          numerator = "unstabilized", include_censor = FALSE, obsdata = obsdata)[[1]]

  weights_trunc <- Weights

  for (j in 1:ncol(Weights)) {
    quantile_value <- na.omit(quantile(Weights[, j], treshold,
                                       na.rm = TRUE))
    weights_trunc[, j] <- ifelse(Weights[, j] > quantile_value, quantile_value,
                                 Weights[, j])
  }


  Hs.all = lapply(1:nregimes, function(x) {
    regime = list.regimes[[x]]
    Hs_temp = as.matrix(rowSums(sapply(1:ntimes_interval, function(i) obsdata[, treatment[i]] == regime[i])) == ntimes_interval) *
      weights_trunc[, ntimes_interval]
    return(Hs_temp %*% t(traj[x, 1:number_traj]))
  })

  # All Hs
  obsdata.all.new$Hs <- as.vector(unlist(do.call(rbind, Hs.all)))

  # Predicted trajectory groups
  obsdata.all.new$class.all <- factor(do.call(rbind,lapply(1:nregimes, function(x)
    rep(class_pred[x], nrow(obsdata)))))

  # Repeat unlist(Qs) x times
  obsdata.all.new$offset_term <- rep(qlogis(unlist(Qs)), number_traj)

  # Simultaneous estimation of all errors
  coef_Es <- c()
  modEs <- glm(as.formula(paste0(outcome, "~ -1 + class.all")),
               weights = obsdata.all.new$Hs,
               offset =  obsdata.all.new$offset_term,
               family = quasibinomial,
               data = obsdata.all.new)

  # Get the coefficients
  coef_Es <- coef(modEs)

  traj_ind <- sapply(1:number_traj, function(i) {
    ifelse(obsdata[, class_var] == i, 1, 0)
  })

  traj_ind[,1] <- 1

  Qstar = lapply(1:nregimes, function(x) plogis(qlogis(Qs[[x]]) +
                                                  as.numeric(t(t(as.matrix(coef_Es)) %*% t(traj_ind)))))

  # Influence curve for each treatment regime
  Ds = lapply(1:nregimes,function(x)sapply(1:number_traj,function(y)
    (Hs.all[[x]][,y]*(obsdata[,outcome] - Qstar[[x]]))))

  D_list[ntimes_interval]<-list(Ds)

  for (i in (ntimes_interval - 1):1) {
    current_treatments = paste(treatment[1:i], collapse = "+")
    current_covariates = paste(unlist(covariates[1:i]), collapse = "+")
    baseline_covariates = paste(baseline, collapse = "+")
    modQs = lapply(1:nregimes, function(x) {
      Qstar_current = Qstar[[x]]
      formula_update = as.formula(paste("Qstar_current ~",
                                        current_treatments, "+", current_covariates,
                                        "+", baseline_covariates))
      glm(formula_update, family = binomial, data = obsdata)
    })
    Qs = lapply(1:nregimes, function(x) predict(modQs[[x]],
                                                newdata = obsdata0.all[[x]], type = "response"))
    Hs.all = lapply(1:nregimes, function(x) {
      regime = list.regimes[[x]]
      Hs = as.matrix(rowSums(sapply(1:i, function(i) obsdata[,treatment[i]] == regime[i])) == i) * weights_trunc[, i]
      return(Hs %*% t(traj[x, 1:number_traj]))
    })
    obsdata.all.new$Hs <- as.vector(unlist(do.call(rbind, Hs.all)))

    current_treatments = paste(treatment[1:i], collapse = "+")
    current_covariates = paste(unlist(covariates[1:i]), collapse = "+")
    baseline_covariates = paste(baseline, collapse = "+")
    modQs = lapply(1:nregimes, function(x) {
      Qstar_current = Qstar[[x]]
      formula_update = as.formula(paste("Qstar_current ~",
                                        current_treatments, "+", current_covariates,
                                        "+", baseline_covariates))
      glm(formula_update, family = binomial, data = obsdata)
    })
    Qs = lapply(1:nregimes, function(x) predict(modQs[[x]],
                                                newdata = obsdata0.all[[x]], type = "response"))
    Hs.all = lapply(1:nregimes, function(x) {
      regime = list.regimes[[x]]
      Hs_temp = as.matrix(rowSums(sapply(1:i, function(i) obsdata[,treatment[i]] == regime[i])) == i) * weights_trunc[, i]
      return(Hs_temp %*% t(traj[x, 1:number_traj]))
    })

    obsdata.all.new$Hs <- as.vector(unlist(do.call(rbind, Hs.all)))
    obsdata.all.new$outcome_repeat <- rep(unlist(Qstar), number_traj)
    obsdata.all.new$offset_term <- rep(qlogis(unlist(Qs)), number_traj)

    # Modify the formula to include the updated offset term
    coef_Es <-c()
    modEs <- glm(as.formula(paste0("outcome_repeat ~ -1 + class.all")),
                 weights = obsdata.all.new$Hs,
                 offset = obsdata.all.new$offset_term,
                 family = quasibinomial,
                 data = obsdata.all.new)

    # Get the coefficients
    coef_Es <- coef(modEs)

    Qstarm1 = lapply(1:nregimes, function(x) plogis(qlogis(Qs[[x]]) +
                                                      as.numeric(t(t(as.matrix(coef_Es)) %*% t(traj_ind)))))
    Ds = lapply(1:nregimes, function(x) sapply(1:number_traj,
                                               function(y) as.matrix(Hs.all[[x]][, y] * (Qstar[[x]] -
                                                                                           Qstarm1[[x]]))))
    D_list[i] <- list(Ds)
    Qstar = Qstarm1
  }
  Q0 = unlist(Qstar)
  obsdata0.all2$Y = Q0
  obsdata0.all2$atraj = apply(obsdata0.all2[, treatment], 1,
                              paste0, collapse = "")
  obsdataT2 = aggregate(Y ~ atraj, data = obsdata0.all2, FUN = mean)
  obsdataT2$Y = as.numeric(as.character(obsdataT2$Y))
  D = lapply(1:nregimes, function(x) as.matrix(Reduce("+",
                                                      lapply(1:number_traj, function(y) {
                                                        comps = as.matrix(D_list[[y]][[x]] + as.numeric(Qstar[[x]] -
                                                                                                          mean(Qstar[[x]])))
                                                        return(comps)
                                                      }))))

  list_pltmle_countermeans = list(counter_means = obsdataT2$Y,
                                  D = D)
  return(list_pltmle_countermeans)
}
