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
#' Obsdata_long = gendata_trajmsm(n = 2000, format = "long", seed = 945)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' time_values <- c(2011,2012,2013)
#' formulaA = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = Obsdata_long, number_traj = 3, formula = formulaA, identifier = "id")
#' datapost = restraj$data_post
#' trajmsm_long <- merge(Obsdata_long, datapost, by = "id")
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
#' res_pltmle = pltmle(formula = formulaY, outcome = outcome,treatment = treatment_var,
#'                   covariates = covar, baseline = baseline_var, ntimes_interval = 3, number_traj = 3,
#'                  time =  "Time",time_values = time_values,identifier = "id",obsdata = trajmsm_wide,traj=traj_indic, treshold = 0.99)
#' res_pltmle$counter_means
#' @author Awa Diop, Denis Talbot


pltmle <- function(formula, outcome, treatment, covariates, baseline, ntimes_interval, number_traj,
                   time, time_values, identifier, obsdata, traj, total_followup, treshold = treshold) {
  # Initialize variables
  D = NULL #Influence curve
  D_list = list() #To store all influence curves
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
                                          total_followup = total_followup, numerator = "unstabilized",
                                          include_censor = include_censor, censor = censor,obsdata = obsdata)[[1]];

  weights_trunc <- sapply(1:ntimes_interval, function(x){
  weights <- ifelse(quantile(Weights[, x], treshold, na.rm = TRUE)> Weights[, x], quantile(Weights[, x], treshold, na.rm = TRUE), Weights[, x])
  return(weights)
  })
  # Compute Hs for each regime
  Hs.all = lapply(1:nregimes, function(x) {
    regime = list.regimes[[x]]
    Hs = as.matrix(rowSums(sapply(1:ntimes_interval, function(i) obsdata[, treatment[i]] == regime[i])) == ntimes_interval) *  weights_trunc[, ntimes_interval]
    return(Hs %*% t(traj[x, 1:number_traj]))
  })

  Hs = do.call(rbind, Hs.all)

  # Update the risk for each regime of treatment
  modEs = glm(as.formula(paste0(outcome, "~", "-1 + offset(qlogis(unlist(Qs))) + Hs")), family = binomial, data = obsdata.all2);;

  coef_Es = ifelse(is.na(coef(modEs)), 0, coef(modEs))

  Qstar = lapply(1:nregimes,function(x)plogis(qlogis(Qs[[x]]) +
                                                as.numeric(t(t(as.matrix( coef_Es))%*%t(Hs.all[[x]])))));

  # Influence curve for each treatment regime
  Ds = lapply(1:nregimes,function(x)sapply(1:number_traj,function(y)
    (Hs.all[[x]][,y]*(obsdata[,outcome] - Qstar[[x]]))))

  D_list[ntimes_interval]<-list(Ds)

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
      Hs = as.matrix(rowSums(sapply(1:i, function(i) obsdata[, treatment[i]] == regime[i])) == i) * weights_trunc[, i]
      return(Hs %*% t(traj[x, 1:number_traj]))
    })

    Hs = do.call(rbind, Hs.all)

    modEs = glm(unlist(Qstar) ~ -1 + offset(qlogis(unlist(Qs))) + Hs, family = binomial(), data = obsdata.all2)
    coef_Es = ifelse(is.na(coef(modEs)), 0, coef(modEs))
    Qstarm1 = lapply(1:nregimes,function(x)plogis(qlogis(Qs[[x]]) +
                                                  as.numeric(t(t(as.matrix( coef_Es))%*%t(Hs.all[[x]])))));
    Ds = lapply(1:nregimes,function(x) sapply(1:number_traj,function(y)
      as.matrix(Hs.all[[x]][,y]*(Qstar[[x]]-Qstarm1[[x]]))))
    D_list[i]<-list(Ds)
    Qstar = Qstarm1;
  }



  # Aggregate the results
  Q0 = unlist(Qstar)
  obsdata0.all2$Y = Q0
  obsdata0.all2$atraj = apply(obsdata0.all2[, treatment], 1, paste0, collapse = "")
  obsdataT2 = aggregate(Y ~ atraj, data = obsdata0.all2, FUN = mean)
  obsdataT2$Y = as.numeric(as.character(obsdataT2$Y))

  #Influence curves
  D = lapply(1:nregimes,function(x) as.matrix(Reduce('+',lapply(1:number_traj, function(y){
    comps = as.matrix(D_list[[y]][[x]] +  as.numeric(Qstar[[x]] - mean(Qstar[[x]])))
    return(comps)
  }))))
  list_pltmle_countermeans = list(counter_means = obsdataT2$Y, D = D)
  return(list_pltmle_countermeans)
}


