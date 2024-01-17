#' @title Pooled LTMLE
#' @description function to estimate parameters of a MSM-LCGA using pooled LTMLE
#'  with influence functions to estimate standard errors.
#' @name trajmsm_pltmle
#' @param formula specification of the model for the outcome to be fitted.
#' @param identifier  name of the column for unique identifiant.
#' @param covariates covariates.
#' @param treatment time-varying treatment.
#' @param time name of the time variable.
#' @param total_followup number of measuring times per interval.
#' @param number_traj an integer to choose the number of trajectory groups.
#' @param trajmodel trajectory model built with the observed treatment.
#' @param obsdata observed data in wide format.
#' @param ref the reference group.
#' @export trajmsm_pltmle
#' @examples
#' Obsdata_long = gendata_trajmsm(n = 2000, format = "long", seed = 945)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' time = "Time"
#' time_values <- c(2011,2012,2013)
#' formulaA = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = Obsdata_long, number_traj = 3, formula = formulaA, identifier = "id")
#' datapost = restraj$data_post
#' trajmsm_long <- merge(Obsdata_long, datapost, by = "id")
#'     AggFormula <- as.formula(paste("statins", "~", "time", "+", "class"))
#'     AggTrajData <- aggregate(AggFormula, data = trajmsm_long, FUN = mean)
#'     AggTrajData
#' trajmsm_wide = reshape(data = trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper"), timevar = "time", sep ="")
#' formulaY =  as.formula(" y ~ statins2011 + statins2012 + statins2013 + hyper2011 + bmi2011 + hyper2012 + bmi2012 +
#'                                     hyper2013 + bmi2013 + age + sex ")
#'trajmsm_pltmle(formula = formulaY, identifier = identifier, baseline = baseline,
#'               covariates = covariates, treatment = treatment, outcome = "y",
#'              time = "time", time_values = time_values, number_traj = 3, total_followup = 3,
#'              trajmodel = restraj$traj_model, ref = "3", obsdata = trajmsm_wide, treshold = 0.99)
#' @return \item{results_msm_pooledltmle}{Estimates of a LCGA-MSM with pooled LTMLE.}
#' @author Awa Diop, Denis Talbot


trajmsm_pltmle <- function(formula = formula,identifier,baseline,covariates,treatment,outcome,
                                 number_traj,total_followup, time,time_values,trajmodel,ref,obsdata,treshold){

 stopifnot(!is.null(identifier));
 stopifnot(!is.null(baseline));
 stopifnot(!is.null(covariates));
 stopifnot(!is.null(treatment));
 stopifnot(!is.null(outcome));
 stopifnot(!is.null(total_followup));
 stopifnot(!is.null(obsdata));
 stopifnot(!is.null(trajmodel));
 stopifnot(!is.null(time));
 stopifnot(!is.null(time_values));
  stopifnot(!is.null(treshold));
     nregimes = 2^total_followup;  #number of treatment regimes

     treatment_names <- sub("\\d+", "", treatment)
     treatment_name <- unique(treatment_names)[1]
     class = factor(predict_traj(identifier = identifier, total_followup = total_followup,
                                 treatment = treatment_name, time = time, time_values = time_values,
                                 trajmodel = trajmodel)$post_class);

     traj_indic=t(sapply(1:nregimes,function(x)sapply(1:number_traj,function(i) ifelse(class[x]==i,1,0))))
     traj_indic[,1]=1 #Intercept
    #Create the obsdata under all the different regime of treatment
     res_pltmle = pltmle(formula = formula, outcome = outcome,treatment = treatment,
                 covariates = covariates, baseline = baseline, ntimes_interval = total_followup, number_traj = number_traj,
                 time =  time, time_values = time_values,identifier = identifier,obsdata = obsdata,traj=traj_indic, treshold = treshold);

    obsdata_pool= data.frame(Y=res_pltmle$counter_means);
    D=res$D; #Influence functions

    # Estimation
    obsdata_pool$TMLE_group = class
    obsdata_pool$TMLE_group <- relevel(as.factor(obsdata_pool$TMLE_group), ref = ref)
    mod = glm(Y~TMLE_group, family = quasibinomial, data = obsdata_pool);
    coefs = summary(mod)$coefficients[,1];
    X = model.matrix(mod)[1:nregimes,]
    B = matrix(coefs, nrow = number_traj);

    CQ=lapply(1:nregimes,function(i)(as.matrix(X[i,]))%*%as.matrix((plogis(as.matrix(t(X[i,]))%*%B))/(1+exp(as.matrix(t(X[i,]))%*%B))^2)%*%t(as.matrix(((X[i,])))));
    CQ=Reduce('+',CQ)

    Db = matrix(0, nrow = nrow(D[[1]]), ncol = number_traj);
    for(l in 1:nregimes){
      Db = Db+as.matrix(D[[l]])%*%solve(CQ);
    }

    se = sqrt(diag(var(Db)/nrow(Db)));
    pvalue <- 2*pnorm(-abs(coefs)/se)
    IClo = coefs - 1.96*se ;
    ICup = coefs + 1.96*se

    results_msm_pltmle = cbind(coefs,se,pvalue,IClo, ICup);
    colnames(results_msm_pltmle) = c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")
    return(results_msm_pltmle)
  }

