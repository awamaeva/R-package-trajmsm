#' @title Pooled LTMLE
#' @description Estimate parameters of LCGA-MSM using pooled LTMLE
#'  with influence functions to estimate standard errors.
#' @name trajmsm_pltmle
#' @param formula Specification of the model for the outcome to be fitted.
#' @param identifier  Name of the column for unique identifiant.
#' @param baseline Names of the baseline covariates.
#' @param covariates Names of the time-varying covariates (should be a list).
#' @param treatment Name of the time-varying treatment.
#' @param outcome Name of the outcome variable.
#' @param time Name of the time variable.
#' @param time_values Measuring times.
#' @param total_followup Total length of follow-up.
#' @param number_traj An integer to choose the number of trajectory groups.
#' @param trajmodel Trajectory model built with the observed treatment.
#' @param treshold For weight truncation.
#' @param obsdata Observed data in wide format.
#' @param ref The reference group.
#' @param class_var Name of the trajectory group variable.
#' @importFrom stats na.omit rbinom plogis qlogis  reshape glm
#' binomial coef as.formula ave aggregate relevel pnorm sd quantile model.matrix
#' quasibinomial var
#' @importFrom utils combn
#' @return Provides a matrix of estimates for LCGA-MSM, obtained using the pooled ltlmle method.
#' @export trajmsm_pltmle
#' @examples
#' \donttest{
#' obsdata_long = gendata(n = 1000, format = "long", total_followup = 6, seed = 845)
#' years <- 2011:2016
#' baseline_var <- c("age","sex")
#' variables <- c("hyper", "bmi")
#' covariates <- lapply(years, function(year) {
#' paste0(variables, year)})
#' treatment_var <- paste0("statins", 2011:2016)
#' formula_treatment = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = obsdata_long, number_traj = 3,
#' formula = formula_treatment, identifier = "id")
#' datapost = restraj$data_post
#' trajmsm_long <- merge(obsdata_long, datapost, by = "id")
#'     AggFormula <- as.formula(paste("statins", "~", "time", "+", "class"))
#'     AggTrajData <- aggregate(AggFormula, data = trajmsm_long, FUN = mean)
#' trajmsm_wide = reshape(data = trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper"), timevar = "time", sep ="")
#'formula = paste0("y ~", paste0(treatment_var,collapse = "+"), "+",
#'                 paste0(unlist(covariates), collapse = "+"),"+",
#'                 paste0(baseline_var, collapse = "+"))
#' resmsm_pltmle <- trajmsm_pltmle(formula = formula, identifier = "id",
#'  baseline = baseline_var,
#'  covariates = covariates, treatment = treatment_var,
#'  outcome = "y", time = "time", time_values = years,
#'  number_traj = 3, total_followup = 6,
#'  trajmodel = restraj$traj_model, ref = "1", obsdata = trajmsm_wide,
#'   treshold = 1,class_var = "class")
#'  resmsm_pltmle
#'  }
#' @return \item{results_msm_pooledltmle}{Estimates of a LCGA-MSM with pooled LTMLE.}
#' @author Awa Diop, Denis Talbot


trajmsm_pltmle <- function(formula = formula,identifier,baseline,covariates,treatment,outcome,
                                 number_traj,total_followup, time, time_values , trajmodel,ref,
                                treshold = 0.99, obsdata, class_var){

 stopifnot(!is.null(identifier));
 stopifnot(!is.null(baseline));
 stopifnot(!is.null(covariates));
 stopifnot(!is.null(treatment));
 stopifnot(!is.null(outcome));
 stopifnot(!is.null(total_followup));
 stopifnot(!is.null(obsdata));
 stopifnot(!is.null(trajmodel));
 stopifnot(!is.null(time));
     nregimes = 2^total_followup;  #number of treatment regimes

     treatment_names <- sub("\\d+", "", treatment)
     treatment_name <- unique(treatment_names)[1]
     class = factor(predict_traj(identifier = identifier, total_followup = total_followup,
                                 treatment = treatment_name, time = time, time_values = time_values,
                                 trajmodel = trajmodel)$post_class);
     if(length(unique(class)) < number_traj){
       stop(paste("Number of predicted trajectory groups inferior to", number_traj))
     }

     if(length(unique(class)) == number_traj){
     traj_indic=t(sapply(1:nregimes,function(x)sapply(1:number_traj,function(i) ifelse(class[x]==i,1,0))))
     traj_indic[,1]=1 #Intercept
    #Create the obsdata under all the different regime of treatment
     res_pltmle = pltmle(formula = formula, outcome = outcome,treatment = treatment,
                 covariates = covariates, baseline = baseline, ntimes_interval = total_followup, number_traj = number_traj,
                 time =  time, identifier = identifier,obsdata = obsdata,
                 traj=traj_indic, treshold = treshold, class_pred = class, class_var = class_var);

    obsdata_pool= data.frame(Y=res_pltmle$counter_means);
    D=res_pltmle$D; #Influence functions

    # Estimation
    obsdata_pool$ptmle_group = class
    obsdata_pool$ptmle_group <- relevel(as.factor(obsdata_pool$ptmle_group), ref = ref)
    mod = glm(Y ~ ptmle_group, family = quasibinomial, data = obsdata_pool);
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

}
