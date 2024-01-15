#' @title Parametric g-formula
#' @description function to estimate parameters of a MSM-LCGA using g-formula
#'  and bootstrap to get standard errors.
#' @name trajmsm_gform
#' @param formula specification of the model for the outcome to be fitted.
#' @param identifier name of the column for unique identifiant.
#' @param baseline  vector of names of the baseline covariates.
#' @param covariates list of names of the time-varying covariates.
#' @param treatment vector of names of the time-varying treatment.
#' @param outcome name of the outcome of interest.
#' @param total_followup of measuring times.
#' @param time name of the variable time.
#' @param time_values values of the time variable.
#' @param rep number of repetitions for the bootstrap.
#' @param trajmodel trajectory model built with the observed treatment.
#' @param ref the reference trajectory group.
#' @param obsdata observed data in wide format.
#' @return \item{results_msm_gform}{Estimates of a LCGA-MSM with g-formula.}
#' @export trajmsm_gform
#' @examples
#' Obsdata_long = gendata_trajmsm(n = 1000, include_censor = TRUE, format = "long", seed = 945)
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
#'trajmsm_long[ , "traj_group"] <- trajmsm_long[ , "class"]
#' trajmsm_wide = reshape(trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper","censor"), timevar = "time", sep ="")
#' formulaY =  as.formula(" y ~ statins2011 + statins2012 + statins2013 + hyper2011 + bmi2011 + hyper2012 + bmi2012 +
#'                                     hyper2013 + bmi2013 + age + sex ")
#'trajmsm_gform(formula = formulaY, identifier = "id", baseline = baseline, covariates = covariates,
#'                                     treatment = treatment_var, outcome = "y", total_followup = 3,time = time,
#'                                     time_values = time_values, trajmodel = restraj$traj_model,
#'                                     ref = "3", obsdata = trajmsm_wide)
#' @author Awa Diop Denis Talbot

trajmsm_gform <- function(formula = formula, rep = 50,
                                  identifier,baseline,covariates,treatment,outcome, total_followup,time = time,
                                  time_values = time_values, trajmodel,ref, obsdata){


  stopifnot(!is.null(identifier));
  stopifnot(!is.null(baseline));
  stopifnot(!is.null(covariates));
  stopifnot(!is.null(treatment));
  stopifnot(!is.null(outcome));
  stopifnot(!is.null(total_followup));
  stopifnot(!is.null(rep));
  stopifnot(!is.null(obsdata));
  stopifnot(!is.null(trajmodel));
  stopifnot(!is.null(time));

    bootf=function(df,x=index){
      #Echantillons bootstrap
      df=obsdata[x,];
      res = gformula(formula = formula,outcome = outcome, treatment = treatment, covariates = covariates,baseline = baseline,
                     ntimes_interval = total_followup, obsdata = df)$counter_means
      colnames(res) <- "Y";

      #Counterfactual means + trajectory groups
      obsdataG = data.frame(res);
      treatment_names <- sub("\\d+", "", treatment)
      treatment_name <- unique(treatment_names)[1]
      obsdataG$gform_group = factor(predict_traj(identifier = identifier, total_followup = total_followup,
                                                 treatment = treatment_name, time = time, time_values = time_values,
                                                 trajmodel = trajmodel)$post_class);

      obsdataG$gform_group <- relevel(as.factor(obsdataG$gform_group), ref = ref)
      #Estimation
      mod = summary(glm(Y ~ gform_group, family = quasibinomial, data = obsdataG));
      return(coef(mod)[,1]);
    }

    #Bootstrap
    list_res <- list()
    for(b in 1:rep){
      index <- sample(1:nrow(obsdata), nrow(obsdata) ,replace = T)
      list_res[b] <- list(bootf(df,x = index))
      cat('Replication', b, 'of', rep,'\n')
    }

    result.coef.boot <- do.call(rbind,list_res)
    #results
    res = gformula(formula = formula,outcome = outcome, treatment = treatment, covariates = covariates,baseline = baseline,
                    ntimes_interval = total_followup, obsdata = obsdata)$counter_means
    colnames(res) <- "Y";

    #Counterfactual means + trajectory groups
    obsdataG = data.frame(res);
    obsdataG$gform_group = factor(predict_traj(identifier = identifier, total_followup = total_followup,
                                               treatment = treatment_name, time = time, time_values = time_values,
                                               trajmodel = trajmodel)$post_class);
    obsdataG$gform_group <- relevel(as.factor(obsdataG$gform_group), ref = ref)
    #Estimation
    mod = summary(glm(Y ~ gform_group, family = quasibinomial, data = obsdataG));
    coefs.mean = mod$coefficients[,1]
    se=apply(result.coef.boot,2,sd)
    pvalue <- 2*pnorm(-abs(coefs.mean)/se)
    lo.ci = coefs.mean - 1.96*se
    up.ci = coefs.mean + 1.96*se
    results_msm_gform = rbind(coefs.mean,se,pvalue,lo.ci, up.ci);
    rownames(results_msm_gform) = c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")
    return(t(results_msm_gform));
}
