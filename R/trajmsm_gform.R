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
#' @param rep number of repetitions for the bootstrap.
#' @param trajmodel trajectory model built with the observed treatment.
#' @param ref the reference trajectory group.
#' @param obsdata observed data in wide format.
#' @return \item{results_msm_gform}{Estimates of a LCGA-MSM with g-formula.}
#' @export trajmsm_gform
#' @examples
#' obsdata_long = gendata(n = 2000, format = "long", total_followup = 6, seed = 1945)
#' years <- 2011:2016
#' baseline_var <- c("age","sex")
#' variables <- c("hyper", "bmi")
#' covariates <- lapply(years, function(year) {
#' paste0(variables, year)})
#' treatment_var <- paste0("statins", 2011:2016)
#' formula_treatment = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = obsdata_long, number_traj = 3, formula = formula_treatment, identifier = "id")
#' datapost = restraj$data_post
#' trajmsm_long <- merge(obsdata_long, datapost, by = "id")
#'     AggFormula <- as.formula(paste("statins", "~", "time", "+", "class"))
#'     AggTrajData <- aggregate(AggFormula, data = trajmsm_long, FUN = mean)
#'     AggTrajData
#' trajmsm_wide = reshape(data = trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper"), timevar = "time", sep ="")
#'formula = paste0("y ~", paste0(treatment_var,collapse = "+"), "+",
#'                 paste0(unlist(covariates), collapse = "+"),"+",
#'                 paste0(baseline_var, collapse = "+"))
#'trajmsm_gform(formula = formula, identifier = "id",rep = 2, baseline = baseline_var, covariates = covariates,
#'                                     treatment = treatment_var, outcome = "y", total_followup = 6,time = "time",
#'                                     trajmodel = restraj$traj_model,ref = "1", obsdata = trajmsm_wide)
#' @author Awa Diop Denis Talbot

trajmsm_gform <- function(formula = formula, rep = 50,
                                  identifier,baseline,covariates,treatment,outcome, total_followup,time = time,
                                  trajmodel,ref, obsdata){


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
                                                 treatment = treatment_name, time = time,
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
                                               treatment = treatment_name, time = time,
                                               trajmodel = trajmodel)$post_class);
    obsdataG$gform_group <- relevel(as.factor(obsdataG$gform_group), ref = ref)
    #Estimation
    mod = summary(glm(Y ~ gform_group, family = quasibinomial, data = obsdataG));
    coefs = mod$coefficients[,1]
    se = apply(result.coef.boot,2,sd)
    pvalue <- 2*pnorm(-abs(coefs)/se)
    lo.ci = coefs - 1.96*se
    up.ci = coefs + 1.96*se
    results_msm_gform = rbind(coefs,se,pvalue,lo.ci, up.ci);
    rownames(results_msm_gform) = c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")
    return(t(results_msm_gform));
}
