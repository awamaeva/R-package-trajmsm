#' @title Estimating parameters of LCGA-HRMSM using g-formula
#' @description function to estimate parameters of a LCGA-HRMSM using g-formula
#'  and bootstrap to get standard errors.
#' @name trajhrmsm_gform
#' @param formula specification of the model for the outcome to be fitted for a binomial or gaussian distribution.
#' @param family specification of the error distribution and link function to be used in the model.
#' @param degree_traj to specify the polynomial degree for modelling the time-varying treatment.
#' @param identifier  name of the column for unique identifiant.
#' @param baseline name of baseline covariates.
#' @param covariates names of time-varying covariates in a wide format.
#' @param treatment name of time-varying treatment.
#' @param var_cov names of the time-varying covariates in a long format.
#' @param ntimes_interval length of a time-interval (s).
#' @param total_followup total length of follow-up.
#' @param censor name of the censoring variable.
#' @param number_traj number of trajectory groups.
#' @param rep number of repetition for the bootstrap.
#' @param obsdata data in a long format.
#' @param time name of the time variable.
#' @param time_values measuring times.
#' @return \item{trajMSM_gform }{g-formula}
#' @export
#' @author Awa Diop Denis Talbot
#' @examples
#' obsdata_long = gendata(n = 1000, format = "long", total_followup = 8,
#' timedep_outcome = TRUE,  seed = 945)
#' baseline_var <- c("age","sex")
#' years <- 2011:2018
#' variables <- c("hyper", "bmi")
#' covariates <- lapply(years, function(year) {
#' paste0(variables, year)})
#' treatment_var <- paste0("statins", 2011:2018)
#' var_cov <- c("statins","hyper", "bmi")
#' reshrmsm_gform = trajhrmsm_gform(degree_traj = "linear", rep=5 ,
#' treatment = treatment_var,covariates = covariates, baseline = baseline_var,
#' outcome = "y",var_cov = var_cov, ntimes_interval = 6, total_followup = 8,
#'  time = "time",time_values = years, identifier = "id",
#' number_traj = 3, family = "poisson", obsdata = obsdata_long)
#'reshrmsm_gform$results_hrmsm_gform



trajhrmsm_gform <- function(degree_traj = c("linear","quadratic","cubic"),
                            rep=50,treatment,covariates,baseline,outcome, ntimes_interval,
                            total_followup, time,time_values, identifier, var_cov,
                            number_traj = 3, family = "poisson",obsdata){

  nb_sub = total_followup - ntimes_interval + 1
  list_obsdataG <- list()
  list_obsdata = split_data(obsdata = obsdata, total_followup = total_followup,
                            ntimes_interval = ntimes_interval,
                            time = time,time_values = time_values, identifier = identifier)
  #Trajectory model and identification of the reference group
  dat_sub <- data.frame(do.call(rbind, list_obsdata))
  treatment_names <- sub("\\d+", "", treatment)
  treatment_name <- unique(treatment_names)[1]
  #Choice of degree for the polynomial form to buil the trajectory groups
  if(degree_traj == "linear"){

    restraj = build_traj(obsdata  = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]),
                          number_traj = number_traj,formula = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2")),
                          identifier = "identifier2")
  }

  if(degree_traj == "quadratic"){

    restraj = build_traj(obsdata = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]), number_traj = number_traj,
                          formula  = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2 + I(time2^2)")),
                          identifier = "identifier2")
  }

  if(degree_traj == "cubic"){

    restraj = build_traj(obsdata = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]), number_traj = number_traj,
                          formula  = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2 + I(time2^2) + I(time2^3)")),
                          identifier = "identifier2")
  }

  if(length(unique(restraj$data_post$class)) < number_traj){stop("number of trajectory groups identified is inferior to the target number.")}

  if(length(unique(restraj$data_post$class)) == number_traj){

  dclass <- data.frame(traj = factor(restraj$data_post[,"class"]), identifier2 = restraj$data_post[,"identifier2"])
  colnames(dclass)[1] <- "gform_group"

  dat_final <- merge(dat_sub, dclass, by = "identifier2")
  mean_adh <- aggregate(as.formula(paste0(treatment_name, "~", "gform_group")), FUN = mean, data = dat_final)
  ord_adh<- order(-mean_adh[,2])
  ref <- ord_adh[length(ord_adh)]

  bootf=function(df,x, ref){

    for(i in 1:nb_sub){
      #sampling in each time interval
      df = list_obsdata[[i]];
      df = df[df[,identifier]%in% x[[i]],]
      df_l = reshape(df, direction = "wide", idvar = identifier, v.names = var_cov, timevar = time, sep ="")

      outcome_up <- outcome[outcome %in% colnames(df_l)]
      treatment_up <- treatment[treatment %in% colnames(df_l)]
      cov_up <- lapply(covariates, function(x)x[x %in% colnames(df_l)])
      # Remove elements that are character(0)
      cov_up  <-  cov_up [sapply( cov_up , length) > 0]


      form = paste0(outcome_up, "~", paste0(treatment_up,collapse = "+"), "+",
                    paste0(unlist(cov_up), collapse = "+"),"+",
                    paste0(baseline, collapse = "+"))

      res = gformula(formula = form,outcome = outcome_up, treatment = treatment_up, covariates = cov_up,baseline = baseline,
                     ntimes_interval = ntimes_interval, obsdata = df_l)$counter_means
      colnames(res) <- "Y"

      #Predict trajectory groups based on each treatment regime
      obsdataG=data.frame(res)
      obsdataG$gform_group = factor(predict_traj(identifier = "identifier2", total_followup = ntimes_interval,
                                              treatment = treatment_name, time = "time2", time_values = 1:ntimes_interval,
                                              trajmodel =   restraj$traj_model)$post_class);
      obsdataG$Interv <- i
      list_obsdataG[i] <- list(obsdataG)
    }

    all_datG = data.frame(do.call(rbind, list_obsdataG))
    all_datG$gform_group <- relevel(factor(all_datG$gform_group), ref = ref)
    all_datG$Interv <- factor(all_datG$Interv)
    # Estimation
    mod = summary(glm(Y~factor(gform_group) + factor(Interv), family = family,data= all_datG));
    return(coef(mod)[,1]);
  }

  #Bootstrap
  result_coef_boot <- matrix(numeric(0), nrow = rep, ncol = number_traj)

  for(b in 1:rep){
    list_indices_temp <- list()
    list_indices <- list()

    for(t in 1:nb_sub){
      list_indices_temp[t] <- list(sample(unique(list_obsdata[[t]][,identifier]), length(unique(list_obsdata[[t]][,identifier])) ,replace = T))
    }
 # to ensure to sample same individuals
    for(t1 in 1:nb_sub){
      list_indices[t1] = list(list_obsdata[[t1]][,identifier][list_obsdata[[t1]][,identifier] %in% Reduce(intersect,  list_indices_temp)])
      }


    result_coef_boot[b,] <- bootf(df,x = list_indices, ref= ref)[1:number_traj]
    cat('Replication', b, 'of', rep,'\n')
  }

   for(i in 1:nb_sub){
    #sampling in each time interval
    df = list_obsdata[[i]];
    df_l = reshape(df, direction = "wide", idvar = identifier, v.names = var_cov, timevar = time, sep ="")

    outcome_up <- outcome[outcome %in% colnames(df_l)]
    treatment_up <- treatment[treatment %in% colnames(df_l)]
    cov_up <- lapply(covariates, function(x)x[x %in% colnames(df_l)])
    # Remove elements that are character(0)
    cov_up  <-  cov_up [sapply( cov_up , length) > 0]


    form = paste0(outcome_up, "~", paste0(treatment_up,collapse = "+"), "+",
                  paste0(unlist(cov_up), collapse = "+"),"+",
                  paste0(baseline, collapse = "+"))

    res = gformula(formula = form,outcome = outcome_up, treatment = treatment_up, covariates = cov_up,baseline = baseline,
                   ntimes_interval = ntimes_interval, obsdata = df_l)$counter_means
    colnames(res) <- "Y"

    #Predict trajectory groups based on each treatment regime
    obsdataG=data.frame(res)
    obsdataG$gform_group = factor(predict_traj(identifier = "identifier2", total_followup = ntimes_interval,
                                           treatment = treatment_name, time = "time2", time_values = 1:ntimes_interval,
                                           trajmodel = restraj$traj_model)$post_class);
    obsdataG$Interv <- i
    list_obsdataG[i] <- list(obsdataG)
  }

  all_datG = data.frame(do.call(rbind, list_obsdataG))
  all_datG$gform_group <- relevel(factor(all_datG$gform_group), ref = ref)
  all_datG$Interv <- factor(all_datG$Interv)

  # Estimation
  mod_glm = summary(glm(Y~factor(gform_group) + factor(Interv), family = family,data= all_datG));
  coefs = mod_glm$coef[1:number_traj,1]
  se = apply(result_coef_boot,2,sd)
  pvalue <- 2*pnorm(-abs(coefs)/se)
  lo.ci = coefs - 1.96*se
  up.ci = coefs + 1.96*se
  results = cbind(coefs ,se,pvalue, lo.ci, up.ci);
  colnames(results) =  c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI");
  results_hrmsm_gform = results
  return(list(  results_hrmsm_gform = results_hrmsm_gform, result_coef_boot = result_coef_boot, restraj = restraj, mean_adh = mean_adh));
  }

}