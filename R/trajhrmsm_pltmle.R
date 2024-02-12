#' @title History Restricted MSM and Latent Class of Growth Analysis estimated with a Pooled LTMLE.
#' @description Estimate parameters of LCGA-HRMSM using a Pooled LTMLE.
#' @name trajhrmsm_pltmle
#' @param family Specification of the error distribution and link function to be used in the model.
#' @param degree_traj To specify the polynomial degree for modelling the time-varying treatment.
#' @param identifier  Name of the column for unique identifiant.
#' @param baseline Names of baseline covariates.
#' @param covariates Names of time-varying covariates (should be a list).
#' @param treatment Name of time-varying treatment.
#' @param outcome Name of the outcome variable.
#' @param var_cov Names of the time-varying covariates.
#' @param ntimes_interval Length of a time-interval (s).
#' @param total_followup Total length of follow-up.
#' @param number_traj Number of trajectory groups.
#' @param treshold For weight truncation.
#' @param obsdata Data in a long format.
#' @param time Name of the time variable.
#' @param time_values Measuring times.
#' @importFrom stats na.omit rbinom plogis qlogis  reshape glm
#' binomial coef as.formula ave aggregate relevel pnorm sd quantile model.matrix
#' quasibinomial var
#' @importFrom utils combn
#' @export
#' @return A list containing the following components:
#'   \describe{
#'   \item{results_hrmsm_pltmle}{ Matrix of estimates for LCGA-HRMSM, obtained using the pooled ltlmle method.}
#'   \item{restraj}{ Fitted trajectory model.}
#'   \item{mean_adh}{ Matrix of the mean adherence per trajectory group.}
#'   }
#' @author Awa Diop Denis Talbot
#' @examples
#' \donttest{
#' obsdata_long = gendata(n = 1000, format = "long",
#' total_followup = 8, timedep_outcome = TRUE,  seed = 945)
#' baseline_var <- c("age","sex")
#' years <- 2011:2018
#' variables <- c("hyper", "bmi")
#' covariates <- lapply(years, function(year) {
#'   paste0(variables, year)})
#' treatment_var <- paste0("statins", 2011:2018)
#' var_cov <- c("statins","hyper", "bmi","y")
#' respltmle = trajhrmsm_pltmle(degree_traj = "linear", treatment = treatment_var,
#' covariates = covariates, baseline = baseline_var,
#' outcome = paste0("y", 2016:2018),var_cov = var_cov, ntimes_interval = 6,
#' total_followup = 8, time = "time",time_values = years, identifier = "id",
#' number_traj = 3, family = "poisson", obsdata = obsdata_long)
#' respltmle$results_hrmsm_pltmle
#' }



trajhrmsm_pltmle <-  function(degree_traj = c("linear","quadratic","cubic"),
                              treatment,covariates,baseline,outcome, ntimes_interval,
                              total_followup, time, time_values,identifier, var_cov,
                              number_traj = 3, family = "poisson",obsdata, treshold = 0.99){

  nb_sub = total_followup - ntimes_interval + 1
  nregimes = 2^ntimes_interval
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

  dclass <- data.frame(ptmle_group = factor(restraj$data_post[,"class"]), identifier2 = restraj$data_post[,"identifier2"])
  dat_final <- merge(dat_sub, dclass, by = "identifier2")
  mean_adh <- aggregate(as.formula(paste0(treatment_name, "~", "ptmle_group")), FUN = mean, data = dat_final)
  ord_adh<- order(-mean_adh[,2])
  ref <- as.character(ord_adh[length(ord_adh)])
  nregimes = 2^ntimes_interval #number of treatment regimes

  # Prediction of trajectory groups for each treatment regime
  treatment_names <- sub("\\d+", "", treatment)
  treatment_name <- unique(treatment_names)[1]
  class = factor(predict_traj(identifier = "identifier2", total_followup = ntimes_interval,
                              treatment = treatment_name, time = "time2", time_values = 1:ntimes_interval,
                              trajmodel = restraj$traj_model)$post_class);

  if(length(unique(class)) < number_traj){stop("number of trajectory groups identified is inferior to the target number.")}

  if(length(unique(class)) == number_traj){

    traj_indic=t(sapply(1:nregimes,function(x)sapply(1:number_traj,function(i) ifelse(class[x]==i,1,0))))
    traj_indic[,1]=1 #Intercept

    list_obsdata_pool <- list()
    list_D <- list()
    df <- data.frame()

    for(i in 1:nb_sub){
      #Create the data under all the different regime of treatment
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

      res_pltmle = pltmle(formula = form, outcome = outcome_up,treatment = treatment_up,
                          covariates = cov_up, baseline = baseline, ntimes_interval = ntimes_interval, number_traj = number_traj,
                          time =  time,identifier = identifier,obsdata = df_l,traj=traj_indic, treshold = treshold);

      obsdata_pool= data.frame(Y=res_pltmle$counter_means);
      D=res_pltmle$D; #Influence functions
      obsdata_pool$tmle_group = class
      obsdata_pool$Interv <- i

      list_obsdata_pool[i] <- list(obsdata_pool)
      list_D[i] = list(D)

    }

    all_obsdata_pool <- data.frame(do.call(rbind, list_obsdata_pool))
    all_obsdata_pool$tmle_group <- relevel(factor(all_obsdata_pool$tmle_group), ref = ref)
    # Estimation
    mod = glm(Y ~ factor(tmle_group) + factor(Interv), data =  all_obsdata_pool, family = family);
    coefs = summary(mod)$coefficients[1:number_traj,1];


    #Influence functions
    Db_list  <- list()
    Xall = t(model.matrix(mod))[1:number_traj,]
    X = Xall[,1:nregimes]
    B  = matrix(coefs, nrow = number_traj);

    len = seq(nregimes,ncol(Xall),nregimes-1)
    # We have to select the corresponding model matrix for each time-interval (nregimes*nb_sub)
    CQ_list <- list()
    for(t in 1:(nb_sub-1)){
      #x = pairs[,t]
      Db = matrix(0, nrow = nrow(list_D[[t]][[1]]), ncol = number_traj);
      CQ = lapply(1:nregimes,function(i)(as.matrix(X[,i]))%*%(t(exp(as.matrix(t(X[,i]))%*%B)))%*%t(as.matrix(((X[,i])))));
      CQ = Reduce('+',CQ);
      CQ_list[t] <- list(CQ)
      for(l in 1:nregimes){
        Db = Db+as.matrix(list_D[[t]][[l]])%*%solve(CQ);
      }

      Db_list[t] <- list(Db)
      X  = Xall[,len[t]:len[t+1]]
    }

    #Last Window
    X = Xall[, len[nb_sub-1]:ncol(Xall)]
    CQ = lapply(1:nregimes,function(i)(as.matrix(X[,i]))%*%(t(exp(as.matrix(t(X[,i]))%*%B)))%*%t(as.matrix(((X[,i])))));
    CQ = Reduce('+',CQ);
    CQ_list[3] <- list(CQ)
    Db = matrix(0, nrow = nrow(list_D[[nb_sub]][[1]]), ncol = number_traj);
    for(l in 1:nregimes){
      Db = Db+as.matrix(list_D[[nb_sub]][[l]])%*%solve(CQ);
    }

    Db_list[nb_sub] <- list(as.matrix(Db))

    #Computation

    #All pairs of 2 without repetition
    pairs = combn(nb_sub, 2)
    list_df = Db_list
    #sample nregimes for each dataframe of influences functions
    n_df <- lapply(list_obsdata, function(x) nrow(na.omit(data.frame(x))))

    #variances

    var <- lapply(1:nb_sub,function(i){
      vr = diag(var(data.frame(list_df[[i]]),na.rm = TRUE))
      return(vr)}
    )

    #covariance
    cov <- lapply(1:ncol(pairs),function(i){
      x = pairs[,i]
      temp_df1 = data.frame(list_df[[x[1]]])
      temp_df2 = data.frame(list_df[[x[2]]])

      temp_df1$id1 <- list_obsdata[[x[1]]][1:nrow(temp_df1),identifier]
      temp_df2$id2 <- list_obsdata[[x[2]]][1:nrow(temp_df2),identifier]

      vrc = sapply(1:number_traj, function(j){
        res = mean(temp_df1[temp_df1$id1%in%temp_df2$id2,j]*temp_df2[temp_df2$id2%in%temp_df1$id1,j],na.rm = TRUE)-
          mean(temp_df1[temp_df1$id1%in%temp_df2$id2,j],na.rm = TRUE)*mean(temp_df2[temp_df2$id2%in%temp_df1$id1,j],na.rm = TRUE)
        return(res)})
      return(vrc)}
    )

    # minimum sample nregimes per pairs of values
    min_ndf <- lapply(1:ncol(pairs),function(i){
      x = pairs[,i]
      min_ndf = min(n_df[[x[1]]], n_df[[x[2]]])
      return(min_ndf)
    })

    all_cov <- lapply(1:length(min_ndf),function(i){

      cov_temp = 2*min_ndf[[i]]*cov[[i]]
    })

    all_var = lapply(1:nb_sub,function(i){
      var_temp= n_df[[i]]*as.numeric(var[[i]])
      return(var_temp)}
    )

    se = sqrt((Reduce('+',all_cov) + Reduce('+',all_var))/(Reduce('+',n_df)**2))
    pvalue <- 2*pnorm(-abs(coefs)/se)
    #Results
    lo.ci = coefs - 1.96*se ;
    up.ci = coefs + 1.96*se
    results_hrmsm_pltmle = cbind(coefs, se, pvalue, lo.ci, up.ci)
    colnames(results_hrmsm_pltmle) = c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI");
    return(list(results_hrmsm_pltmle= results_hrmsm_pltmle, restraj = restraj, mean_adh = mean_adh))
  }

}
