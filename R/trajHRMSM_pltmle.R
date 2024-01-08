#' @title Pooled LTMLE
#' @description function to estimate parameters of a HRMSM-LCGA using Pooled LTMLE
#' @name trajHRMSM_pltmle
#' @param L time-varying covariates.
#' @param V baseline covariates.
#' @param A time-varying treatment.
#' @param id name of the id column variable.
#' @param Rep number of repetitions for the bootstrap.
#' @param s number of measuring times per interval.
#' @param K total length of follow-up.
#' @param time measuring times.
#' @param timevar name of the variable time.
#' @param obsdata observed data in wide format.
#' @noRd
#' @return \item{results_hrmsm_pltmle}{Results from the LCGA-HRMSM with pooled ltmle}
#' @author Awa Diop Denis Talbot


trajHRMSM_pltmle <- function(obsdata,
                                 degree_traj = c("linear","quadratic","cubic"),
                                 A,L,V,Y, s,K, timevar, idvar,
                                 J = 3, family = "poisson"){

  stopifnot(!is.null(list_obsdata));
  stopifnot(!is.null(s));
  stopifnot(!is.null(A));
  stopifnot(!is.null(idvar));

  size = 2^s
  nb_sub = K-s +1
  list_obsdataG <- list()
  list_obsdata = split_data(obsdata = obsdata, K = K, s = s, timevar = timevar, idvar = idvar)
  #Trajectory model and identification of the reference group
  dat_sub <- data.frame(do.call(rbind, list_obsdata))

  #Choice of degree for the polynomial form to build the trajectory groups
  if(degree_traj == "linear"){
    res_traj = buildtraj(Rdat = na.omit(dat_sub[,c(A,"time2","id2")]), J=J,formula = cbind(A,1-A) ~ time2, id="id2")
  }

  if(degree_traj == "quadratic"){
    res_traj = buildtraj(Rdat = na.omit(dat_sub[,c(A,"time2","id2")]), J=J,formula = cbind(A,1-A) ~ time2 + I(time2^2), id="id2")
  }

  if(degree_traj == "cubic"){
    res_traj = buildtraj(Rdat = na.omit(dat_sub[,c(A,"time2","id2")]), J=J,formula = time2 + I(time2^2) + I(time2^3), id="id2")
  }

  trajmodel = res_traj$model
  dclass <- data.frame(traj = factor(res_traj$dpost[,"class"]), id2 = res_traj$dpost[,"id2"])
  dat_final <- merge(dat_sub, dclass, by = "id2")
  mean_adh <- aggregate(as.formula(paste0(A, "~", "traj")), FUN = mean, data = dat_final)
  ord_adh<- order(-mean_adh[,2])
  ref <- ord_adh[length(ord_adh)]
  nb_regimes = 2^s #number of treatment regimes

  # Prediction of trajectory groups for each treatment regime
  class = predicTraj(t = s, trajmodel = trajmodel,
                     trt = A, time_name = "time2", id = "id2")$postclass

  Traj=t(sapply(1:nb_regimes,function(x)sapply(1:J,function(i) ifelse(class[x]==i,1,0))))
  Traj[,1]=1 #Intercept

  list_daTpool <- list()
  list_D <- list()
  df = list_obsdata[[1]];
  df1 = longtowide(df, v.names = c(A,L,Y,"time", "Interv","time2","id2"))
  cnames <- colnames(df1)
  for(i in 1:nb_sub){
    form = paste0(paste0(Y,".",nb_sub, "~"), paste0(A,".", 1:s,collapse = "+"), "+",
                  paste0(L,".", i:s,collapse = "+"),"+",
                  V, collapse = "+")
  #Create the data under all the different regime of treatment
  df = list_obsdata[[i]];
  df_l = longtowide(df, v.names = c(A,L,Y,"time", "Interv","time2","id2"))
  colnames(df_l)<-cnames
  res = sub_pltmle2(obsdata = df_l,Traj = Traj,formula = form,
                        Y=Y,A=A,L=L,V=V,s=s, time = 1:s,timevar = timevar)

  datTpool = data.frame(Y=res$counter.means)
  datTpool$TMLE.group = class
  datTpool$time <- i

  list_daTpool[i] <- list(datTpool)
  list_D[i] = list(res$D)

  }

  all_datTpool <- data.frame(do.call(rbind, list_daTpool))
  all_datTpool$TMLE.group <- relevel(factor(all_datTpool$TMLE.group), ref = ref)
  # Estimation
  mod = glm(Y ~ factor(TMLE.group) + factor(time), data =  all_datTpool, family = family);
  coefs = summary(mod)$coefficients[1:J,1];


  #Influence functions
  Db_list  <- list()
  Xall = t(model.matrix(mod))[1:J,]
  X = Xall[,1:size]
  B  = matrix(coefs, nrow = J);

  len = seq(size,ncol(Xall),size-1)
  # We have to select the corresponding model matrix for each time-interval (size*nb_sub)
  CQ_list <- list()
  for(t in 1:(nb_sub-1)){
    #x = pairs[,t]
    Db = matrix(0, nrow = nrow(list_D[[t]][[1]]), ncol = J);
    CQ = lapply(1:size,function(i)(as.matrix(X[,i]))%*%(t(exp(as.matrix(t(X[,i]))%*%B)))%*%t(as.matrix(((X[,i])))));
    CQ = Reduce('+',CQ);
    CQ_list[t] <- list(CQ)
    for(l in 1:size){
      Db = Db+as.matrix(list_D[[t]][[l]])%*%solve(CQ);
    }

    Db_list[t] <- list(Db)
    X  = Xall[,len[t]:len[t+1]]
  }

  #Last Window
  X = Xall[, len[nb_sub-1]:ncol(Xall)]
  CQ = lapply(1:size,function(i)(as.matrix(X[,i]))%*%(t(exp(as.matrix(t(X[,i]))%*%B)))%*%t(as.matrix(((X[,i])))));
  CQ = Reduce('+',CQ);
  CQ_list[3] <- list(CQ)
  Db = matrix(0, nrow = nrow(list_D[[nb_sub]][[1]]), ncol = J);
  for(l in 1:size){
    Db = Db+as.matrix(list_D[[nb_sub]][[l]])%*%solve(CQ);
  }

  Db_list[nb_sub] <- list(as.matrix(Db))

  #Computation

  #All pairs of 2 without repetition
  pairs = combn(nb_sub, 2)
  list_df = Db_list
  #sample size for each dataframe of influences functions
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

    temp_df1$id1 <- list_obsdata[[x[1]]][1:nrow(temp_df1),idvar]
    temp_df2$id2 <- list_obsdata[[x[2]]][1:nrow(temp_df2),idvar]

    vrc = sapply(1:J, function(j){
      res = mean(temp_df1[temp_df1$id1%in%temp_df2$id2,j]*temp_df2[temp_df2$id2%in%temp_df1$id1,j],na.rm = TRUE)-
        mean(temp_df1[temp_df1$id1%in%temp_df2$id2,j],na.rm = TRUE)*mean(temp_df2[temp_df2$id2%in%temp_df1$id1,j],na.rm = TRUE)
      return(res)})
    return(vrc)}
  )

  # minimum sample size per pairs of values
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
  #Results
  lo.ci = coefs - 1.96*se ;
  up.ci = coefs + 1.96*se
  results_hrmsm_pltmle = cbind(coefs, se, lo.ci, up.ci)
  colnames(results_hrmsm_pltmle) = c("estimate","std.error","lower CI", "Upper CI")
  return(list(res_trajHRMSM = results_hrmsm_pltmle[2:J,], res_traj = res.traj, mean_adh = mean_adh))
}


