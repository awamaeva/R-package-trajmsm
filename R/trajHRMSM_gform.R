#' @title g-formula without SuperLearner
#' @description function to estimate parameters of a HRMSM-LCGA using g-formula
#'  and bootstrap to get standard errors.
#' @name trajHRMSM_gform
#' @param id  name of the column for unique identifiant.
#' @param L time-varying covariates.
#' @param V baseline covariates.
#' @param A time-varying treatment.
#' @param Rep number of repetitions for the bootstrap.
#' @param s number of measuring times per interval.
#' @param K total length of follow-up.
#' @param time measuring times.
#' @param time_name name of the variable time.
#' @param list_obsdata the data for each time interval.
#' @param obsdata observed data.
#' @noRd
#' @return \item{trajMSM_gform }{g-formula}
#' @author Awa Diop Denis Talbot

trajHRMSM_gform <- function(obsdata, degree_traj = c("linear","quadratic","cubic"),
                            Rep=50,A,L,V,Y, s,K, timevar = timevar, idvar,
                            J = 3, family = "poisson"){
  stopifnot(!is.null(list_obsdata));
  stopifnot(!is.null(s));
  stopifnot(!is.null(K));
  stopifnot(!is.null(Rep));
  stopifnot(!is.null(A));
  stopifnot(!is.null(L));
  stopifnot(!is.null(Y));
  stopifnot(!is.null(idvar));


  nb_sub = K-s+1
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
  bootf=function(df,x, ref){
    df = list_obsdata[[1]];
    df = df[df[,idvar]%in% x[[1]],]
    df1 = longtowide(df, v.names = c(A,L,Y,timevar, "Interv","time2","id2"))
    cnames <- colnames(df1)

    for(i in 1:nb_sub){
      #sampling in each time interval
      df = list_obsdata[[i]];
      df = df[df[,idvar]%in%x[[i]],]
      df_l = longtowide(df, v.names = c(A,L,Y,timevar, "Interv","time2","id2"))
      colnames(df_l)<-cnames
      form = paste0(paste0(Y,".",nb_sub, "~"), paste0(A,".", 1:s,collapse = "+"), "+",
                    paste0(L,".", 1:s,collapse = "+"),"+",
                    V, collapse = "+")
      res = sub_gform(obsdata = df_l, formula = form, Y=paste0(Y,"."),
                      A=paste0(A,"."),L=paste0(L,"."),V=V, s=s)$counter_means
      colnames(res) <- "Y"

      #Predict trajectory groups based on each treatment regime
      datG=data.frame(res)
      datG$Gform.group = predicTraj(t = s, trajmodel = trajmodel,
                                    trt = A,
                                    time_name = "time2", id = "id2")$postclass
      datG$time <- i
      list_obsdataG[i] <- list(datG)
    }

    all_datG = data.frame(do.call(rbind, list_obsdataG))
    all_datG$Gform.group <- relevel(factor(all_datG$Gform.group), ref = ref)
    all_datG$time <- factor(all_datG$time)
    # Estimation
    mod = summary(glm(Y~factor(Gform.group) + factor(time), family = family,data= all_datG));
    return(coef(mod)[,1]);
  }

  #Bootstrap
  result.coef.boot <- matrix(numeric(0), nrow = Rep, ncol = J)

  for(b in 1:Rep){
    list_indices_temp <- list()
    list_indices <- list()

    for(t in 1:nb_sub){
      list_indices_temp[t] <- list(sample(list_obsdata[[t]][,idvar], nrow(list_obsdata[[t]]) ,replace = T))
    }
 # to ensure to sample same individuals
    for(t1 in 1:nb_sub){
      list_indices[t1] = list(list_obsdata[[t1]][,idvar][list_obsdata[[t1]][,idvar] %in% Reduce(intersect,  list_indices_temp)])
    }


    result.coef.boot[b,] <- bootf(df,x = list_indices, ref= ref)[1:J]
  }

  #results
  df = list_obsdata[[1]];
  df = df[df[,idvar]%in% list_indices[[1]],]
  df1 = longtowide(df, v.names = c(A,L,Y,timevar, "Interv","time2","id2"))
  cnames <- colnames(df1)
  for(i in 1:nb_sub){
    #sampling in each time interval
    df = list_obsdata[[i]];
    df_l = longtowide(df, v.names = c(A,L,Y,timevar, "Interv","time2","id2"))
    colnames(df_l) <- cnames
    form = paste0(paste0(Y,".",nb_sub, "~"), paste0(A,".", 1:s,collapse = "+"), "+",
                  paste0(L,".", 1:s,collapse = "+"),"+",
                  V, collapse = "+")
    res = sub_gform(obsdata = df_l, formula = form, Y=paste0(Y,"."),
                    A=paste0(A,"."),L=paste0(L,"."),V=V, s=s)$counter_means
    colnames(res) <- "Y"

    #Predict trajectory groups based on each treatment regime
    datG=data.frame(res)
    datG$Gform.group = predicTraj(t = s, trajmodel = trajmodel,
                                  trt = A,
                                  time_name = "time2", id = "id2")$postclass
    datG$time <- i
    list_obsdataG[i] <- list(datG)
  }

  all_datG = data.frame(do.call(rbind, list_obsdataG))
  all_datG$Gform.group <- relevel(factor(all_datG$Gform.group), ref = ref)
  all_datG$time <- factor(all_datG$time)

  # Estimation
  mod_glm = summary(glm(Y~factor(Gform.group) + factor(time), family = family,data= all_datG));
  coefs = mod_glm$coef[1:J,1]
  se = apply(result.coef.boot,2,sd)
  lo.ci = coefs - 1.96*se
  up.ci = coefs + 1.96*se
  results = cbind(coefs ,se, lo.ci, up.ci);
  colnames(results) = c("estimate","std.error","lower CI", "Upper CI")
  results_hrmsm_gform = results[2:J,]
  return(list(res_trajHRMSM = results_hrmsm_gform, result.coef.boot = result.coef.boot, res_traj = res.traj, mean_adh = mean_adh));
}
