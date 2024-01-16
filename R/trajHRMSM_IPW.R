#' @title Marginal Structural Model and Latent Class of Growth Analysis estimated with IPW
#' @description Combine Marginal Structural Model and Latent Class of Growth Analysis
#' @name trajhrmsm_ipw
#' @param formulaY specification of the model for the outcome to be fitted for a binomial or gaussian distribution.
#' @param formula_traj specification of the model to buil the trajectory groups.
#' @param family specification of the error distribution and link function to be used in the model.
#' @param idvar  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param s number of measuring times per interval.
#' @param K total length of follow-up.
#' @param censor name of the censoring variable.
#' @param J number of trajectory groups.
#' @param weights a vector of estimated weights. If NULL, the weights are computed by the function.
#' @param v.names
#' @author Awa Diop, Denis Talbot
#' @return \item{IPW}{Stabilized and unstabilized inverse of probabilities}
#' @export trajHRMSM_IPW
#' @importFrom sandwich
#' @import flexmix
#' @examples
#' obsdata <- genobsdat(n=5000)
#' V <- "V"
#' L <- "L"
#' A <- "A"
#' Y <- "Y"
#' v.names = c("A","L", "Y","time")
#' trajHRMSM_IPW(formulaY = as.formula("Y ~ factor(traj) + factor(Interv)"), numerator = "stabilized",
#'           idvar = "id", V = V,L = L,A = A,Y=Y,
#'           C = FALSE,s = 3, K = 5, J=3, family = "poisson",
#'           obsdata = obsdata, v.names = v.names)


trajHRMSM_IPW <- function(formulaY = as.formula("Y ~ factor(traj) + factor(Interv)"),
                          degree_traj = c("linear","quadratic","cubic"),
                          numerator = c("stabilized", "unstabilized"),
                          idvar, V,L,A,Y, C = FALSE, s,K,time,timevar,
                          family = "poisson",censor = censor,
                          J, obsdata, weights = NULL, v.names, level = 0.999){

  if(is.null(weights)){
    obsdata_wide <- longtowide(obsdata, v.names = v.names)
    wts = IPW(numerator = numerator, id = idvar , V = V,L = paste0(L,"."),A = paste0(A,"."), C = C,
              censor = censor, K = K,time = time, obsdata = obsdata_wide);
    colnames(wts)<- paste0("w",1:K)
    }

  if(!is.null(weights)){
    wts <- weights
  }

  #Level of truncation
  trunc <- function(x,level){
    x = ifelse(x>quantile(x,level, na.rm=TRUE),quantile(x,level,na.rm=TRUE), x)
    return(x)
  }
  obsdata$IPW <- trunc(as.vector(wts),level=level)
  dat_sub = data.frame(do.call(rbind, split_data(obsdata = obsdata, K = K, s = s, timevar = timevar, idvar = idvar)))

  #Choice of degree for the polynomial form to buil the trajectory groups
  if(degree_traj == "linear"){
  res_traj = build_traj(Rdat = na.omit(dat_sub[,c(A,"time2","id2")]), J=J,formula = cbind(A,1-A) ~ time2, id="id2")
  }

  if(degree_traj == "quadratic"){
    res_traj = build_traj(Rdat = na.omit(dat_sub[,c(A,"time2","id2")]), J=J,formula = cbind(A,1-A) ~ time2 + I(time2^2), id="id2")
  }

  if(degree_traj == "cubic"){
    res_traj = build_traj(Rdat = na.omit(dat_sub[,c(A,"time2","id2")]), J=J,formula = time2 + I(time2^2) + I(time2^3), id="id2")
  }

     dclass <- data.frame(traj = factor(res_traj$dpost[,"class"]), id2 = res_traj$dpost[,"id2"])
    dat_final <- merge(dat_sub, dclass, by = "id2")

  mean_adh <- aggregate(as.formula(paste0(A, "~", "traj")), FUN = mean, data = dat_final)
  ord_adh<- order(-mean_adh[,2])
  ref <- ord_adh[length(ord_adh)]
  dat_final$traj <- relevel(factor(dat_final$traj), ref = ref)

    mod_glm = glm(formula = formulaY, weights = IPW,family = family,data=dat_final);
    coefs <- summary(mod_glm)$coefficients[1:J,1];
     se <- summary(mod_glm)$coefficients[1:J,2];
    #se = sqrt(diag(vcovCL(mod_glm, cluster = ~ "id2" )))[1:J];
    IClo = coefs- 1.96*se ;
    ICup = coefs + 1.96*se;
    res_trajHRMSM = cbind(coefs,se,IClo, ICup);
    colnames(res_trajHRMSM) = c("estimate","std.error","lower CI", "Upper CI");
    res_trajHRMSM = res_trajHRMSM[2:J,]
return(list(res_trajHRMSM = res_trajHRMSM, res_traj = res_traj, mean_adh = mean_adh))
}
