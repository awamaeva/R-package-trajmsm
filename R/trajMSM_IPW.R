#' @title Marginal Structural Model and Latent Class of Growth Analysis estimated with IPW
#' @description Combine Marginal Structural Model and Latent Class of Growth Analysis
#' @name trajMSM_IPW
#' @param formula1 specification of the model for the outcome to be fitted for a binomial or gaussian distribution.
#' @param formula2 specification of the model for the outcome to be fitted for a survival outcome.
#' @param family specification of the error distribution and link function to be used in the model.
#' @param id  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param K number of measuring times.
#' @param time measuring times.
#' @param censor name of the censoring variable.
#' @param J number of trajectory groups.
#' @param weights a vector of estimated weights. If NULL, the weights are computed by the function.
#' @author Awa Diop, Denis Talbot
#' @return \item{IPW}{Stabilized and unstabilized inverse of probabilities}
#' @export trajMSM_IPW
#' @importFrom sandwich
#' @importFrom survival coxph
#' @import flexmix
#' @examples
#' obsdata_msm <- gen_obsdata_msm(n=5000,t=10, seed = 345)
#' V <- "V"
#' L <- "L"
#' time <- 1:10
#' K = 10
#' A <- "A"
#' Y <- "Y"
#'Adat_long <- widetolong(obsdata = obsdata_msm[, paste0(A,time)], varying = 1:10)
#'res.traj = build_traj(Rdat = Adat_long, J=3,formula = cbind(A,1-A) ~ time + I(time^2), id= "id")
#'table(res.traj$dpost$class)
#'obsdata_msm2 = merge(obsdata_msm,res.traj$dpost, by = "id")
#'obsdata_msm2[,"class"] <- relevel(as.factor(obsdata_msm2[,"class"]), ref = 3)
#'trajMSM_IPW(formula1 = as.formula("Y ~ class"), numerator = "stabilized",
#'            id = "id", V = V,L = L,A = A,Y=Y,J=3,
#'            C = FALSE,K = K, time = time, family = "binomial",
#'            obsdata = obsdata_msm2)


trajMSM_IPW <- function(formula1 = as.formula("Y ~ class"),
                        formula2 = as.formula("Surv(Y,event) ~ class"),
                        numerator = c("stabilized", "unstabilized"),
                        id, V,L,A,Y, C = FALSE,censor = NULL, K, time,
                        family = c("binomial","gaussian","survival"), J, obsdata, weights = NULL){

if(is.null(weights)){
wts = IPW(numerator = numerator, id = id , V = V,L = L,A = A, C = C,
              censor = censor, K = K,time = time,obsdata = obsdata)[,K];
obsdata$wts <- wts ;}

  if(!is.null(weights)){
    obsdata$wts <- weights
  }
if(family == "gaussian"){

  mod_glm <- glm(formula = formula1,
                                 weights = wts,family = gaussian,data = obsdata);
  coefs <- summary(mod_glm)$coefficients[,1];
  se <- sqrt(diag(vcovCL(mod_glm, cluster = ~id)))[1:J];
  pvalue <- summary(mod_glm)$coefficients[,3];
  IClo = coefs- 1.96*se ;
  ICup = coefs + 1.96*se;

  res_trajMSM = cbind(coefs,se,pvalue,IClo, ICup);
  colnames(res_trajMSM) = c("estimate","std.error","pvalue","lower CI", "Upper CI");
}

  if(family == "binomial"){
      mod_glm = glm(formula = formula1, weights = wts,family = binomial,data=obsdata);
      coefs <- summary(mod_glm)$coefficients[,1];
      se <- sqrt(diag(vcovCL(mod_glm, cluster = ~ id)))[1:J];
      pvalue <- summary(mod_glm)$coefficients[,3];
      IClo = coefs- 1.96*se ;
      ICup = coefs + 1.96*se;

      res_trajMSM = cbind(coefs,se,pvalue,IClo, ICup);
      colnames(res_trajMSM) = c("estimate","std.error","pvalue","lower CI", "Upper CI");
  }

  if(family == "survival"){
        mod_cox = coxph(formula = formula2, cluster=eval(id),weights = wts,data = obsdata)
        coefs <-  mod_cox$coefficients[,1];
        se <-     mod_cox$coefficients[,2];
        pvalue <- mod_cox$coefficients[,3];
        IClo = coefs- 1.96*se ;
        ICup = coefs + 1.96*se;
        res_trajMSM = cbind(coefs,se,pvalue,IClo, ICup);
        colnames(res_trajMSM) = c("estimate","std.error","pvalue","lower CI", "Upper CI");
  }

return(res_trajMSM)
}

