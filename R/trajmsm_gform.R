#' @title Parametric g-formula
#' @description function to estimate parameters of a MSM-LCGA using g-formula
#'  and bootstrap to get standard errors.
#' @name trajmsm_gform
#' @param formula specification of the model for the outcome to be fitted.
#' @param id name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param K number of measuring times.
#' @param time measuring times.
#' @param time_name name of the variable time.
#' @param rep number of repetitions for the bootstrap.
#' @param trajmodel trajectory model built with the observed treatment.
#' @param ref the reference group.
#' @param obsdata observed data in wide format.
#' @return \item{results_msm_gform}{Estimates of a LCGA-MSM with g-formula.}
#' @export trajMSM_gform
#' @examples
#' obsdata_msm <- gen_obsdata_msm(n=2500)
#' V <- "V"
#' L <- "L"
#' time <- 1:3
#' A <- "A"
#' Y <- "Y"
#' Adat_long <- widetolong(obsdata = obsdata_msm[, c("A1","A2","A3")], varying = 1:3)
#' res.traj = buildtraj(Rdat = Adat_long, J=3,formula = cbind(A,1-A) ~ time, id= "id")
#' obsdata_msm$class = res.traj$dpost$class
#'trajMSM_gform(formula = paste0("Y~", paste0("A", 1:K,collapse = "+"), "+",
#'                               paste0("L", 1:K,collapse = "+"),"+",
#'                               V, collapse = "+"),rep = 50,
#'              id = id, V = V,L = L,A = A,Y=Y,K=K,time = time,timevar = "time",J=3,
#'              trajmodel = res.traj$model, ref = "2", obsdata = obsdata)
#' @author Awa Diop Denis Talbot

trajMSM_gform <- function(formula = paste0("Y~", paste0("A", 1:K,collapse = "+"), "+",
                                          paste0("L", 1:K,collapse = "+"),"+",
                                          V, collapse = "+"), rep = 50,
                                  id,V,L,A,Y, K, time, timevar,
                          J = 2,trajmodel,obsdata, ref){


  stopifnot(!is.null(id));
  stopifnot(!is.null(V));
  stopifnot(!is.null(L));
  stopifnot(!is.null(A));
  stopifnot(!is.null(Y));
  stopifnot(!is.null(K));
  stopifnot(!is.null(J));
  stopifnot(!is.null(rep));
  stopifnot(!is.null(obsdata));
  stopifnot(!is.null(trajmodel));
  stopifnot(!is.null(time));



    bootf=function(df,x=indices){
      #Echantillons bootstrap
      df=obsdata[x,];
      res = sub_gform(formula = formula, Y = Y, A = A,L = L,V = V,
                      s = K, time = time,obsdata = df)$counter_means
      colnames(res) <- "Y";

      #Counterfactual means + trajectory groups
      obsdataG = data.frame(res);
      obsdataG$Gform_group = factor(predicTraj(t = K, trajmodel = trajmodel, trt = A,
                                        time_name = timevar, id = id)$postclass);
      obsdataG$Gform_group <- relevel(as.factor(obsdataG$Gform_group), ref = ref)
      #Estimation
      mod = summary(glm(Y ~ Gform_group, family = quasibinomial, data = obsdataG));
      return(coef(mod)[,1]);
    }

    #Bootstrap
    #result.coef.boot <- matrix(numeric(0), nrow = rep, ncol = J);
    list_res <- list()
    for(b in 1:rep){
      indices <- sample(1:nrow(obsdata), nrow(obsdata) ,replace = T)
      list_res[b] <- list(bootf(df,x = indices))
      cat('Replication', b, 'of', rep,'\n')
    }

    result.coef.boot <- do.call(rbind,list_res)
    #results
    res = sub_gform(formula = formula, Y = Y, A = A,L = L,V = V,
                    s = K, time = time,obsdata = obsdata)$counter_means
    colnames(res) <- "Y";

    #Counterfactual means + trajectory groups
    obsdataG = data.frame(res);
    obsdataG$Gform_group = factor(predicTraj(t = K, trajmodel = trajmodel, trt = A,
                                             time_name = timevar, id = id)$postclass);
    obsdataG$Gform_group <- relevel(as.factor(obsdataG$Gform_group), ref = ref)
    #Estimation
    mod = summary(glm(Y ~ Gform_group, family = quasibinomial, data = obsdataG));
    coefs.mean = mod$coefficients[,1]
    se=apply(result.coef.boot,2,sd)
    pvalue <- 2*pnorm(-abs(coefs.mean)/se)
    lo.ci = coefs.mean - 1.96*se
    up.ci = coefs.mean + 1.96*se
    results_msm_gform = rbind(coefs.mean,se,pvalue,lo.ci, up.ci);
    rownames(results_msm_gform) = c("estimate","std.error","pvalue","lower CI", "Upper CI")
    return(t(results_msm_gform));
}
