#' @title Pooled LTMLE
#' @description function to estimate parameters of a MSM-LCGA using pooled LTMLE
#'  with influence functions to estimate standard errors.
#' @name trajmsm_pltmle
#' @param formula specification of the model for the outcome to be fitted.
#' @param id  name of the column for unique identifiant.
#' @param L covariates.
#' @param A time-varying treatment.
#' @param time measuring times.
#' @param timevar name of the time variable.
#' @param K number of measuring times per interval.
#' @param J an integer to choose the number of trajectory groups.
#' @param trajmodel trajectory model built with the observed treatment.
#' @param obsdata observed data in wide format.
#' @param ref the reference group.
#' @export trajmsm_pltmle
#' @examples
#' obsdata_msm <- genobsdata_msm(n=2500)
#' V <- "V"
#' L <- "L"
#' K <- 10
#' time <- 1:10
#' A <- "A"
#' Y <- "Y"
#' Adat_long <- widetolong(obsdata = obsdata_msm[, c("A1","A2","A3")], varying = 1:3)
#' res.traj = buildtraj(Rdat = Adat_long, J=3,formula = cbind(A,1-A) ~ time, id= "id")
#' obsdata_msm$class = res.traj$dpost$class
#'trajMSM_pltmle(formula = paste0("Y~", paste0("A", 1:K,collapse = "+"), "+",
#'                               paste0("L", 1:K,collapse = "+"),"+",
#'                               V, collapse = "+"),
#'              id = id, V = V,L = L,A = A,Y=Y,K = K,time = time,time_name = "time",J=3,
#'             trajmodel = res.traj$model, obsdata = obsdata_msm)
#' @return \item{results_msm_pooledltmle}{Estimates of a LCGA-MSM with pooled LTMLE.}
#' @author Awa Diop, Denis Talbot


trajMSM_pltmle <- function(formula = paste0("Y~", paste0("A", 1:K,collapse = "+"), "+",
                                                  paste0("L", 1:K,collapse = "+"),"+",
                                                  paste0("V", 1:K,collapse = "+"), collapse = "+"),
                                 id,V,L,A,Y, K, time, timevar,J,
                                 trajmodel,ref,obsdata){

  stopifnot(!is.null(id));
  stopifnot(!is.null(V));
  stopifnot(!is.null(L));
  stopifnot(!is.null(A));
  stopifnot(!is.null(Y));
  stopifnot(!is.null(K));
  stopifnot(!is.null(J));
  stopifnot(!is.null(ref));
  stopifnot(!is.null(obsdata));
  stopifnot(!is.null(trajmodel));
  stopifnot(!is.null(time));



     nb_regimes = 2^K;  #number of treatment regimes
     expit = plogis;

     class = predicTraj(t = K, trajmodel = trajmodel,
                        trt = A, time_name = timevar, id = id)$postclass;
     Traj=t(sapply(1:nb_regimes,function(x)sapply(1:J,function(i) ifelse(class[x]==i,1,0))))
     Traj[,1]=1 #Intercept
    #Create the obsdataa under all the different regime of treatment
    res = sub_pltmle(formula = formula, id, V = V, L = L, A = A, Y = Y,
                           s = K,J=J,
                          time = time,
                          timevar= timevar, Traj = Traj,obsdata = obsdata);

    obsdata_pool= data.frame(Y=res$counter.means);
    D=res$D; #Influence functions

    # Estimation
    obsdata_pool$TMLE_group = class
    obsdata_pool$TMLE_group <- relevel(as.factor( obsdata_pool$TMLE_group), ref = ref)
    mod = glm(Y~factor(TMLE_group), family = quasibinomial, data = obsdata_pool);
    coefs = summary(mod)$coefficients[,1];
    X = model.matrix(mod)[1:nb_regimes,]
    B = matrix(coefs, nrow = J);

    CQ=lapply(1:nb_regimes,function(i)(as.matrix(X[i,]))%*%as.matrix((expit(as.matrix(t(X[i,]))%*%B))/(1+exp(as.matrix(t(X[i,]))%*%B))^2)%*%t(as.matrix(((X[i,])))));
    CQ=Reduce('+',CQ)

    Db = matrix(0, nrow = nrow(D[[1]]), ncol = J);
    for(l in 1:nb_regimes){
      Db = Db+as.matrix(D[[l]])%*%solve(CQ);
    }

    se = sqrt(diag(var(Db)/nrow(Db)));
    pvalue <- 2*pnorm(-abs(coefs)/se)
    IClo = coefs- 1.96*se ;
    ICup = coefs + 1.96*se

    results_msm_pltmle = cbind(coefs,se,pvalue,IClo, ICup);
    colnames(results_msm_pltmle) = c("estimate","std.error","pvalue","lower CI", "Upper CI")
    return(results_msm_pltmle)
  }

