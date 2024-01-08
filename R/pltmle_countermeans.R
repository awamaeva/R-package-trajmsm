#' @title Counterfactual means for a Pooled LTMLE
#' @description function to estimate counterfactual means for a pooled LTMLE.
#' @name sub_pltmle
#' @param formula specification of the model for the outcome to be fitted.
#' @param id  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param Y outcome variable.
#' @param s number of measuring times.
#' @param J number of trajectory groups.
#' @param time measuring times.
#' @param obsdata observed data in wide format.
#' @returns  \item{list_pltmle_countermeans}{Counterfactual means and influence functions with the pooled ltmle.}
#' \item{D}{Influence functions}
#' @import e1071
#' @examples
#' obsdata_msm <- genobsdata_msm(n=2500)
#' V <- "V"
#' L <- "L"
#' time <- 1:3
#' A <- "A"
#' Y <- "Y"
#' @author Awa Diop, Denis Talbot


sub_pltmle<-function(formula,Y,A,L,V,s,J, time,timevar,idvar,obsdata,Traj){
  varA <-  paste0(A,time[1:s])
  D = NULL
  j = 0;
  obsdata0.all=list()
  obsdata.all=list()
  nb_regimes = 2^s #number of treatment regimes
  logit = qlogis;
  expit = plogis;

  if(nb_regimes > 5000) {
    cat("With that number of periods, there is more than 5,000 thousand potential treatment regimes:",nb_regimes);
    cont <- readline("Continue? y/n: ");
    if(cont == "n" | cont == "no" | cont == "NO") stop("Aborted.")
  }
  #Treatment regimes
  dat.combn = bincombinations(s);
  list.regimes=lapply(1:nb_regimes,function(x){;
    regime = dat.combn[x,];
    return(regime);
  })

  j = 0
  for(regimes in list.regimes){
    j=j+1
    obsdata0 = obsdata;
    obsdata0[,paste0(A,time[1:s])] <- do.call(cbind,lapply(1:s, function(x){
      rgm  <- rep(regimes[x],nrow(obsdata0))
      return(rgm)}));

    obsdata0.all[j]=list(obsdata0)
    obsdata.all[j]=list(obsdata)
  }
  obsdata.all2  = data.frame(do.call(rbind,obsdata.all))
  obsdata0.all2 = data.frame(do.call(rbind,obsdata0.all))

  #Model
  modQs = glm(formula = formula, family = binomial, data = obsdata);
  #for t=3 predict the outcome for all different regimes of treatment
  Qs = lapply(1:nb_regimes,function(i)predict(modQs, newdata = obsdata0.all[[i]], type = "res"));

  #Compute the weights for all different regimes of treatment
  W = IPW(numerator = "unstabilized", id = idvar , V = V,L =L,A = A, C = FALSE,
            censor = censor, K = s,time = time,obsdata = obsdata);

  Hs.all=list()
  j=0;
  for(regimes in list.regimes){
    j=j+1
    Hs.all[j] = list(as.matrix(1*(rowSums(sapply(1:s,function(x)
      1*(obsdata[,varA[x]] == unlist(regimes)[x])))==s)*W[,s])%*%t(as.matrix(Traj[j,1:J])))
  }

  Hs=do.call(rbind,Hs.all)
  #  Update the risk for each regime of treatment
  ## First estimate the vector of fluctuation error

  modEs = glm(as.formula(paste0(Y, "~", "-1 + offset(logit(unlist(Qs))) + Hs")), family = binomial, data = obsdata.all2);
  coef_Es <- coef(modEs)
  coef_Es<- ifelse(is.na(coef_Es),0, coef_Es)
  Qstar = lapply(1:nb_regimes,function(x) expit(logit(Qs[[x]]) +
                                                  as.numeric(t(t(as.matrix(coef_Es))%*%t(Hs.all[[x]])))));

  ##Influence curve for each regime of treatment
  D_list <- list()
  Ds = lapply(1:nb_regimes,function(x)sapply(1:J,function(y)
    (Hs.all[[x]][,y]*(obsdata[,Y] - Qstar[[x]]))))

  D_list[s]<-list(Ds)
  # initialization of a progress bar
  pb <- txtProgressBar(min = 1, max = s-1, style = 3)

  for(i in (s-1):1){
    #For t=s-1 repeat the same process as for t=s
    sf1 <- paste0(A,time[1:i], collapse  = "+");
    sf2 <- paste0(unlist(lapply(1:length(L), function(x) paste0(L[x],time[1:i], collapse  = "+"))),collapse  = "+");
    sf3 <- paste0(V,collapse = "+")


    modQs = lapply(1:nb_regimes,function(x){
      qstar = Qstar[[x]];
      form.up <- as.formula(paste("qstar ~",paste0(sf1, "+", sf2, "+", sf3)));
      mod_temp = glm(form.up, family = binomial, data = obsdata);
      return(mod_temp);
    });

    Qs = lapply(1:nb_regimes,function(x)predict(modQs[[x]], newdata = obsdata0.all[[x]], type = "res"));

    Hs.all=list()
    j=0
    for(regimes in list.regimes){
      j=j+1
      Hs.all[j] = list(as.matrix(1*(rowSums(sapply(1:i,function(x) 1*(obsdata[,varA[x]] == unlist(regimes)[x])))==i)*W[,i])%*%t(as.matrix(Traj[j,1:J])))
    }

    Hs = do.call(rbind,Hs.all)

    modEs = glm(unlist(Qstar) ~-1+offset(logit(unlist(Qs)))+Hs, family = binomial,data = obsdata.all2);
    coef_Es <- coef(modEs)
    coef_Es<- ifelse(is.na(coef_Es),0, coef_Es)
    Qstarm1 = lapply(1:nb_regimes,function(x)expit(logit(Qs[[x]]) +
                                                     as.numeric(t(t(as.matrix( coef_Es))%*%t(Hs.all[[x]])))));

    ##Influence curve for each regime of treatment
    Ds = lapply(1:nb_regimes,function(x) sapply(1:J,function(y)
      as.matrix(Hs.all[[x]][,y]*(Qstar[[x]]-Qstarm1[[x]]))))
    D_list[i]<-list(Ds)
    Qstar = Qstarm1;
    setTxtProgressBar(pb, i)
  }
 close(pb)
  Q0 = unlist(Qstar);

  D = lapply(1:nb_regimes,function(x) as.matrix(Reduce('+',lapply(1:J, function(y){
    comps = as.matrix(D_list[[y]][[x]] +  as.numeric(Qstar[[x]] - mean(Qstar[[x]])))
    return(comps)
  }))))

  obsdata0.all2$Y = Q0
  obsdata0.all2$atraj <- apply(obsdata0.all2[,varA],1,paste0,collapse="")
  obsdataT2 = aggregate(as.formula(paste0(Y,"~", "atraj")), FUN = mean,data = obsdata0.all2)
  obsdataT2[,Y] = as.numeric(as.character(obsdataT2[,Y]))


  list_pltmle_countermeans = list(counter.means = obsdataT2[,Y],D=D)
  return(list_pltmle_countermeans)
}
