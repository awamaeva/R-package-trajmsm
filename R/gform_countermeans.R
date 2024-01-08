#' @title Counterfactual means for g-formula.
#' @description Get the counterfactual means for the g-formula.
#' @name gform_countermeans
#' @param formula specification of the model for the outcome to be fitted.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param Y outcome variable.
#' @param s number of measuring times.
#' @param time measuring times.
#' @param obsdata observed data in wide format.
#' @returns  \item{list_gform_countermeans}{Counterfactual means obtained with g-formula.}
#' @import e1071
#' @author Awa Diop, Denis Talbot
#' @examples
#' Obswidedata = longtowide(Obsdata = gendatrajMSM(n=500), idvar = "ID", timevar = "Time");
#' formula = paste0("Y.2013~", paste0("Statins.", c(2011:2013),collapse = "+"), "+",
#' paste0("BMI.", c(2011:2013),collapse = "+"),"+",
#' paste0("Hyper.", c(2011:2013),collapse = "+"),"+",
#' "Age.2011 + Sex.2011", collapse = "+")
#' Y = "Y.2013 "
#' A = "Statins."
#' L = c("Hyper.", "BMI.")
#' V = c("Age.", "Sex.")
#' s=3
#' time = c(2011,2012,2013)
#' res.gform = sub_gform(dat=Obswidedata, formula = formula, Y=Y, A=A,L=L,V=V,s=3, time=time)
#' res.gform$counter_means

sub_gform <- function(formula,V,A,L,Y,s, time, obsdata){

  nb_regimes = 2^s; #number of treatment regimes
  # if(nb_regimes > 5000){;
  #   cat("With that number of periods, there is more than 1,000 thousand potential treatment regimes:",nb_regimes);
  #   cont <- readline("Continue? y/n: ");
  #   if(cont == "n" | cont == "no" | cont == "NO") stop("Aborted.");}

  #Treatment regimes
  dat.combn = bincombinations(s);
  list.regimes=lapply(1:nb_regimes,function(x){;
    regime=dat.combn[x,];
    return(regime);
  })

  #To put all individuals under the same treatment regime.
  j = 0; #index for treatment regimes
  #Matrix to store counterfactual means
  res_cm <- matrix(0, nrow = nb_regimes, ncol= 1);
  modQsp1 = glm(formula = formula, family = binomial, data = obsdata);

  # For each treatment regime
  for(regimes in list.regimes){
    j=j+1
    dat0 = obsdata;
    dat0[,paste0(A,1:s)] <- do.call(cbind,lapply(1:s, function(x){
      rgm  <- rep(regimes[x],nrow(dat0))
      return(rgm)}));

    obsdata$Qs = predict(modQsp1, newdata = dat0, type = "res");

    for(i in (s-1):1){

      # Update formula
      sf1 <- paste0(A,1:i, collapse  = "+");
      sf2 <- paste0(unlist(lapply(1:length(L), function(x) paste0(L[x],time[1:i], collapse  = "+"))),collapse  = "+");
      sf3 <- paste0(V,collapse = "+")
      form.up <- as.formula(paste0("Qs~", paste0(sf1, "+", sf3, "+", sf3)));

      #Estimation
      modQs = glm(form.up, family = binomial, data = obsdata);

      # Prediction under a fix treatment regime
      dat0 = obsdata;
      dat0[,paste0(A,1:i)] <- data.frame(do.call(cbind,lapply(1:i, function(x){
        rgm  <- rep(regimes[x],nrow(dat0))
        return(rgm)})));


      obsdata$Qs = predict(modQs, newdata = dat0, type = "res");
    }
    Q0 = mean(obsdata$Qs, na.rm = TRUE);

    res_cm[j,] <- Q0;
  }

  list_gform_countermeans = list(treatment_regimes = list.regimes, counter_means = res_cm);

  return(list_gform_countermeans);
}

