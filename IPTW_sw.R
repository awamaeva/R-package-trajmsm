#' @title Unstabilized Inverse of Probability of Treatment Weighting
#' @description Compute stabilized IPTW
#' @name IPTW_sw
#' @param id  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param K number of measuring times.
#' @param obsdata observed data in wide format.
#' @return \item{SW}{Stabilized inverse of treatment probabilities}
#' @noRd
#' @export IPTW_sw
#' @author Awa Diop, Denis Talbot

IPTW_sw<-function(id,A,L,V,K,time,obsdata){
  args <- match.call()
  if(class(obsdata) != "data.frame") stop("Convert obsdata into a data.frame")
  if(!("id" %in% names(args)))       stop("id not specified")
  if(!("V" %in% names(args)))        stop("Baseline variables V not specified")
  if(!("L" %in% names(args)))        stop("Time-varying variables L not specified")
  if(!("A" %in% names(args)))        stop("Time-varying Exposure A not specified")
  if(!("K" %in% names(args)))        stop("Number of measuring times not specified")
  if(!("time" %in% names(args)))        stop("Time not specified")
  if(!("obsdata" %in% names(args)))  stop("Observed Data not specified")

  #varA <- if(length(A)==1){paste0("A.",time[1:K])}else{sapply(1:length(A), function(x) paste0("A.",time[1:K]))}
  #varL <- if(length(L)==1){paste0(L,time[1:K])}else{sapply(1:length(L), function(x) paste0(L[x],time[1:K]))}

  varA <- c("A.1", "A.2", "A.3", "A.4", "A.5")
  varL <- c("L.1", "L.2", "L.3", "L.4", "L.5")
  varV <- V


  SW_temp <- matrix(numeric(0),nrow=nrow(obsdata),ncol=length(varA))

  for(i in 2:length(varA)){
    sf1 <- if(length(A)==1){paste0(varA[1:(i-1)],collapse  = "+")}else{paste0(varA[1:(i-1)],collapse  = "+")}
    sf2 <- if(length(L)==1){paste0(varL[1:i],collapse  = "+")}else{paste0(unlist(varL[1:i,]),collapse  = "+")}
    sf3 <-  paste0(varV,collapse  = "+")

    form_num<- as.formula(paste0(varA[i],"~",sf1))
    form_denom <- as.formula(paste(paste0(varA[i],"~",sf1),sf2,sf3, sep = "+"))

    fit_num <- glm(form_num, family=binomial(link="logit"), data = obsdata)
    fit_denom <- glm(form_denom, family=binomial(link="logit"), data = obsdata)

    ps_num  <- ifelse(!is.na(obsdata[,varA[i]]),predict(fit_num, type = "response"),NA)
    ps_denom  <- ifelse(!is.na(obsdata[,varA[i]]),predict(fit_denom, type = "response"),NA)

    SW_temp[,i] <- ((1-obsdata[,varA[i]])*(1-ps_num))/(1-ps_denom) + ((obsdata[,varA[i]])*ps_num)/ps_denom
  }
  #t=1
  if(length(L)==1){varL1= paste0(varL[1],collapse  = "+")}else{varL1 = paste0(unlist(varL[1,,]),collapse  = "+")}

  form_numt1 <- as.formula(paste0(varA[1],"~", 1))
  form_denomt1 <- as.formula(paste0(varA[1],"~",
                                    paste0(varL1,collapse  = "+"),
                                    "+",paste0(unlist(varV),collapse  = "+")))

  fit_denomt1 <- glm(form_denomt1, family=binomial(link="logit"),data=obsdata)
  fit_numt1 <- glm(form_numt1, family=binomial(link="logit"),data=obsdata)

  ps_numt1 <-  ifelse(!is.na(obsdata[,varA[1]]),predict(fit_numt1, type = "response"),NA)
  ps_denomt1 <-  ifelse(!is.na(obsdata[,varA[1]]),predict(fit_denomt1, type = "response"),NA)

  SW_temp[,1] <- ((1-obsdata[,varA[1]])*(1-ps_numt1))/(1-ps_denomt1)+((obsdata[,varA[1]])*ps_numt1)/ps_denomt1
  SW <- t(apply(SW_temp,1,cumprod))
  return(SW)
}
