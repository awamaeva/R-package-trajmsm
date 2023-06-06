#' @title Inverse Probability of Censoring Weights
#' @description Compute stabilized IPCW
#' @name IPCW_sw
#' @param id  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param K total length of follow-up.
#' @param censor name of the censoring variable.
#' @param time time variable.
#' @param obsdata observed data in wide format.
#' @author Awa Diop, Denis Talbot
#' @return \item{SW_censor_trt}{Stabilized Inverse Probability of Censoring Weights}
#' @noRd
#' @export IPCW_sw
#' @examples dat = gendatTrajMSM(n = 500, Censor = TRUE)
#' V <- c("Age","Sex")
#' L <- c("Hyper","BMI")
#' time <- "Time"
#' A<- "Statins"
#' Censor = "C"
#' swc = IPCW_sw(id= "ID", V = c("Age","Sex"),
#'         L = c("Hyper","BMI"),
#'         time ="Time",
#'         A=c("Statins"),Censor = "C",obsdata = dat)
#' summary(swc)


IPCW_sw<-function(id,A,L,V,K,Censor,obsdata){
  args <- match.call()
  if(class(obsdata) != "data.frame") stop("Convert obsdata into a data.frame")
  if(!("id" %in% names(args)))       stop("id not specified")
  if(!("V" %in% names(args)))        stop("Baseline variables V not specified")
  if(!("L" %in% names(args)))        stop("Time-varying variables L not specified")
  if(!("A" %in% names(args)))        stop("Time-varying Exposure A not specified")
  if(!("K" %in% names(args)))        stop("Number of measuring times not specified")
  if(!("censor" %in% names(args)))   stop("censor variable not specified")
  if(!("obsdata" %in% names(args)))  stop("Observed Data not specified")

  varA <- sapply(1:length(A), function(x) paste0(A[x],time[1:K]))
  varL <- sapply(1:length(L), function(x) paste0(L[x],time[1:K]))
  varV <- V
  varC <- sapply(1:length(censor), function(x) paste0(L[x],time[1:K]))

  SW_trt_temp <- matrix(numeric(0),nrow=nrow(obsdata),ncol=length(varA))
  SW_censor_temp <- matrix(numeric(0),nrow=nrow(obsdata),ncol=length(varA))

  for(i in 2:length(varA)){

    sf1 <- paste0(varA[1:(i-1)],collapse  = "+")
    sf2 <- paste0(unlist(varL[1:i,]),collapse  = "+")
    sf3 <- paste0(unlist(varV),collapse  = "+")


    #numerator
    form_num_censor <- as.formula(paste(paste0(varC[i],"~",1)))
    form_num_trt <- as.formula(paste(paste0(varA[i],"~",sf1)))

    C.name <- as.name(varC[i-1])
    fit_num_censor <- glm(form_num_censor, family=binomial(link="logit"), data = obsdata)
    ps_num_censor = rep(NA,nrow(obsdata))
    ps_num_censor <- predict(fit_num_censor, type = "response")

    ps_num_trt  <- rep(NA,nrow(obsdata))
    fit_num_trt <- glm(form_num_trt, family=binomial(link="logit"), data = subset(obsdata, eval(C.name)==0) )
    ps_num_trt[obsdata[,varC[i-1]]==0]  <- predict(fit_num_trt, type = "response")

    #denominator
    form_denom_censor <- as.formula(paste(paste0(varC[i],"~",sf1),paste0(unlist(varL[1:(i-1),]),collapse  = "+"),sf3, sep = "+"))
    form_denom_trt <- as.formula(paste(paste0(varA[i],"~",sf1),sf2,sf3, sep = "+"))

    fit_denom_censor <- glm(form_denom_censor, family=binomial(link="logit"), data =obsdata )
    ps_denom_censor = rep(NA,nrow(obsdata))
    ps_denom_censor <- predict(fit_denom_censor, type = "response")

    fit_denom_trt <- glm(form_denom_trt, family=binomial(link="logit"), data = subset(obsdata, eval(C.name)==0))
    ps_denom_trt = rep(NA,nrow(obsdata))
    ps_denom_trt[obsdata[,varC[i-1]]==0]  <- predict(fit_denom_trt, type = "response")

    #stabilized weights
    SW_censor_temp[,i] <- (1-ps_num_censor)/(1-ps_denom_censor)
    SW_trt_temp[,i] <- ((1-obsdata[,varA[i]])*(1-ps_num_trt))/(1-ps_denom_trt) +
      (obsdata[,varA[i]]*ps_num_trt)/ps_denom_trt
  }
  C.name <- as.name(varC[1])
  #t=1
  #numerator
  form_numt1_censor <- as.formula(paste0(paste0(varC[1],"~",1)))
  form_numt1_trt <- as.formula(paste0(varA[1],"~",1))

  fit_numt1_censor <- glm(form_numt1_censor, family=binomial(link="logit"),data = obsdata)
  ps_numt1_censor <- predict(fit_numt1_censor, type = "response")

  fit_numt1_trt <- glm(form_numt1_trt, family=binomial(link="logit"),data=subset(obsdata, eval(C.name)==0))
  ps_numt1_trt = rep(NA,nrow(obsdata))
  ps_numt1_trt[obsdata[,varC[1]]==0]  <- predict(fit_numt1_trt, type = "response")

  #denominator
  form_denomt1_censor <- as.formula(paste0(paste0(varC[1],"~",varA[1],"+"),
                                           paste0(unlist(varL[1,]),collapse  = "+"),
                                           "+",paste0(unlist(varV),collapse  = "+")))
  form_denomt1_trt <- as.formula(paste0(varA[1],"~",
                                        paste0(unlist(varL[1,]),collapse  = "+"),
                                        "+",paste0(unlist(varV),collapse  = "+")))

  fit_denomt1_censor <- glm(form_denomt1_censor, family=binomial(link="logit"),data = obsdata)
  ps_denomt1_censor <- predict(fit_denomt1_censor, type = "response")

  fit_denomt1_trt <- glm(form_denomt1_trt, family=binomial(link="logit"),data = subset(obsdata, eval(C.name)==0))
  ps_denomt1_trt = rep(NA,nrow(obsdata))
  ps_denomt1_trt[obsdata[,varC[1]]==0] <- predict(fit_denomt1_trt, type = "response")

  #t=1 stabilized weights
  SW_censor_temp[,1] <- 1*(1-ps_numt1_censor)/(1-ps_denomt1_censor)
  SW_trt_temp[,1] <- ((1-obsdata[,varA[1]])*(1-ps_numt1_trt))/(1-ps_denomt1_trt)+
    (obsdata[,varA[1]]*ps_numt1_trt)/ps_denomt1_trt

  SW_censor_trt <- matrix(numeric(0),nrow=nrow(obsdata),ncol = length(varA))
  SW_censor_trt[,1]<-  SW_trt_temp[,1]*SW_censor_temp[,1]
  for(i in 2:length(varA)){
    SW_censor_trt[,i]<-SW_trt_temp[,i-1]*SW_trt_temp[,i]*SW_censor_temp[,i]*SW_censor_temp[,i-1]
  }
  return(SW_censor_trt)
}



