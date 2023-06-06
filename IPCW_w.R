#' @title Inverse Probability of Censoring Weights
#' @description Compute unstabilized IPCW
#' @name IPCW_w
#' @param id  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param K number of measuring times.
#' @param censor name of the censoring variable.
#' @param time time variable.
#' @param obsdata observed data in wide format.
#' @author Awa Diop, Denis Talbot
#' @return \item{W_censor_trt}{Unstabilized Inverse Probability of Censoring Weights}
#' @noRd
#' @export IPCW_w
#' @examples
#'dat = gendatTrajMSM(n = 500, Censor = TRUE)
#' V <- c("Age","Sex")
#' L <- c("Hyper","BMI")
#' time <- "Time"
#' A<- "Statins"
#' Censor = "C"
#' wc = IPCW_w(id= "ID", V <- c("Age","Sex"),
#'        L = c("Hyper","BMI"),
#'        time ="Time",
#'        A=c("Statins"),Censor = "C",obsdata = dat)
#'  summary(wc)

IPCW_w<-function(id,V,L,A,K,censor,time,obsdata){
  args <- match.call()
  if(class(obsdata) != "data.frame") stop("Convert obsdata into a data.frame")
  if(!("id" %in% names(args)))       stop("id not specified")
  if(!("V" %in% names(args)))        stop("Baseline variables V not specified")
  if(!("L" %in% names(args)))        stop("Time-varying variables L not specified")
  if(!("A" %in% names(args)))        stop("Time-varying Exposure A not specified")
  if(!("K" %in% names(args)))        stop("Number of measuring times not specified")
  if(!("time" %in% names(args)))     stop("time variable not specified")
  if(!("censor" %in% names(args)))   stop("censor variable not specified")
  if(!("obsdata" %in% names(args)))  stop("Observed Data not specified")

  time <- 1:K
  varA <- sapply(1:length(A), function(x) paste0(A[x],time[1:K]))
  varL <- sapply(1:length(L), function(x) paste0(L[x],time[1:K]))
  varV <- V
  varC <- sapply(1:length(censor), function(x) paste0(L[x],time[1:K]))

  W_trt_temp <- matrix(numeric(0),nrow=nrow(obsdata),ncol=length(varA))
  W_censor_temp <- matrix(numeric(0),nrow=nrow(obsdata),ncol=length(varA))

  for(i in 2:length(varA)){

    sf1 <- paste0(varA[1:(i-1)],collapse  = "+")
    sf2 <- paste0(unlist(varL[1:i,]),collapse  = "+")
    sf3 <- paste0(unlist(varV),collapse  = "+")

    #denominator
    form_denom_censor <- as.formula(paste(paste0(varC[i],"~",sf1),paste0(unlist(varL[1:(i-1),]),collapse  = "+"),sf3, sep = "+"))
    form_denom_trt <- as.formula(paste(paste0(varA[i],"~",sf1),sf2,sf3, sep = "+"))

    C.name <- as.name(varC[i-1])
    fit_denom_censor <- glm(form_denom_censor, family=binomial(link="logit"), data = subset(obsdata, eval(C.name)==0))
    ps_denom_censor = rep(NA,nrow(obsdata))
    ps_denom_censor[obsdata[,varC[i-1]]==0] <- predict(fit_denom_censor, type = "response")

    fit_denom_trt <- glm(form_denom_trt, family=binomial(link="logit"), data = subset(obsdata, eval(C.name)==0))
    ps_denom_trt = rep(NA,nrow(obsdata))
    ps_denom_trt[obsdata[,varC[i-1]]==0]  <- predict(fit_denom_trt, type = "response")

    #unstabilized weights
    W_censor_temp[,i] <- 1/(1-ps_denom_censor)
    W_trt_temp[,i] <- (1-obsdata[,varA[i]])/(1-ps_denom_trt) +
      (obsdata[,varA[i]])/ps_denom_trt
  }
  #t=1

  #denominator
  form_denomt1_censor <- as.formula(paste0(paste0(varC[1],"~",varA[1],"+"),
                                           paste0(unlist(varL[1,]),collapse  = "+"),
                                           "+",paste0(unlist(varV),collapse  = "+")))
  form_denomt1_trt <- as.formula(paste0(varA[1],"~",
                                        paste0(unlist(varL[1,]),collapse  = "+"),
                                        "+",paste0(unlist(varV),collapse  = "+")))

  fit_denomt1_censor <- glm(form_denomt1_censor, family=binomial(link="logit"),data=obsdata)
  ps_denomt1_censor <- predict(fit_denomt1_censor, type = "response")

  fit_denomt1_trt <- glm(form_denomt1_trt, family=binomial(link="logit"),data=obsdata)
  ps_denomt1_trt <- predict(fit_denomt1_trt, type = "response")

  #t=1 unstabilized weights
  W_censor_temp[,1] <- 1/(1-ps_denomt1_censor)
  W_trt_temp[,1] <- (1-obsdata[,varA[1]])/(1-ps_denomt1_trt)+
    (obsdata[,varA[1]])/ps_denomt1_trt

  W_censor_trt <- matrix(numeric(0),nrow=nrow(obsdata),ncol = length(varA))
  W_censor_trt[,1]<-  W_trt_temp[,1]*W_censor_temp[,1]
  for(i in 2:length(varA)){
    W_censor_trt[,i]<-W_trt_temp[,i-1]*W_trt_temp[,i]*W_censor_trt*W_censor_temp[,i-1]
  }

  return(W_censor_trt)
}


