#' @title Unstabilized Inverse of Probability of Treatment Weighting
#' @description Compute unstabilized IPTW
#' @name IPTW_w
#' @param id  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param K number of measuring times.
#' @param obsdata observed data in wide format.
#' @return \item{W}{Unstabilized inverse of treatment probabilities}
#' @noRd
#' @export IPTW_w
#' @examples
#' obsdata_msm <- genobsdata_msm(n=2500)
#' id = "id"
#' V <- "V"
#' L <- "L"
#' A <- "A"
#' K = 3
#' time = 1:3
#' w = IPTW_w(id= "ID", V = V,
#'      A = A,
#'      L = L,
#'      K = K,
#'      time = time,
#'      obsdata = obsdata_msm)
#' summary(w)
#' @author Awa Diop, Denis Talbot

IPTW_w<-function(id,V,L,A,K,obsdata, time){

  args <- match.call()
  if(class(obsdata) != "data.frame") stop("Convert obsdata into a data.frame")
  if(!("id" %in% names(args)))       stop("id not specified")
  if(!("V" %in% names(args)))        stop("Baseline variables V not specified")
  if(!("L" %in% names(args)))        stop("Time-varying variables L not specified")
  if(!("A" %in% names(args)))        stop("Time-varying Exposure A not specified")
  if(!("K" %in% names(args)))        stop("Number of measuring times not specified")
  if(!("obsdata" %in% names(args)))  stop("Observed Data not specified")

  varA <- sapply(1:length(A), function(x) paste0(A[x],time[1:K]))
  varL <- sapply(1:length(L), function(x) paste0(L[x],time[1:K]))
  varV <- V

  W_temp <- matrix(numeric(0),nrow=nrow(obsdata),ncol=length(varA))
  for(i in 2:length(varA)){
    sf1 <- paste0(sapply(1:length(A), function(x) paste0(A[x],time[1:(i-1)])), collapse  = "+")
    sf2 <- paste0(sapply(1:length(L), function(x) paste0(L[x],time[1:i])), collapse  = "+")
    sf3 <-  paste0(V,collapse="+")

    form_denom <- as.formula(paste(paste0(varA[i],"~",sf1),sf2,sf3, sep = "+"))
    fit_denom <- glm(form_denom, family=binomial(link="logit"), data = obsdata)
    ps_denom  <- predict(fit_denom, type = "response")

    W_temp[,i] <- (1-obsdata[,varA[i]])/(1-ps_denom) + (obsdata[,varA[i]])/ps_denom
  }
  #t=1
  form_denomt1 <- as.formula(paste0(varA[1],"~",
                                    paste0(unlist(varL[1,]),collapse  = "+"),
                                    "+",paste0(unlist(varV),collapse  = "+")))

  fit_denomt1 <- glm(form_denomt1, family=binomial(link="logit"),data=obsdata)
  ps_denomt1 <- predict(fit_denomt1, type = "response")
  W_temp[,1] <- (1-obsdata[,varA[1]])/(1-ps_denomt1)+(obsdata[,varA[1]])/ps_denomt1

  W <- t(apply(W_temp,1,cumprod))

  return(W)
}
