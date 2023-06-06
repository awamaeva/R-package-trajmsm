#' @title Data Simulation for the LCGA-MSM.
#' @description Example of longitudinal data with three hidden subgroups to perform LCGA/GBTM.
#' @name gen_obsdata_msm
#' @param seed to set a seed for replicability.
#' @param n sample size.
#' @param t follow-up time.
#' @export gen_obsdata_msm
#' @return \item{dat}{data generated for the LCGA-MSM}
#' @noRd

gen_obsdata_msm <- function(n, t, seed){

  set.seed = seed
  expit = plogis ;
    V <- rnorm(n)
    L1 =rbinom(n,1,expit(V));
    A1 = rbinom(n, 1, expit(.75 - 0.25*L1 -0.5*V));

    #Starting Values of A, L et Y
    L<-matrix(0,n,t)
    A<-matrix(0,n,t)
    A[,1]<-A1
    L[,1]<-L1
    #Generate data for time varying exposure and covariates
    for(i in 2:(t)){
      L[,i]<-rbinom(n,1,expit(0.25*L[,i-1]-(.5)*A[,i-1]));
      A[,i]<- rbinom(n, 1, expit(.5 +(.25^(t))*A[,i-1] + (.5^(t))*L[,i] -(.5^(t))*A[,i-1]*L[,i]))
    }
    colnames(A)<-paste0("A",1:t)
    colnames(L)<-paste0("L",1:t)

    #Generate the binary outcome
    Y<- rbinom(n,1,expit(0.05-A[,t] - 0.1*rowSums(A[,2:(t-1),drop=FALSE]) - 0.2*A[,1] + L[,t] + 0.5*rowSums(L[,2:(t-1),drop=FALSE]) + 0.2*L[,1]));
   return(dat = data.frame(id = 1:n,V,A,L,Y))
}
