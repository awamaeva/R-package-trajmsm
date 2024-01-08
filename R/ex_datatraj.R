#' @title Data Simulation for Trajectory Analysis
#' @description Example of longitudinal data with three hidden subgroups to perform LCGA/GBTM.
#' @name gendatatraj
#' @param n1 sample size of the first subgroup.
#' @param n2 sample size of the second subgroup.
#' @param n3 sample size of the third subgroup.
#' @param beta01 intercept for the first subgroup.
#' @param beta02 intercept for the second subgroup.
#' @param beta03 intercept for the third subgroup.
#' @param beta11 slope for the first subgroup.
#' @param beta12 slope for the second subgroup.
#' @param beta13 slope for the third subgroup.
#' @param set.seed to add a seed.
#' @return \item{Obsdata}{Wide format data.}
#' @author Awa Diop, Denis Talbot
#' @export gendatatraj
#' @examples Obsdata = gendatatraj()
#' head(Obsdata)

gendatatraj <- function(n1=250, n2=350, n3=400, beta01= -5, beta02=-.15, beta03=-0.01,
                         beta11=-1, beta12=0.5, beta13=5, set.seed=355){
 args <- match.call()

 #Input Check
 if(n1<=0|n2<=0|n3<=0) stop("Specify a number greater than 0")

#Temps 1
A11 <- rbinom(n1,1,1)
A21 <- rbinom(n2,1,1)
A31 <- rbinom(n3,1,1)

#Temps 2
A12 <- rbinom(n1,1,plogis(beta01 +  beta11*A11))
A22 <- rbinom(n2,1,plogis(beta02 +  beta12*A21))
A32 <- rbinom(n3,1,plogis(beta03 +  beta13*A31))

#Temps 3
A13 <- rbinom(n1,1,plogis(beta01 + beta11*A12))
A23 <- rbinom(n2,1,plogis(beta02 + beta12*A22))
A33 <- rbinom(n3,1,plogis(beta03 + beta13*A32))

#Temps 4
A14 <- rbinom(n1,1,plogis(beta01 + beta11*A13))
A24 <- rbinom(n2,1,plogis(beta02 + beta12*A23))
A34 <- rbinom(n3,1,plogis(beta03 + beta13*A33))

#Temps 5
A15 <- rbinom(n1,1,plogis(beta01 + beta11*A14))
A25 <- rbinom(n2,1,plogis(beta02 + beta12*A24))
A35 <- rbinom(n3,1,plogis(beta03 + beta13*A34))

Obsdata <- data.frame(A1= c(A11,A21,A31), A2 = c(A12, A22, A32),
                      A3 = c(A13, A23, A33), A4 = c(A14, A24, A34),
                      A5 = c(A15, A25, A35))

return(Obsdata)
}
