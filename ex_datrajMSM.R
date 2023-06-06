gendatrajMSM <-function(n, Censor = FALSE, format = c("long","wide")){
  options(warn = -1)
  if(Censor == FALSE & format == "long"){

    ID = rep(1:n, each = 3)
    Time=rep(c(2011,2012,2013),n)
    Age = rbinom(n,1,0.5)
    Sex = rbinom(n,1,0.5)
    Statins1 = rbinom(n,1,plogis(-0.5+0.2*Age+0.1*Sex))
    Hyper1 = rbinom(n,1,0.5)
    BMI1 = rbinom(n,1,0.5)

    Statins2 = rbinom(n,1,plogis(-1+0.2*Age+0.1*Sex + 0.5*Statins1 + 0.2*Hyper1 + 0.1*BMI1))
    Hyper2 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.5*Statins1 + Hyper1 + 0.1*BMI1))
    BMI2 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.1*Statins1 + 0.2*Hyper1 + BMI1))

    Statins3 = rbinom(n,1,plogis(-1.5+0.2*Age+0.1*Sex + 0.5*Statins2 + 0.2*Hyper2 + 0.1*BMI2))
    Hyper3 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.5*Statins2 + Hyper2 + 0.1*BMI2))
    BMI3 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.1*Statins2 + 0.2*Hyper2 + BMI2))

   Y = rbinom(n,1, plogis(0.05-0.1*Statins1 + 0.2*Hyper1   +    0.1*BMI1
                          - 0.15*Statins2 + 0.3*Hyper2 +  0.2*BMI2
                          - 0.2*Statins3 + 0.4*Hyper3  +  0.3*BMI3))
    dat <- data.frame(ID=ID,
                      Statins = c(Statins1,Statins2,Statins3),
                      Hyper = c(Hyper1,Hyper2,Hyper3),
                      BMI = c(BMI1,BMI2,BMI3),
                      Age = Age,
                      Sex = Sex,
                      Time=Time,Y=Y)
  }
if(Censor == TRUE & format == "long"){
ID = rep(1:n, each = 3)
Time=rep(c(2011,2012,2013),n)

Statins1 = rbinom(n,1,plogis(0.2*Age+0.1*Sex))
Hyper1 = rbinom(n,1,0.5)
BMI1 = rbinom(n,1,0.5)

Statins2 = rbinom(n,1,plogis(0.2*Age+0.1*Sex + 0.5*Statins1 + 0.2*Hyper1 + 0.1*BMI1))
Hyper2 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.5*Statins1 + Hyper1 + 0.1*BMI1))
BMI2 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.1*Statins1 + 0.2*Hyper1 + BMI1))

Statins3 = rbinom(n,1,plogis(0.2*Age+0.1*Sex + 0.5*Statins2 + 0.2*Hyper2 + 0.1*BMI2))
Hyper3 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.5*Statins2 + Hyper2 + 0.1*BMI2))
BMI3 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.1*Statins2 + 0.2*Hyper2 + BMI2))

Y = rbinom(n,1, plogis(-0.1*Statins1 + 0.2*Hyper1   +    0.1*BMI1
                       - 0.15*Statins2 + 0.3*Hyper2 +  0.2*BMI2
                       - 0.2*Statins3 + 0.4*Hyper3  +  0.3*BMI3))

C1= rbinom(n,1,plogis(-5-0.5*Statins1 + 0.5*Hyper1 + 0.2*BMI1))
C2= rbinom(n,1,plogis(-5-0.5*Statins2 + 0.5*Hyper2 + 0.2*BMI2))
C3= rbinom(n,1,plogis(-5-0.5*Statins3 + 0.5*Hyper3 + 0.2*BMI3))

C2[C1==1]<-1
C3[C2==1]<-1

Statins2[C1==1]<-NA
Statins3[C2==1]<-NA

Hyper2[C1==1]<-NA
Hyper3[C2==1]<-NA

BMI2[C1==1]<-NA
BMI3[C2==1]<-NA


dat <- data.frame(ID=ID,
                  Statins = c(Statins1,Statins2,Statins3),
                  Hyper = c(Hyper1,Hyper2,Hyper3),
                  BMI = c(BMI1,BMI2,BMI3),
                  Age = Age,
                  Sex = Sex,
                  Time=Time,
                  C = c(C1,C2,C3), Y = Y)

}

  if(Censor == FALSE & format == "wide"){

    ID = rep(1:n, each = 3)
    Time=rep(c(2011,2012,2013),n)
    Age = rbinom(n,1,0.5)
    Sex = rbinom(n,1,0.5)
    Statins1 = rbinom(n,1,plogis(0.2*Age+0.1*Sex))
    Hyper1 = rbinom(n,1,0.5)
    BMI1 = rbinom(n,1,0.5)

    Statins2 = rbinom(n,1,plogis(0.2*Age+0.1*Sex + 0.5*Statins1 + 0.2*Hyper1 + 0.1*BMI1))
    Hyper2 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.5*Statins1 + Hyper1 + 0.1*BMI1))
    BMI2 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.1*Statins1 + 0.2*Hyper1 + BMI1))

    Statins3 = rbinom(n,1,plogis(0.2*Age+0.1*Sex + 0.5*Statins2 + 0.2*Hyper2 + 0.1*BMI2))
    Hyper3 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.5*Statins2 + Hyper2 + 0.1*BMI2))
    BMI3 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.1*Statins2 + 0.2*Hyper2 + BMI2))

    Y = rbinom(n,1, plogis(-0.1*Statins1 + 0.2*Hyper1   +    0.1*BMI1
                           - 0.15*Statins2 + 0.3*Hyper2 +  0.2*BMI2
                           - 0.2*Statins3 + 0.4*Hyper3  +  0.3*BMI3))
    dat <- data.frame(ID=ID,
                      Statins1 = Statins1,Statins2 = Statins2, Statins3 = Statins3,
                      Hyper1 = Hyper1, Hyper2 = Hyper2, Hyper3 = Hyper3,
                      BMI1 = BMI1, BMI2 = BMI2, BMI3 = BMI3,
                      Age = Age,
                      Sex = Sex,
                      Time=Time,Y=Y)
  }
  if(Censor == TRUE & format == "wide"){
    ID = rep(1:n, each = 3)
    Time=rep(c(2011,2012,2013),n)

    Statins1 = rbinom(n,1,plogis(0.2*Age+0.1*Sex))
    Hyper1 = rbinom(n,1,0.5)
    BMI1 = rbinom(n,1,0.5)

    Statins2 = rbinom(n,1,plogis(0.2*Age+0.1*Sex + 0.5*Statins1 + 0.2*Hyper1 + 0.1*BMI1))
    Hyper2 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.5*Statins1 + Hyper1 + 0.1*BMI1))
    BMI2 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.1*Statins1 + 0.2*Hyper1 + BMI1))

    Statins3 = rbinom(n,1,plogis(0.2*Age+0.1*Sex + 0.5*Statins2 + 0.2*Hyper2 + 0.1*BMI2))
    Hyper3 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.5*Statins2 + Hyper2 + 0.1*BMI2))
    BMI3 = rbinom(n,1,plogis(0.5*Age+0.7*Sex - 0.1*Statins2 + 0.2*Hyper2 + BMI2))

    Y = rbinom(n,1, plogis(-0.1*Statins1 + 0.2*Hyper1   +    0.1*BMI1
                           - 0.15*Statins2 + 0.3*Hyper2 +  0.2*BMI2
                           - 0.2*Statins3 + 0.4*Hyper3  +  0.3*BMI3))

    C1= rbinom(n,1,plogis(-5-0.5*Statins1 + 0.5*Hyper1 + 0.2*BMI1))
    C2= rbinom(n,1,plogis(-5-0.5*Statins2 + 0.5*Hyper2 + 0.2*BMI2))
    C3= rbinom(n,1,plogis(-5-0.5*Statins3 + 0.5*Hyper3 + 0.2*BMI3))

    C2[C1==1]<-1
    C3[C2==1]<-1

    Statins2[C1==1]<-NA
    Statins3[C2==1]<-NA

    Hyper2[C1==1]<-NA
    Hyper3[C2==1]<-NA

    BMI2[C1==1]<-NA
    BMI3[C2==1]<-NA


    dat <- data.frame(ID=ID,
                      Statins1 = Statins1,Statins2 = Statins2, Statins3 = Statins3,
                      Hyper1 = Hyper1, Hyper2 = Hyper2, Hyper3 = Hyper3,
                      BMI1 = BMI1, BMI2 = BMI2, BMI3 = BMI3,
                      C1 = C1, C2 = C2, C3 = C3,
                      Age = Age,
                      Sex = Sex,
                      Time=Time,Y=Y)

  }
  return(dat)
}
Adat_long <- widetolong(obsdata = obsdata_msm[, c("A1","A2","A3")], varying = 1:3)
