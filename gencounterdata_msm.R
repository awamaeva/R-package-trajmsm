gencounterdata_msm <- function(n){
ds.counter = list();
j=0
for(regimes in list(c(0,0,0), c(0,0,1), c(0,1,0), c(0,1,1),
                    c(1,0,0), c(1,0,1), c(1,1,0), c(1,1,1))){
  j=j+1
  #Window 1
  Y1<-rep(0,n)#1=CVD, 0= no event;no CVD
  V=rbinom(n=n,size=1,prob=0.5)
  A1 = rep(regimes[1],n); #We suppose there is a second intervention that prevents the event of happening.
  mu1<-expit(1+V+0.5*A1);
  L1 <-rbinom(n,1,mu1);


  A2 = rep(regimes[2],n);
  mu2<-expit(1+L1+0.5*A2)
  L2<-rbinom(n=n,size=1,prob=mu2)


  A3 = rep(regimes[3],n);
  mu3<-expit(1+L2+0.5*A3)
  L3<-rbinom(size=1,prob=mu3,n=n)

  Y4<-rep(0,n)
  r4<-expit(1+V-0.7*L3-0.5*A3)
  Y4<-rbinom(n=n,size=1,prob=r4)

  ds.counter[j] =list(cbind(A1=A1, A2=A2, A3=A3, Y = Y4));
}
return(ds.counter)
}
