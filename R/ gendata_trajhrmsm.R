
 gendata_trajhrmsm<-function(n,seed = 945){
 set.seed = seed;
 expit = plogis;
  #Window 1
  Y1<-rep(0,n)#1=CVD, 0= no event;no CVD
  V=rbinom(n=n,size=1,prob=0.5)
  p1<-plogis(-0.5+2.5*V);
  A1 <-rbinom(n, 1, p1); #We suppose there is a second intervention that prevents the event of happening.
  mu1<-expit(1+V+0.5*A1);
  L1 <-rbinom(n,1,mu1);
  Y2<-rep(0,n);

  p2<-expit(-.5+0.5*A1+V+1.2*L1);
  A2<-rbinom(n=n,size=1,prob=p2)
  mu2<-expit(1+L1+0.5*A2)
  L2<-rbinom(n=n,size=1,prob=mu2)
  Y3<-rep(0,n)

  A3<-rep(0,length=n)
  p3<-expit(-.5+0.5*A2+V+1.2*L2)
  A3<-rbinom(n=n,size=1,prob=p3)

  mu3<-expit(1+L2+0.5*A3);
  L3<-rbinom(size=1,prob=mu3,n=n)

  Y4<-rep(0,n)
  r4<-expit(0.5*V-0.5*L3-0.25*A3);
  Y4[Y3==0]<-rbinom(n=sum(Y3==0),size=1,prob=r4[Y3==0])
  Y4[Y3==1]<-1

  #Window 2
  p4<-expit(-.5+0.5*A3+V+1.2*L3)
  A4<-rbinom(n=n,size=1,prob=p4)
  mu4<-expit(1+L3+0.5*A4);
  L4<-rbinom(prob=mu4,size=1,n=n)

  Y5<-rep(0,n)
  r5<-expit(0.5*V-0.5*L4-0.25*A4);
  Y5[Y4==0]<-rbinom(n=sum(Y4==0),size=1,prob=r5[Y4==0])
  Y5[Y4==1]<-1

#Window 3
p5<-expit(-.5+0.5*A4+V+1.2*L4)
A5<-rbinom(n=n,size=1,prob=p5)
mu5<-expit(1+L4+0.5*A5)
L5<-rbinom(prob=mu5,size=1,n=n)

Y6<-rep(0,n)
r6<-expit(0.5*V-0.5*L5-0.25*A5)
Y6[Y5==0]<-rbinom(n=sum(Y5==0),size=1,prob=r6[Y5==0])
Y6[Y5==1]<-1

A4[Y4==1]<-NA
L4[Y4==1]<-NA
A5[Y5==1]<-NA
L5[Y5==1]<-NA

return(data.frame(id = 1:n, A = c(A1,A2,A3,A4,A5), L = c(L1,L2,L3,L4,L5), Y = c(Y2,Y3,Y4,Y5,Y6),V, time = rep(1:5, each = n)))
}
