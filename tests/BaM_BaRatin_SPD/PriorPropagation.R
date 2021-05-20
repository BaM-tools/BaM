continuity<-function(x,a1,b1,c1,a2,b2,c2){
  z=a1*(x-b1)^c1-a2*(x-b2)^c2
  return(z)
}

P1<-function(nsim,a1.m,a1.s,a2.m,a2.s,c1.m,c1.s,c2.m,c2.s, # fixed parameters
             b1.m,b1.s,b2.m,b2.s,dg.s,dl.s, # varying parameters
             dg.m=0*dg.s,dl.m=0*dl.s){
  n=length(dg.s)+1
  b1=matrix(NA,nsim,n)
  b2=b1
  k2=b1
  for(i in 1:nsim){
    a1=rnorm(1,a1.m,a1.s)
    a2=rnorm(1,a2.m,a2.s)
    c1=rnorm(1,c1.m,c1.s)
    c2=rnorm(1,c2.m,c2.s)
    b1[i,1]=rnorm(1,b1.m,b1.s)
    b2[i,1]=rnorm(1,b2.m,b2.s)
    dg=rnorm(n-1,dg.m,dg.s)
    dl=rnorm(n-1,dl.m,dl.s)
    b1[i,2:n]=b1[i,1]-(cumsum(dg)+cumsum(dl))
    b2[i,2:n]=b2[i,1]-cumsum(dg)
    # find kappas
    for(j in 1:n){
      w=tryCatch(uniroot(continuity,c(max(b1[i,j],b2[i,j]),10),a1,b1[i,j],c1,a2,b2[i,j],c2),
                 error=function(e) list(root=NA))
      k2[i,j]=w$root
    }    
  }
  return(list(b1=b1,b2=b2,k1=b1,k2=k2))
}

############
# Meyras case study
nperiod=5
a1.m=14.17;a1.s=5.01
a2.m=26.52;a2.s=8.4
c1.m=1.5;c1.s=0.025
c2.m=1.67;c2.s=0.025
b1.m=-0.6;b1.s=0.5
b2.m=-0.6;b2.s=0.5
dg.s=rep(0.25,nperiod-1);dl.s=rep(0.25,nperiod-1)
# Propagation
w=P1(10000,a1.m,a1.s,a2.m,a2.s,c1.m,c1.s,c2.m,c2.s,b1.m,b1.s,b2.m,b2.s,dg.s,dl.s)
k=cbind(w$k1,w$k2)
apply(k,2,mean,na.rm=T)
apply(k,2,sd,na.rm=T)
nInfer=2+2+3+nperiod*2+2
C=diag(nInfer)
W=cor(k,use="complete.obs")
C[1:5,1:5]=W[1:5,1:5]
C[8:12,8:12]=W[6:10,6:10]
C[1:5,8:12]=W[1:5,6:10]
C[8:12,1:5]=W[6:10,1:5]
image(C)
write(C,"PriorCorrelation.txt",ncolumns=nInfer)