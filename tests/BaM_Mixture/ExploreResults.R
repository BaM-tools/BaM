nS=5
nEvt=27
nElt=13
missingSource=T

if(missingSource){
  folder='WithMissingSource/'
  prefix='Missing_'
} else {
  folder='WithoutMissingSource/'
  prefix=''
}
mcmc=read.table(paste(folder,'Results_MCMC_Cooked.txt',sep=''),header=T)
res=read.table(paste(folder,'Results_Residuals.txt',sep=''),header=T)

# Contributions
pdf(paste(folder,prefix,'contributions.pdf',sep=''),width=12,height=12)
par(mfrow=c(6,5),mar=c(2,2,1,1))
for(i in 1:nEvt){
  nam=c('Ai','Ar','Bo','Fi','Gu')
  if(missingSource){
    col=( ((i-1)*nS+1) : (i*nS) )
    col=c(col,nS*nEvt+nElt+2+1+i)
    nam=c(nam,'Mi')
  } else {
    col=( ((i-1)*(nS-1)+1) : (i*(nS-1)) )
    col=c(col,(nS-1)*nEvt+2+1+i)
  }
  boxplot(mcmc[,col],col=1+(1:(nS+1)),
          names=nam,
          ylim=c(0,1))
}
dev.off()

# Fit
pdf(paste(folder,prefix,'fit.pdf',sep=''),width=12,height=12)
par(mfrow=c(6,5),mar=c(2,2,1,1))
for(i in 1:nEvt){
  row=( ((i-1)*nElt+1) : (i*nElt) )
  plot(res$Y1_obs[row],type='b',pch=19)
  lines(res$Y1_sim[row],col='red')
 }
dev.off()

# Residuals
pdf(paste(folder,prefix,'residuals.pdf',sep=''),width=12,height=3)
par(mfrow=c(1,3))
plot(res$Y1_res,xlab='index',ylab='residual')
plot(res$Y1_obs,res$Y1_stdres,xlab='C_Jons',ylab='residual')
if(missingSource){
  col=nS*nEvt+nElt+(1:2)
} else {
  col=(nS-1)*nEvt+(1:2)
}
boxplot(mcmc[,col],names=c('gamma0','gamma1'))
dev.off()


# profiles
if(missingSource){np=nS+1} else {np=nS}
pdf(paste(folder,prefix,'profiles.pdf',sep=''),width=12,height=3)
par(mfrow=c(1,np))
for(i in 1:nS){
  plot(res[1:nElt,2+i],type='b',pch=19,
       xlab='Elt',ylab='C',ylim=c(0.5,1.7),main=nam[i])
}

if(missingSource){
  col=nS*nEvt+(1:nElt)
  boxplot(mcmc[,col],names=1:nElt,main=nam[nS+1],xlab='Elt',ylab='C')
}

dev.off()
