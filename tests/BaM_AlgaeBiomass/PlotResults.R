# MCMC
mcmc=read.table('Results_MCMC.txt',header=T)
p=NCOL(mcmc);n=NROW(mcmc)
X11(width=12);par(mfrow=c(ceiling(2*p/8),8))
for(i in 1:p){
  hist(mcmc[ceiling(n/2):n,i],xlab=names(mcmc)[i],main=names(mcmc)[i])
  plot(mcmc[,i],type='l',ylab=names(mcmc)[i])
}

# Predictions
files=c('PredictedBiomass','PredictedFb','PredictedFt',
        'PredictedFi','PredictedFn','PredictedRhyd')
n=length(files)
X11(width=12);par(mfrow=c(ceiling(n/3),3))
for(i in 1:n){
  f=paste0(files[i],'.env')
  pred=read.table(f,header=T)
  plot(pred[,1],type='l',main=files[i],
       ylim=c(min(pred[,2]),max(pred[,3])),
       xlab='time step',ylab=files[i])
  lines(pred[,2],lty=2,col='red')
  lines(pred[,3],lty=2,col='red')
}
