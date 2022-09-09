library(RBaM);library(gridExtra);library(ggplot2)
mcmc=readMCMC('Results_MCMC_Cooked.txt')
mcmc$cpMean=apply(mcmc[,1:5],1,mean)
res=read.table('Results_Residuals.txt',header=TRUE)
Qsim=read.table('Qsim_totalU.env',header=TRUE)
Qsim_param=read.table('Qsim_paramU.env',header=TRUE)
n=NROW(Qsim)

# MCMC plots
cp=tracePlot(mcmc[,c(paste0('cp',1:5),'cpMean')])
grid.arrange(grobs=cp)
violinPlot(mcmc[,c(paste0('cp',1:5),'cpMean')])

# Residuals
mask=res$Y1_obs>=0
g1=ggplot(res[mask,])+geom_point(aes(x=Y1_obs,y=Y1_sim))+theme_bw()
g2=ggplot(res[mask,])+geom_point(aes(x=Y1_obs,y=Y1_res))+theme_bw()
g3=ggplot(res[mask,])+geom_point(aes(x=Y1_obs,y=Y1_stdres))+theme_bw()
grid.arrange(g1,g2,g3,nrow=1)

# Qsim
g=ggplot()+
  geom_ribbon(data=Qsim,aes(x=1:n,ymin=Q_q2.5,ymax=Q_q97.5),fill='pink')+
  geom_ribbon(data=Qsim_param,aes(x=1:n,ymin=Q_q2.5,ymax=Q_q97.5),fill='red')+
  geom_line(data=res,aes(x=1:n,y=Y1_sim))+
  geom_point(data=res,aes(x=1:n,y=Y1_obs))+
  theme_bw()+
  coord_cartesian(ylim=c(0,200))
plotly::ggplotly(g)
