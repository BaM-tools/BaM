setwd('C:\\BEN\\FORGE\\BaRatin\\trunk\\test\\BaRatin_IVF\\Tests\\BaM_SFD')

# h1-h2 grid
h1=seq(1.5,6,0.1)
h2=seq(1,3,0.1)
grid=expand.grid(h1,h2)
write(grid[,1],'H1grid.txt',1)
write(grid[,2],'H2grid.txt',1)

# dummy h1-h2 limnis
nspag=500
h1=c(seq(2.1,6,0.01),seq(6,2.1,-0.005))
h2=rep(2,length(h1))
h1.spag=h1+matrix(rnorm(nspag*length(h1),0,0.1),length(h1),nspag)
h2.spag=h2+matrix(rnorm(nspag*length(h2),0,0.1),length(h2),nspag)
write(t(h1.spag),'H1limni.txt',nspag)
write(t(h2.spag),'H2limni.txt',nspag)
