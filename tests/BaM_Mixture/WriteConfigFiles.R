nS=5
nEvt=27
nElt=13
missingSource=T

writeConfig<-function(txt,file){
  write.table(matrix(txt, ncol = 1, nrow = length(txt)), file = file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Write Config_Options
txt=c(
  paste(nS,"! number of sources"),
  paste(nEvt,"! number of events"),  
  paste(nElt,"! number of elements")
  )
if(missingSource){foo=".true."} else {foo=".false."}
foo=paste(foo,"! add missing source?")
txt=c(txt,foo)
writeConfig(txt,"Config_Options.txt")

# Write Config_Model
if(missingSource){npar=nS*nEvt+nElt} else {npar=(nS-1)*nEvt}
txt=c(
  paste('"Mixture"','! model ID'),
  paste(nS+2,"! number of input variables"),  
  paste(1,"! number of output variables"),
  paste(npar,"! number of parameters")
)
for (i in 1:nEvt){
  for (j in 1:(nS-1)){
    foo=c(
      paste(paste("W",j,"_evt",i,sep=""),'! parameter name'),
      paste(1/(nS+1),'! initial guess'),
      paste('"Uniform"','! prior distribution'),
      paste('0,1','! prior parameters')
      )
    txt=c(txt,foo)
  }
  if(missingSource){
    foo=c(
      paste(paste("W",nS,"_evt",i,sep=""),'! parameter name'),
      paste(1/(nS+1),'! initial guess'),
      paste('"Uniform"','! prior distribution'),
      paste('0,1','! prior parameters')
    )
    txt=c(txt,foo)
  }
}
if(missingSource){
  for (i in 1:nElt){
    foo=c(
      paste(paste("Missing",i,sep=""),'! parameter name'),
      paste(1,'! initial guess'),
      paste('"FlatPrior+"','! prior distribution'),
      paste(' ','! prior parameters')
    )
    txt=c(txt,foo)    
  }  
}

writeConfig(txt,"Config_Model.txt")



