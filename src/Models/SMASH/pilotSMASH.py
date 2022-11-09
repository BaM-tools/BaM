import os
import sys
import smash
import numpy
import rasterio

# Make sense of command-line arguments
projectDir=sys.argv[1]
precipDir=sys.argv[2]
petDir=sys.argv[3]
commDir=sys.argv[4]

# Initialize SMASH
os.chdir(projectDir)
model=smash.Model(configuration="configuration.txt",dir_ph1=precipDir,dir_etp=petDir)
model.to_pickle("model.pkl")

# let BaM know initialization is done
open(os.path.join(commDir,'loadingDone.smash'),'a').close()

# start infinite loop running SMASH when BaM asks for it
keepGoing=True
while keepGoing:
	if os.path.exists(os.path.join(commDir,'run.bam')): # BaM is asking for a run
		# Define distributed parameters
		cp = rasterio.open(os.path.join(commDir,"cp.asc")).read(1).transpose()
		model.parameters.cp=cp.copy()
		ctr = rasterio.open(os.path.join(commDir,"ctr.asc")).read(1).transpose()
		model.parameters.ctr=ctr.copy()
		cr = rasterio.open(os.path.join(commDir,"cr.asc")).read(1).transpose()
		model.parameters.cr=cr.copy()
		ml = rasterio.open(os.path.join(commDir,"ml.asc")).read(1).transpose()
		model.parameters.ml=ml.copy()
		model2=model.direct_run() 
		# print('OK I ran, deleting trigger file...')
		os.remove(os.path.join(commDir,'run.bam')) # cleanup
		open(os.path.join(commDir,'runDone.smash'),'a').close() # tell BaM it's done
	if os.path.exists(os.path.join(commDir,'stop.bam')): # BaM says it's all done
		keepGoing=False
		print('BaM told me to stop')
		os.remove(os.path.join(commDir,'stop.bam'))  # cleanup
