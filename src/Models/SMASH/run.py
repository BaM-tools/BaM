import os
import sys
import smash
import numpy
import rasterio
import matplotlib.pyplot as plt

# Interpret command-line arguments
projectDir=sys.argv[1]
thetaDir=sys.argv[2]
# go into project directory
os.chdir(projectDir)
# load model from pickle (avoids re-reading inputs etc.)
model=smash.read_pickle("model.pkl")
# Define distributed parameters
cp = rasterio.open(os.path.join(thetaDir,"cp.asc")).read(1).transpose()
model.parameters.cp=cp.copy()
ctr = rasterio.open(os.path.join(thetaDir,"ctr.asc")).read(1).transpose()
model.parameters.ctr=ctr.copy()
cr = rasterio.open(os.path.join(thetaDir,"cr.asc")).read(1).transpose()
model.parameters.cr=cr.copy()
ml = rasterio.open(os.path.join(thetaDir,"ml.asc")).read(1).transpose()
model.parameters.ml=ml.copy()
# Run
model2=model.direct_run()
# Plot
model2.plot_hydrograph()
plt.show()

