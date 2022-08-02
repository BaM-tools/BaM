import os
import sys
import smash
import numpy
import rasterio
import matplotlib.pyplot as plt


projectDir=sys.argv[1]
os.chdir(projectDir)

# load model from pickle (avoids re-reading inputs etc.)
model=smash.read_pickle("model.pkl")

# Define parameter
# uniform 
# model.parameters.cp=180. 
# distributed
myCP = rasterio.open("myCP.asc").read(1).transpose()
model.parameters.cp=myCP.copy()

# Run
model2=model.direct_run()
model2.plot_hydrograph()
plt.show()

