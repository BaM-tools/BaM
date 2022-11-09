import os
import sys
import smash

projectDir=sys.argv[1]
precipDir=sys.argv[2]
petDir=sys.argv[3]

os.chdir(projectDir)
model=smash.Model(configuration="configuration.txt",dir_ph1=precipDir,dir_etp=petDir)
model.to_pickle("model.pkl")

