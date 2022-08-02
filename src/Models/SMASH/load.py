import os
import sys
import smash

projectDir=sys.argv[1]
precipDir=sys.argv[2]
petDir=sys.argv[3]
startDate=sys.argv[4]
endDate=sys.argv[5]

os.chdir(projectDir)
model=smash.Model(configuration="configuration.txt",dir_ph1=precipDir,dir_etp=petDir,date_deb=startDate,date_prv=endDate)
model.to_pickle("model.pkl")

