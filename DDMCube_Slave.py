#
#  DDMCube_Slave.py
#  DDMCubeTeraGrid
#
#  Created by nicain on 4/29/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#

################################################################################
# Preamble
################################################################################

# Import packages:
import DDMCube, pickle
import analysisTools as at
import pbsTools as pt

# Grab the name of the settings file, for this run:
execfile('nameOfCurrentSettingsFile.dat')

# Load up the settings:
fIn = open(nameOfCurrentSettingsFile,'r')
(settings, FD, numberOfJobs, gitVersion) = pickle.load(fIn)
fIn.close()

# Run the sims:
(resultsArray, crossTimesArray) = DDMCube.DDMOU(settings, FD, numberOfJobs[0])

# Save output:
pt.saveToFile({'resultsArray':resultsArray,'crossTimesArray':crossTimesArray})



