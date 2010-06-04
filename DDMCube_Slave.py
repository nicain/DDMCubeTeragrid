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
import DDMCube
import analysisTools as at
import pbsTools as pt

# Grab the name of the settings for the run:
settings, FD, numberOfJobs, gitVersion = pt.getFromPickleJar(loadDir = './', fileNameSubString = '.settings')[0]

# Run the sims:
(resultsArray, crossTimesArray) = DDMCube.DDMOU(settings, FD, numberOfJobs[0])

# Save output:
pt.pickle((resultsArray, crossTimesArray),saveFileName = 'simResults.dat')



