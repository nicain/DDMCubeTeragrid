#
#  recoverLostJob.py
#  DDMCubeTeraGrid
#
#  Created by nicain on 5/5/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#
import pbsTools as pt
import sys, numpy, os, pickle
from os.path import join as join, abspath as abspath

# Get details from input:
outputDir = sys.argv[1]
settingsFile = sys.argv[2]

# Read back recovery details to user:
print 'Recovering job: `' + str(settingsFile)
print '  From directory: ' + str(outputDir)

# Split off uuid fomr settingFile:
jobQuickNameIn, currentSuffix = settingsFile.split('_')
myUUID, trash = currentSuffix.split('.')

# Import settings:
execfile('DDMCube_Settings.py')
quickNameSuffix = os.environ['JOBLOCATION']
saveResultDir = 'savedResults-' + quickNameSuffix
if quickNameSuffix == 'Booboo' or localRun == 1:
	user='nicain'
	runLocation = 'local'
elif quickNameSuffix == 'Abe':
	user='ncain'
	runLocation = 'abe'
elif quickNameSuffix == 'Steele':
	user='cainn'
	runLocation = 'steele'
quickName = quickNamePrefix + '-' + quickNameSuffix
numberOfJobs = [simsPerRep, simsPerRep*repsPerProc*procsPerNode*nodes]

# Match incoming quickName against that from the settings file:
if not jobQuickNameIn == quickName:
	print 'Quicknames from input file and settings file dont match; aborting'
	sys.exit()

# Grab the name of the settings for the run:
settings, FD, numberOfJobs, gitVersion = pt.unpickle(os.path.join(os.getcwd(),saveResultDir,settingsFile))

# Save the day:
resultList = pt.getFromPickleJar(loadDir = outputDir, fileNameSubString = 'simResults.dat')
print resultList
arrayLength = len(resultList[0][0])

resultsArray = numpy.zeros(arrayLength, dtype=float)
crossTimesArray = numpy.zeros(arrayLength, dtype=float)
for i in range(len(resultList)):
	resultsArray = resultsArray + resultList[i][0]
	crossTimesArray = crossTimesArray + resultList[i][1]

crossTimesArray = crossTimesArray/numberOfJobs[1]
resultsArray = resultsArray/numberOfJobs[1]
			
# Reshape results and save to output:	
params = settings.keys()
params.sort()
newDims = [len(settings[parameter]) for parameter in params]
crossTimesArray = numpy.reshape(crossTimesArray,newDims)
resultsArray = numpy.reshape(resultsArray,newDims)
fOut = open(join(os.getcwd(),saveResultDir,quickName + '_' + str(myUUID) + '.dat'),'w')
pickle.dump((crossTimesArray, resultsArray, params),fOut)
fOut.close()
