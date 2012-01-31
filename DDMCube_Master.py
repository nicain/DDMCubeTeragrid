#       
#  DDMCubeTeraGrid_Master.py
#  DDMCubeTeraGrid
#
#  Created by nicain on 4/28/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#

################################################################################
# Preamble
################################################################################

# Import packages:
from subprocess import call as call
from os.path import join as join, abspath as abspath
import numpy, random, time, os, pickle, uuid
import analysisTools as at
import pbsTools as pt

# Compile cython extension DDMCube.pyx:
call('python setup.py build_ext --inplace', shell=True)

################################################################################
# Main function: 
################################################################################

# Import settings:
execfile('DDMCube_Settings.py')
outputDir = '.batchSimResults'
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
elif quickNameSuffix == 'AMath':
	user='nicain'
	runLocation = 'cluster'
elif quickNameSuffix == 'Hyak':
	user='nicain'
	runLocation = 'hyak'

if quickNameSuffix != 'Hyak':
	pythonPath = ''
else:
	pythonPath = '/usr/lusers/nicain/epd-7.0-2-rh5-x86_64/bin/'

# Beginning computation:
quickName = quickNamePrefix + '-' + quickNameSuffix
numberOfJobs = [simsPerRep, simsPerRep*repsPerProc*procsPerNode*nodes]

# Write a "settings" file:
myUUID = uuid.uuid4()
gitVersion = 'None'
totalLength = 1
for parameter in settings: 
	thisSetting = settings[parameter]
	totalLength *= len(thisSetting)
settingsFileName = join(os.getcwd(), saveResultDir, quickName + '_' + str(myUUID) + '.settings')
pt.pickle((settings, FD, numberOfJobs, gitVersion), saveFileName = settingsFileName)

# Display settings:
at.printSettings(quickName, saveResultDir = saveResultDir)

# Run the job:
tBegin = time.mktime(time.localtime())
pt.runPBS(pythonPath + 'python DDMCube_Slave.py',
          fileList = ['DDMCube_Slave.py', settingsFileName, 'DDMCube.so'],
          nodes=nodes,
          ppn=procsPerNode,
		  repspp=repsPerProc,
		  outputDir=outputDir,
          runLocation=runLocation,
		  runType=runType,
		  waitForSims=waitForSims,
          wallTime=wallTime,
		  dryRun=dryRun,
		  queue=queue,
		  wallTimeEstCount=wallTimeEstCount)
tEnd = time.mktime(time.localtime())

if not dryRun == 1 and not runType == 'wallTimeEstimate' and not waitForSims == 0:
	# Collect results:
	resultList = pt.getFromPickleJar(loadDir = outputDir, fileNameSubString = 'simResults.dat')
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

# Display Computation Time:
print 'Total Computation Time: ', time.strftime("H:%H M:%M S:%S",time.gmtime(tEnd - tBegin))
if simsPerRep < 1000:
	for NN in [2000,5000]: print ' Time to complete ' + str(NN) +  ' sims: ', time.strftime("H:%H M:%M S:%S",time.gmtime(NN*totalLength*(tEnd - tBegin)/(totalLength*simsPerRep)))
