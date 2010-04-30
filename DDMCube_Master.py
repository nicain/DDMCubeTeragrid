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
import scipy, random, time, os, pickle, uuid
import analysisTools as at
import pbsTools as pt
DDMCube_cythonDevDir = '/Users/Nick/Library/Python/src/DDMCube_cython'

# Compile cython extension DDMCube.pyx:
call('python setup.py build_ext --inplace', shell=True)

################################################################################
# Dashboard:
################################################################################

# Define job settings:
settings={									# Example values:
'A':list(scipy.linspace(0,0,1)),			# 0
'B':list(scipy.linspace(0,0,1)),			# 0
'beta':list(scipy.linspace(0,.1,1)),		# 0
'chopHat':list(scipy.linspace(0,2,1)),			# 0
'dt':list(scipy.linspace(.1,.1,1)),		# .02
'K':list(scipy.linspace(9,10.5,1)),		# 20
'tMax':list(scipy.linspace(2000,1000,1)),	# 2000, or 400->600 in FD paradigm
'theta':list(scipy.linspace(1,100,41)),		# 10
'C':[6.4],#(3.2*scipy.concatenate([[0],2**scipy.linspace(0,4,5)])),  #[0,3.2,6.4,12.8,25.6,51.2],#list(scipy.linspace(0,20,7)),		# Now units of C
'xStd':list(scipy.linspace(0,15,1)),		# 12.8; 0 means set in DDMCube as 4.46*Mean
'xTau':list(scipy.linspace(20,25,1)),		# 20
'yBegin':list(scipy.linspace(0,40,1)),		# 40
'yTau':list(scipy.linspace(20,10,1)),
'noiseSigma':list(scipy.linspace(14,15,2))		# 10
}

# Define job parameters:
quickName = 'debug'
FD=0

# Define type of job
dryRun = 0
localRun = 0						# Superceded by dryRun
server = 'wallTimeEstimate'			# 'normal' or 'wallTimeEstimate'

# Divide jobs among processing unit settings:
nodes = 1
procsPerNode = 2
repsPerProc = 1
simsPerRep = 5
numberOfJobs = [simsPerRep, simsPerRep*repsPerProc*procsPerNode*nodes]


################################################################################
# Main function: 
################################################################################

# Set up saving directories
saveResultDir = 'savedResults'
outputDir = 'batchSimResults'

# Write a "settings" file:
myUUID = uuid.uuid4()
gitVersion = 'None'
totalLength = 1
for parameter in settings: 
	thisSetting = settings[parameter]
	totalLength *= len(thisSetting)
settingsFileName = join(os.getcwd(), saveResultDir, quickName + '_' + str(myUUID) + '.settings')
fOutSet = open(settingsFileName,'w')	
pickle.dump((settings, FD, numberOfJobs, gitVersion),fOutSet)
fOutSet.close()

# Write file containing the name of the settings file, for the subordinate jobs:
fOut = open('nameOfCurrentSettingsFile.dat','w')
fOut.write('nameOfCurrentSettingsFile = "' + os.path.abspath(settingsFileName) + '"')
fOut.close()

# Display settings:
at.printSettings(quickName, saveResultDir)

# Run the job:
tBegin = time.mktime(time.localtime())
pt.runPBS('python DDMCube_Slave.py',
          ['DDMCube_Slave.py', 'nameOfCurrentSettingsFile.dat', 'DDMCube.so'],
          nodes=nodes,
          ppn=procsPerNode,
		  repspp=repsPerProc,
		  outputDir=outputDir,
          server=server,
          wallTime=10,
		  dryRun=dryRun,
		  localRun=localRun)
tEnd = time.mktime(time.localtime())

# Collect results:
resultList = pt.getSavedVariables(['resultsArray','crossTimesArray'], outputDir = outputDir)
arrayLength = len(resultList[0]['resultsArray'])

resultsArray = scipy.zeros(arrayLength, dtype=float)
crossTimesArray = scipy.zeros(arrayLength, dtype=float)
for i in range(len(resultList)):
	resultsArray = resultsArray + resultList[i]['resultsArray']
	crossTimesArray = crossTimesArray + resultList[i]['crossTimesArray']

crossTimesArray = crossTimesArray/numberOfJobs[1]
resultsArray = resultsArray/numberOfJobs[1]
			
# Reshape results and save to output:	
params = settings.keys()
params.sort()
newDims = [len(settings[parameter]) for parameter in params]
crossTimesArray = scipy.reshape(crossTimesArray,newDims)
resultsArray = scipy.reshape(resultsArray,newDims)
fOut = open(join(os.getcwd(),saveResultDir,quickName + '_' + str(myUUID) + '.dat'),'w')
pickle.dump((crossTimesArray, resultsArray, params),fOut)
fOut.close()

if localRun == 1:
	# Display Computation Time:
	print 'Total Computation Time: ', time.strftime("H:%H M:%M S:%S",time.gmtime(tEnd - tBegin))
	if simsPerRep < 1000:
		for NN in [2000,5000]: print ' Time to complete ' + str(NN) +  ' sims: ', time.strftime("H:%H M:%M S:%S",time.gmtime(NN*totalLength*(tEnd - tBegin)/(totalLength*simsPerRep)))