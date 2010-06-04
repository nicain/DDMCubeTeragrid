#!/usr/bin/env python

# neuroEconomicAnalysis.py
# DDMCubeTeraGrid
#
# Created by Nicholas Cain on 6/3/10.
# Copyright 2010 __MyCompanyName__. All rights reserved.

# Need the data files, the settings file from the save dir, and 
#	the original DDMCube_Settings.py file.

################################################################
# Preamble
################################################################

# Import necessary packages:
import pbsTools as pt
import analysisTools as at
import numpy as np
import numpy, copy, os, pickle, math

# Settings:
saveResultDir = 'NeuroeconomicOptimizeBackupRunData'
outputDir = os.path.abspath(os.path.join(os.curdir,saveResultDir))
settingsFile = 'NeuroeconomicOptimizeBackup-Booboo_30406946-ab98-4bf9-a313-e209ff2e3547.settings'

################################################################
# Gather data:
################################################################

# Split off uuid fomr settingFile:
jobQuickNameIn, currentSuffix = settingsFile.split('_')
myUUID, trash = currentSuffix.split('.')

# Import settings:
execfile(os.path.join(outputDir,'DDMCube_Settings.py'))
quickNameSuffix = os.environ['JOBLOCATION']
quickName = quickNamePrefix + '-' + quickNameSuffix
numberOfJobs = [simsPerRep, simsPerRep*repsPerProc*procsPerNode*nodes]

# Match incoming quickName against that from the settings file:
if not jobQuickNameIn == quickName:
	print 'Quicknames from input file and settings file dont match; aborting'
	sys.exit()

# Grab the name of the settings for the run:
settings, FD, numberOfJobs, gitVersion = pt.unpickle(os.path.join(os.getcwd(),saveResultDir,settingsFile))

# Generate crossTimesArray and resultsArray:
resultList = pt.getFromPickleJar(loadDir = outputDir, fileNameSubString = 'simResults.dat')
arrayLength = len(resultList[0][0])
resultsArray = np.zeros(arrayLength, dtype=float)
crossTimesArray = np.zeros(arrayLength, dtype=float)
for i in range(len(resultList)):
	resultsArray = resultsArray + resultList[i][0]
	crossTimesArray = crossTimesArray + resultList[i][1]
crossTimesArray = crossTimesArray/numberOfJobs[1]
resultsArray = resultsArray/numberOfJobs[1]

# Reshape results:	
params = settings.keys()
params.sort()
newDims = [len(settings[parameter]) for parameter in params]
crossTimesArray = np.reshape(crossTimesArray,newDims)
resultsArray = np.reshape(resultsArray,newDims)

# Export data in familiar format, to outputDir:
fOut = open(os.path.join(outputDir,quickName + '_' + str(myUUID) + '.dat'),'w')
pickle.dump((crossTimesArray, resultsArray, params),fOut)
fOut.close()

################################################################
# Do the marginalization:
################################################################

# Import Data:
crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion =  at.getDataAndSettings(quickName, saveResultDir)

# Reduce Constant Dimensions:
constDims = ['A', 'B', 'K', 'dt', 'noiseSigma', 'tMax', 'xStd', 'xTau', 'yTau']
for collapseDim in constDims:
	crossTimeData, resultData, dims = at.reduce1D(crossTimeData, resultData, dims, collapseDim, settings[collapseDim], settings[collapseDim])
crossTimeData = np.squeeze(crossTimeData)
resultData = np.squeeze(resultData)

# Marginalize across C:
marginalizeDim = 'C'
numOfValues = len(settings[marginalizeDim])

probDist = np.ones(numOfValues)*(1./numOfValues)
marginalVal = settings[marginalizeDim][0]
crossTimeCubeMarginal, resultCubeMarginal, newDimsMarginal = at.reduce1D(crossTimeData, resultData, copy.copy(dims), marginalizeDim, settings[marginalizeDim], marginalVal)
crossTimeCubeMarginal *= probDist[0]
resultCubeMarginal *= probDist[0]
for i in range(1,numOfValues):
	marginalVal = settings[marginalizeDim][i]
	crossTimeCubeTemp, resultCubeTemp, newDimsMarginalTemp = at.reduce1D(crossTimeData, resultData, copy.copy(dims), marginalizeDim, settings[marginalizeDim], marginalVal)
	crossTimeCubeMarginal += probDist[i]*crossTimeCubeTemp
	resultCubeMarginal += probDist[i]*resultCubeTemp
	
# Define the functions for the gaussian integral:
def intErfAB(a, b, mu=0, sigma=1):
	
	# Ill need erf(x)
	def erf(x):
		# save the sign of x
		sign = 1
		if x < 0: 
			sign = -1
		x = abs(x)

		# constants
		a1 =  0.254829592
		a2 = -0.284496736
		a3 =  1.421413741
		a4 = -1.453152027
		a5 =  1.061405429
		p  =  0.3275911

		# A&S formula 7.1.26
		t = 1.0/(1.0 + p*x)
		y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
		return sign*y
	
	a=a*1.
	b=b*1.
	mu=mu*1.
	sigma=sigma*1.
	C1 = (b-mu)/(math.sqrt(2)*sigma)
	C2 = (a-mu)/(math.sqrt(2)*sigma)
	P1 = erf(C1)
	P2 = erf(C2)
	
	return .5*(P1-P2)

# Marginalize across beta:
marginalizeDim = 'beta'
numOfValues = len(settings[marginalizeDim])

print numOfValues

betaStd = [.01,.02,.05]
probDists = [0]*len(betaStd)
betaVals = np.array(settings[marginalizeDim])
delta = betaVals[1]-betaVals[0]
betaVals = betaVals - delta*1./2
L = np.append(-10, betaVals[1:])
R = np.append(betaVals[1:], 10)
for i in range(len(probDists)):
	probDists[i] = [0]*len(L)
	for j in range(len(L)):
		probDists[i][j] = intErfAB(L[j], R[j], sigma=betaStd[i])
		
print len(probDists[0])
print len(L)
print len(R)		
print len(betaVals)

finalMarginalCrossTimeCube = [0]*len(betaStd)
finalMarginalResultCube = [0]*len(betaStd)
for i in range(len(betaStd)):
	marginalVal = settings[marginalizeDim][0]
	finalMarginalCrossTimeCube[i], finalMarginalResultCube[i], finalDimsMarginal = at.reduce1D(crossTimeCubeMarginal, resultCubeMarginal, copy.copy(newDimsMarginal), marginalizeDim, settings[marginalizeDim], marginalVal)
	finalMarginalCrossTimeCube[i] *= probDists[i][0]
	finalMarginalResultCube[i] *= probDists[i][0]
	for j in range(1,numOfValues):
		marginalVal = settings[marginalizeDim][j]
		crossTimeCubeTemp, resultCubeTemp, newDimsMarginalTemp = at.reduce1D(crossTimeCubeMarginal, resultCubeMarginal, copy.copy(newDimsMarginal), marginalizeDim, settings[marginalizeDim], marginalVal)
		finalMarginalCrossTimeCube[i] += probDists[i][j]*crossTimeCubeTemp
		finalMarginalResultCube[i] += probDists[i][j]*resultCubeTemp

################################################################
# Now make Speed/Accuracy Tradeoff Plot:
################################################################






















