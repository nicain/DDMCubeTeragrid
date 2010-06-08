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
import numpy, copy, os, pickle, math, sys
from math import sqrt as sqrt
import pylab as pl

# Main Settings:
DDMCubeDir = '~/Desktop/currentProjects/DDMCubeTeragrid'
saveResultDir = 'NeuroeconomicOptimizeBackupRunData'
outputDir = os.path.abspath(os.path.join(os.curdir,saveResultDir))
settingsFile = 'NeuroeconomicOptimizeBackup-Booboo_30406946-ab98-4bf9-a313-e209ff2e3547.settings'
saveFig = 1
smoothBeta = 100

# Plot Settings:
yBegin = 0
chopHatVals = np.array([0,1,1.5])
betaPlotIndex = 2

# Figure Saving Settings:
saveFileName = 'tradeoffNegMarginalize.eps'
figOutputDir = 'Default' # Change from "Default" to path if you want to specify.
saveType = 'eps'		 #   Default is parent directory.
transparentBackground = True

# Font settings:
fig_width = 3				# width in inches
fontSize = 12
labelSize = 10

# Set up:
currentDir = os.path.abspath(os.getcwd())
DDMCubeDir = os.path.abspath(os.path.expanduser(DDMCubeDir))
if figOutputDir == 'Default':
	figOutputDir = os.path.abspath(os.pardir)
os.chdir(DDMCubeDir)
sys.path.append(os.getcwd())
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
		  'axes.labelsize': fontSize,
		  'text.fontsize': fontSize,
		  'legend.fontsize': fontSize,
		  'xtick.labelsize': labelSize,
		  'ytick.labelsize': labelSize,
		  'text.usetex': True,
		  'figure.figsize': fig_size}
pl.rcParams.update(params)

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

# Define the smoothing function:
def smooth(y,smoothBeta=smoothBeta):
	m = len(y)
	
	p = pl.diff(pl.eye(m),3).transpose()
	A = pl.eye(m)+smoothBeta*pl.dot(p.transpose(),p)
	
	smoothY = pl.solve(A, y)
	
	return smoothY

# Set up plot window:
pl.figure(1)
pl.clf()
LCushion = 0.19
RCushion = .06
BCushion = .21
TCushion = .05
pl.axes([LCushion, BCushion, 1-RCushion-LCushion, 1-TCushion-BCushion])

# Initializations:
sliceDict={'yBegin':yBegin}
FC = [0]*len(chopHatVals)
RT = [0]*len(chopHatVals)

# Get data slices:
for i in range(len(chopHatVals)):
	sliceDict['chopHat'] = chopHatVals[i]
	crossTimeDataThisSlice = finalMarginalCrossTimeCube[betaPlotIndex]
	resultDataThisSlice = finalMarginalResultCube[betaPlotIndex]
	crossDimsTemp = copy.copy(finalDimsMarginal)
	for collapseDim in iter(sliceDict):
		crossTimeDataThisSlice, resultDataThisSlice, crossDimsTemp = at.reduce1D(crossTimeDataThisSlice, resultDataThisSlice, crossDimsTemp, collapseDim, settings[collapseDim], sliceDict[collapseDim])
	RT[i] = np.squeeze(crossTimeDataThisSlice)
	FC[i] = np.squeeze(resultDataThisSlice)

# Plot Data:
for i in range(len(chopHatVals)):
	pl.plot(RT[i], FC[i])

#for i in range(len(chopHatVals)):
#	pl.plot(smooth(RT[i]), smooth(FC[i]))

# Set Axes:
pl.xlabel('RT')
pl.ylabel('FC')


# Save figure:
if saveFig == 1:
	(pl.savefig(os.path.join(figOutputDir, saveFileName),
		format = saveType))#,transparent = transparentBackground))

################################################################################
# Clean Up
################################################################################

# Change back directory:
os.chdir(currentDir)










