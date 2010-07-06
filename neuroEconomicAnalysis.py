#!/usr/bin/env python

# neuroEconomicAnalysis.py
# DDMCubeTeraGrid
#
# Created by Nicholas Cain on 6/3/10.
# Copyright 2010 __MyCompanyName__. All rights reserved.

# Need the data files, the settings file from the save dir, and 
#	the original DDMCube_Settings.py file.

################################################################################
# Preamble
################################################################################

# Import necessary packages:
import pbsTools as pt
import analysisTools as at
import numpy as np
import numpy, copy, os, pickle, math, sys
from math import sqrt as sqrt
import pylab as pl
import sys

# Main Settings:
DDMCubeDir = '~/Desktop/currentProjects/DDMCubeTeragrid'
saveResultDir = 'NeuroEconomicOptimizeBig'
outputDir = os.path.abspath(os.path.join(os.curdir,saveResultDir))
quickName = 'NeuroeconomicOptimizeBig-Steele'

################################################################################
# Settings:
################################################################################

# Plot Settings:
yBegin = 0
chopHatVals = np.array([0,1.5,2])
betaStd = 0.05

#runDists = {
#    'beta':[('Normal',betaStd)],
#    'C':[('Delta',4)]}
    
#runDists = {
#    'beta':[('Delta',6),('Delta',7),('Delta',8),('Delta',9),('Delta',10)],
#    'C':[('Delta',4)]}

#runDists = {
#    'beta':[('Normal',betaStd)],
#    'C':map(lambda x: ('Delta',x), range(11))}
    
runDists = {
    'beta':map(lambda x: ('Delta',x), range(17)),
    'C':[('Delta',4)]}

# Figure Saving Settings:
saveFileNamePrefix = 'BetaDelta'
figOutputSubDir = 'fixedC_4E00'
saveFig = 1

#-------------------------------------------------------------------------------

# General Settings:
smoothBeta = 5

# Figure Saving Settings:
figOutputDir = '/Users/Nick/Desktop/currentProjects/DDMCubeTeragrid/NeuroEcoFigs'
saveType = 'eps'		 #   Default is parent directory.
transparentBackground = True

# Font settings:
fig_width = 3				# width in inches
fontSize = 12
labelSize = 10

################################################################################
# Set up:
################################################################################

# Set up:
currentDir = os.path.abspath(os.getcwd())
DDMCubeDir = os.path.abspath(os.path.expanduser(DDMCubeDir))
if figOutputDir == 'Parent':
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

################################################################################
# Define functions:
################################################################################

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
    
# Define the smoothing function:
def smooth(y,smoothBeta=smoothBeta):
    m = len(y)
    
    p = pl.diff(pl.eye(m),3).transpose()
    A = pl.eye(m)+smoothBeta*pl.dot(p.transpose(),p)
    
    smoothY = pl.solve(A, y)
    
    return smoothY

# Define Main function:
def makePlot(jobDist):

    ############################################################################
    # Do the marginalization:
    ############################################################################

    # Setup:
    saveName = saveFileNamePrefix

    # Import Data:
    crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion =  at.getDataAndSettings(quickName, saveResultDir)

    # Reduce Constant Dimensions:
    constDims = ['A', 'B', 'K', 'dt', 'noiseSigma', 'tMax', 'xStd', 'xTau', 'yTau']
    for collapseDim in constDims:
        crossTimeData, resultData, dims = at.reduce1D(crossTimeData, resultData, dims, collapseDim, settings[collapseDim], settings[collapseDim])
    crossTimeData = np.squeeze(crossTimeData)
    resultData = np.squeeze(resultData)

    # Pick Distribution, C:
    if jobDist['C'][0] == 'Uniform':
        numOfValues = len(settings['C'])
        probDist = np.ones(numOfValues)*(1./numOfValues)
        saveName += '_C_Uniform'
    elif jobDist['C'][0] == 'Delta':
        numOfValues = len(settings['C'])
        probDist = np.zeros(numOfValues)
        probDist[jobDist['C'][1]] = 1
        saveName += '_C_Delta_' + str(jobDist['C'][1])
    else:
        print 'Not Supported: ' + jobDist['C'][0]
        sys.exit()

    # Marginalize across C:
    marginalVal = settings['C'][0]
    crossTimeCubeMarginal, resultCubeMarginal, newDimsMarginal = at.reduce1D(crossTimeData, resultData, copy.copy(dims), 'C', settings['C'], marginalVal)
    crossTimeCubeMarginal *= probDist[0]
    resultCubeMarginal *= probDist[0]
    for i in range(1,numOfValues):
        marginalVal = settings['C'][i]
        crossTimeCubeTemp, resultCubeTemp, newDimsMarginalTemp = at.reduce1D(crossTimeData, resultData, copy.copy(dims), 'C', settings['C'], marginalVal)
        crossTimeCubeMarginal += probDist[i]*crossTimeCubeTemp
        resultCubeMarginal += probDist[i]*resultCubeTemp
    
    # Pick Distribution, beta:
    if jobDist['beta'][0] == 'Normal':
        numOfValues = len(settings['beta'])
        betaVals = np.array(settings['beta'])
        delta = betaVals[1]-betaVals[0]
        betaVals = betaVals - delta*1./2
        L = np.append(-10, betaVals[1:])
        R = np.append(betaVals[1:], 10)
        probDists = [0]*len(L)
        for j in range(len(L)):
            probDists[j] = intErfAB(L[j], R[j], sigma=jobDist['beta'][1])
        saveName += '_beta_Normal_' + str(jobDist['beta'][1])
    elif jobDist['beta'][0] == 'Delta':
        numOfValues = len(settings['beta'])
        probDists = np.zeros(numOfValues)
        probDists[jobDist['beta'][1]] = 1
        saveName += '_beta_Delta_' + str(jobDist['beta'][1])
    else:
        print 'Not Supported: ' + jobDist['beta'][0]
        sys.exit()

    # Marginalize across beta:
    marginalVal = settings['beta'][0]
    finalMarginalCrossTimeCube, finalMarginalResultCube, finalDimsMarginal = at.reduce1D(crossTimeCubeMarginal, resultCubeMarginal, copy.copy(newDimsMarginal), 'beta', settings['beta'], marginalVal)
    finalMarginalCrossTimeCube *= probDists[0]
    finalMarginalResultCube *= probDists[0]
    for j in range(1,numOfValues):
        marginalVal = settings['beta'][j]
        crossTimeCubeTemp, resultCubeTemp, newDimsMarginalTemp = at.reduce1D(crossTimeCubeMarginal, resultCubeMarginal, copy.copy(newDimsMarginal), 'beta', settings['beta'], marginalVal)
        finalMarginalCrossTimeCube += probDists[j]*crossTimeCubeTemp
        finalMarginalResultCube += probDists[j]*resultCubeTemp

    ############################################################################
    # Now make Speed/Accuracy Tradeoff Plot:
    ############################################################################

    # Initializations:
    sliceDict={'yBegin':yBegin}
    FC = [0]*len(chopHatVals)
    RT = [0]*len(chopHatVals)

    # Get data slices:
    for i in range(len(chopHatVals)):
        sliceDict['chopHat'] = chopHatVals[i]
        crossTimeDataThisSlice = finalMarginalCrossTimeCube
        resultDataThisSlice = finalMarginalResultCube
        crossDimsTemp = copy.copy(finalDimsMarginal)
        for collapseDim in iter(sliceDict):
            crossTimeDataThisSlice, resultDataThisSlice, crossDimsTemp = at.reduce1D(crossTimeDataThisSlice, resultDataThisSlice, crossDimsTemp, collapseDim, settings[collapseDim], sliceDict[collapseDim])
        RT[i] = np.squeeze(crossTimeDataThisSlice)
        FC[i] = np.squeeze(resultDataThisSlice)

    # Set up plot window:
    pl.figure(1)
    pl.clf()
    LCushion = 0.19
    RCushion = .06
    BCushion = .21
    TCushion = .05
    pl.axes([LCushion, BCushion, 1-RCushion-LCushion, 1-TCushion-BCushion])

    # Plot Data: Comment out one of the two, for smoothing or not...
    for i in range(len(chopHatVals)):
        if i == 0:
            pl.plot((RT[i]), (FC[i]),'--')
        else:
            pl.plot((RT[i]), (FC[i]))
        
    # Set Axes:
    pl.xlabel('RT')
    pl.ylabel('FC')
    pl.ylim([.5,1])

    # Save figure:
    if saveFig == 1:
        (pl.savefig(os.path.join(figOutputDir, figOutputSubDir, saveName + '.' + saveType),
            format = saveType, transparent = transparentBackground))
    
    return (RT,FC)

################################################################################
# Generate the plots:
################################################################################

# Create jobDist from runDist

# Call the Main function:
for betaSetting in runDists['beta']:
    for CSetting in runDists['C']:
        makePlot({'beta':betaSetting,'C':CSetting})

################################################################################
# Clean Up
################################################################################

# Change back directory:
os.chdir(currentDir)

















################################################################
# Gather data:
################################################################

#settingsFile = 'NeuroeconomicOptimizeBig-Steele_a6f9c67c-fc7e-4bcc-b687-3ebb002d0104.settings'

# Split off uuid fomr settingFile:
#jobQuickNameIn, currentSuffix = settingsFile.split('_')
#myUUID, trash = currentSuffix.split('.')
#quickNamePrefix, quickNameSuffix = jobQuickNameIn.split('-')
#
# Import settings:
#execfile(os.path.join(outputDir,'DDMCube_Settings.py'))
#quickName = quickNamePrefix + '-' + quickNameSuffix
#numberOfJobs = [simsPerRep, simsPerRep*repsPerProc*procsPerNode*nodes]
#
# Grab the name of the settings for the run:
#settings, FD, numberOfJobs, gitVersion = pt.unpickle(os.path.join(os.getcwd(),saveResultDir,settingsFile))
#
# Generate crossTimesArray and resultsArray:
#resultList = pt.getFromPickleJar(loadDir = outputDir, fileNameSubString = 'simResults.dat')
#arrayLength = len(resultList[0][0])
#resultsArray = np.zeros(arrayLength, dtype=float)
#crossTimesArray = np.zeros(arrayLength, dtype=float)
#for i in range(len(resultList)):
#	resultsArray = resultsArray + resultList[i][0]
#	crossTimesArray = crossTimesArray + resultList[i][1]
#crossTimesArray = crossTimesArray/numberOfJobs[1]
#resultsArray = resultsArray/numberOfJobs[1]
#
# Reshape results:	
#params = settings.keys()
#params.sort()
#newDims = [len(settings[parameter]) for parameter in params]
#crossTimesArray = np.reshape(crossTimesArray,newDims)
#resultsArray = np.reshape(resultsArray,newDims)
#
# Export data in familiar format, to outputDir:
#fOut = open(os.path.join(outputDir,quickName + '_' + str(myUUID) + '.dat'),'w')
#pickle.dump((crossTimesArray, resultsArray, params),fOut)
#fOut.close()






