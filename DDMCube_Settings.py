#!/usr/bin/env python

# Define job:
quickNamePrefix = 'OneDDebug'
dryRun = 0
localRun = 1
runType = 'batch'# 'batch' or 'wallTimeEstimate'
waitForSims = 1
wallTime = 30000
wallTimeEstCount = 3
queue = 'default'
FD=0

# Divide jobs among processing unit settings:
nodes = 1
procsPerNode = 1
repsPerProc = 1
simsPerRep = 10000

# Define job settings, master:
settings={# Example values:
    'SNR':[0.3901],#list(numpy.linspace(0,.4,1)),# 0
    'deltaT':[38],
    'chopHat':[0,1],#list(numpy.linspace(0,1.5,7)),# 0
    'tMax':list(numpy.linspace(500,590,1))}#list(numpy.linspace(15.0,15.8,2))}# 10

# Define job settings, master:
#settings={# Example values:
#    'betaSigma':[0],#list(numpy.linspace(0,.4,1)),# 0
#    'betaBar':[0],
#    'chopHat':[0],#list(numpy.linspace(0,1.5,7)),# 0
#    'tMax':list(numpy.linspace(10000,590,1)),# 2000, or 400->600 in FD paradigm
#    'thetaBegin':list(numpy.linspace(10,13.0,1)),# 0
#    'thetaHalf':list(numpy.linspace(500,13.0,1)),# 0
#    'thetaSS':list(numpy.linspace(30,13.0,1)),# 0
#    'C':[12.8],#(3.2*numpy.concatenate([[0],2**numpy.linspace(0,4,5)])),
#    'noiseSigma':[11.7]}#list(numpy.linspace(15.0,15.8,2))}# 10

#45 200
#25 200

#settings={# Example values:
#    'betaSigma':[.1],#list(numpy.linspace(0,.4,1)),# 0
#    'betaBar':[0],
#    'chopHat':[1.25],#list(numpy.linspace(0,1.5,7)),# 0
#    'tMax':list(numpy.linspace(10000,590,1)),# 2000, or 400->600 in FD paradigm
#    'theta':list(numpy.linspace(8,13.0,1)),# 0
#    'C':[12.8],#(3.2*numpy.concatenate([[0],2**numpy.linspace(0,4,5)])),
#    'noiseSigma':[17.2]}#list(numpy.linspace(15.0,15.8,2))}# 10

