#!/usr/bin/env python

# Define job:
quickNamePrefix = 'NeuralDebug'
dryRun = 0
localRun = 1
runType = 'batch'# 'batch' or 'wallTimeEstimate'
waitForSims = 1
wallTime = 5000
wallTimeEstCount = 1
queue = 'default'
FD=0

# Divide jobs among processing unit settings:
nodes = 1
procsPerNode = 1
repsPerProc = 1
simsPerRep = 1000

# Define job settings, effective:
settings={# Example values:
    'C':[12.8],#(3.2*numpy.concatenate([[0],2**numpy.linspace(0,4,5)])),
    'K':[9],
    'betaSigma':[.1],
    'chopHat':[1.25],#list(numpy.linspace(0,1.5,7)),# 0    
    'dt':[.1],#list(numpy.linspace(.1,.1,1)),# .02
    'noiseSigma':[14.0],
    'numberOfNeurons':[250],
    'tMax':list(numpy.linspace(10000,590,1)),# 2000, or 400->600 in FD paradigm
    'theta':list(numpy.linspace(3,20,1))}# 0


# Define job settings, master:
#settings={# Example values:
#    'A':list(numpy.linspace(0,5,1)),# 0
#    'B':list(numpy.linspace(0,0,1)),# 0
#    'beta':[0],#list(numpy.linspace(0,.4,1)),# 0
#    'chopHat':[1.25],#list(numpy.linspace(0,1.5,7)),# 0
#    'dt':[.1],#list(numpy.linspace(.1,.1,1)),# .02
#    'K':list(numpy.linspace(9,10.5,1)),# 20
#    'tMax':list(numpy.linspace(10000,590,1)),# 2000, or 400->600 in FD paradigm
#    'tFrac':list(numpy.linspace(1,1,1)),
#	'tBeginFrac':list(numpy.linspace(0,1,1)),
#    'theta':list(numpy.linspace(7.6,20,1)),# 0
#    'COn':list(numpy.linspace(0.00,51.20,11)),#[12.8],#(3.2*numpy.concatenate([[0],2**numpy.linspace(0,4,5)])),
#    'xStd':list(numpy.linspace(0,15,1)),# 12.8; 0 means set in DDMCube as 4.46*Mean
#    'xTau':list(numpy.linspace(20,25,1)),# 20
#    'yBegin':list(numpy.linspace(0,40,1)),# 40
#    'yTau':list(numpy.linspace(20,10,1)),
#    'noiseSigma':[14.0]}#list(numpy.linspace(12.0,60,1))}# 10

# Define job settings, collapsingBound:
#settings={# Example values:
#    'A':list(numpy.linspace(0,5,1)),# 0
#    'B':list(numpy.linspace(0,0,1)),# 0
#    'betaMu':[0],#list(numpy.linspace(0,.4,1)),# 0
#	'betaSigma':[0,.1],#list(numpy.linspace(0,.4,1)),# 0
#    'chopHat':[0,1,1.5,2],#list(numpy.linspace(2,1,1)),# 0
#    'dt':[.1],#list(numpy.linspace(.1,.1,1)),# .02
#    'K':list(numpy.linspace(9,10.5,1)),# 20
#    'tMax':list(numpy.linspace(10000,1000,1)),# 2000, or 400->600 in FD paradigm
#    'tFrac':list(numpy.linspace(1,1,1)),
#    'thetaBegin':list((numpy.linspace(.1,200,200)**1.5)*(200**(1-1.5))),#[25],#list(numpy.linspace(1,50,7)),# 0
#	'thetaHalf':[500],#list(numpy.linspace(1,50,7)),# 0
#    'thetaSS':[0],#list(numpy.linspace(1,50,7)),# 0
#    'CPre':[12.8],#(3.2*numpy.concatenate([[0],2**numpy.linspace(0,4,5)])),
#    'xStd':list(numpy.linspace(0,15,1)),# 12.8; 0 means set in DDMCube as 4.46*Mean
#    'xTau':list(numpy.linspace(20,25,1)),# 20
#    'yBegin':list(numpy.linspace(0,40,1)),# 40
#    'yTau':list(numpy.linspace(20,10,1)),
#    'noiseSigma':[14.0]}#list(numpy.linspace(12.0,60,1))}# 10

# Define job settings, returnRTCollapsingBound:
#settings={# Example values:
#    'A':list(numpy.linspace(0,5,1)),# 0
#    'B':list(numpy.linspace(0,0,1)),# 0
#    'betaMu':[0],#list(numpy.linspace(0,.4,1)),# 0
#	'betaSigma':[.1],#list(numpy.linspace(0,.4,1)),# 0
#    'chopHat':[2],#list(numpy.linspace(2,1,1)),# 0
#    'dt':[.1],#list(numpy.linspace(.1,.1,1)),# .02
#    'K':list(numpy.linspace(9,10.5,1)),# 20
#    'tMax':list(numpy.linspace(10000,1000,1)),# 2000, or 400->600 in FD paradigm
#    'tFrac':list(numpy.linspace(1,1,1)),
#    'thetaBegin':[25],#list(numpy.linspace(1,50,7)),# 0
#	'thetaHalf':[140],#list(numpy.linspace(1,50,7)),# 0
#    'thetaSS':[0],#list(numpy.linspace(1,50,7)),# 0
#    'CPre':(3.2*numpy.concatenate([[0],2**numpy.linspace(0,4,5)])),
#    'xStd':list(numpy.linspace(0,15,1)),# 12.8; 0 means set in DDMCube as 4.46*Mean
#    'xTau':list(numpy.linspace(20,25,1)),# 20
#    'yBegin':list(numpy.linspace(0,40,1)),# 40
#    'yTau':list(numpy.linspace(20,10,1)),
#    'noiseSigma':[12.0]}#list(numpy.linspace(12.0,60,1))}# 10

# Define job settings, randomDelay:
#settings={# Example values:
#    'A':list(numpy.linspace(0,5,1)),# 0
#    'B':list(numpy.linspace(0,0,1)),# 0
#    'beta':list(numpy.linspace(0,.4,1)),# 0
#    'chopHat':[0,2],#list(numpy.linspace(2,1,1)),# 0
#    'dt':[.1],#list(numpy.linspace(.1,.1,1)),# .02
#    'K':list(numpy.linspace(9,10.5,1)),# 20
#    'tMax':list(numpy.linspace(10000,1000,1)),# 2000, or 400->600 in FD paradigm
#    'tFrac':list(numpy.linspace(1,1,1)),
#	'tBeginFrac':list(numpy.linspace(0,1,1)),
#    'theta':list(numpy.linspace(1,50,7)),# 0
#    'COn':[51.2],#(3.2*numpy.concatenate([[0],2**numpy.linspace(0,4,5)])),
#    'xStd':list(numpy.linspace(0,15,1)),# 12.8; 0 means set in DDMCube as 4.46*Mean
#    'xTau':list(numpy.linspace(20,25,1)),# 20
#    'yBegin':list(numpy.linspace(0,40,1)),# 40
#    'yTau':list(numpy.linspace(20,10,1)),
#    'noiseSigma':list(numpy.linspace(12.0,60,1)),
#	'delayMean':[2000]}# 10

# Define job settings, returnRT:
#settings={# Example values:
#    'A':list(numpy.linspace(0,5,1)),# 0
#    'B':list(numpy.linspace(0,0,1)),# 0
#    'betaMu':[0],#list(numpy.linspace(0,.4,1)),# 0
#	 'betaSigma':[.1],#list(numpy.linspace(0,.4,1)),# 0
#    'chopHat':[2],#list(numpy.linspace(2,1,1)),# 0
#    'dt':[.1],#list(numpy.linspace(.1,.1,1)),# .02
#    'K':list(numpy.linspace(9,10.5,1)),# 20
#    'tMax':list(numpy.linspace(10000,1000,1)),# 2000, or 400->600 in FD paradigm
#    'tFrac':list(numpy.linspace(1,1,1)),
#    'theta':[7.2],#list(numpy.linspace(1,50,7)),# 0
#    'CPre':(3.2*numpy.concatenate([[0],2**numpy.linspace(0,4,5)])),
#    'xStd':list(numpy.linspace(0,15,1)),# 12.8; 0 means set in DDMCube as 4.46*Mean
#    'xTau':list(numpy.linspace(20,25,1)),# 20
#    'yBegin':list(numpy.linspace(0,40,1)),# 40
#    'yTau':list(numpy.linspace(20,10,1)),
#    'noiseSigma':[12.0]}#list(numpy.linspace(12.0,60,1))}# 10