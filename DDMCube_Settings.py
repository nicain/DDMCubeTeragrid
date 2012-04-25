#!/usr/bin/env python

# Define job:
quickNamePrefix = 'DiscreteFDDebug'
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
simsPerRep = 500

# Define job settings, master:
settings={# Example values:
    'SNR':[0.3901],#list(numpy.linspace(0,.4,1)),# 0
    'deltaT':[38],
    'chopHat':[0,.5,1],#list(numpy.linspace(0,1.5,7)),# 0
    'theta':list(numpy.linspace(1,20,11))}#list(numpy.linspace(15.0,15.8,2))}# 10


