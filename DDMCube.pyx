#  DDMCube.pyx
#  Created by nicain on 11/1/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.

# Wrapper for the RNG:
cdef extern from "MersenneTwister.h":
    ctypedef struct c_MTRand "MTRand":
        double randNorm( double mean, double stddev)
        void seed( unsigned long bigSeed[])

# External math functions that are needed:
cdef extern from "math.h":
    float sqrt(float sqrtMe)
    float fabs(float absMe)

################################################################################
######################## Main function, the workhorse:  ########################
################################################################################
def DDMOU(settings, int FD,int perLoc):

    # Import necessary python packages:
    import random, uuid, os, product
    from numpy import zeros
    from math import exp
    
    # C initializations
    cdef float xCurr, tCurr, yCurrP, yCurrN, C, xStd, xNoise, CPost, tieBreak, tBegin, CNull
    cdef float theta, crossTimes, results, chop, beta, betaBar, betaSigma, A, B,chopHat, noiseSigma, SNR
    cdef float dt = .1
    cdef float K = 9.0
    cdef float yTau = 20.0
    cdef float yBegin = 0.0
    cdef float xTau = 20.0
    cdef double mean = 0, std = 1
    cdef unsigned long mySeed[624]
    cdef c_MTRand myTwister
    cdef int i, overTime
    cdef float tMax = 10000
    cdef float deltaT
    
    # Convert settings dictionary to iterator:
    params = settings.keys()
    params.sort()
    settingsList = []
    totalLength = 1
    for parameter in params: 
        settingsList.append(settings[parameter])
        totalLength *= len(settings[parameter])
    settingsIterator = product.product(*settingsList)
    resultsArray = zeros(totalLength, dtype=float)
    crossTimesArray = zeros(totalLength, dtype=float)

    # Initialization of random number generator:
    myUUID = uuid.uuid4()
    random.seed(myUUID.int)
    for i in range(624): mySeed[i] = random.randint(0,int(exp(21)))
    myTwister.seed(mySeed)

    # Parameter space loop:
    counter = 0
    for currentSettings in settingsIterator:
        SNR, chopHat, deltaT, theta = currentSettings        # Must be alphabetized, with capitol letters coming first!

        results = 0
        crossTimes = 0
        for i in range(perLoc):
            tCurr = 0
            yCurrP = 0
            while fabs(yCurrP) < theta:

                xCurr = myTwister.randNorm(SNR,1.0)
                if fabs(xCurr) >= chopHat:
                    yCurrP += xCurr
                
                # Update Time Step
                tCurr=tCurr+deltaT


            if yCurrP > 0:
                results += 1
            elif yCurrP == 0:
                tieBreak = myTwister.randNorm(mean,std)
                if tieBreak > 0:
                    results += 1
        
            crossTimes += tCurr
    
        resultsArray[counter] = results
        crossTimesArray[counter] = crossTimes
        counter += 1
        
##        chop = sqrt(xStd*xStd + noiseSigma*noiseSigma)*chopHat
#        tBegin = 0
#        if FD:
#            theta = 1000000000
#        crossTimes = 0
#        results = 0
#        for i in range(perLoc):
#            beta = myTwister.randNorm(betaBar, betaSigma)
#            xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
#            overTime = 0
#            tCurr = 0
#            xCurr = myTwister.randNorm(C*.6,xStd)
#            xNoise = myTwister.randNorm(0,noiseSigma)
#            yCurrP = yBegin
#            yCurrN = yBegin
#            while yCurrP - yBegin < theta and yCurrN - yBegin < theta:
#                
#                # Create Input Signal
#                xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
#                xCurr = xCurr+dt*(C*.6 - xCurr)*xTau**(-1) + xStd*sqrt(2*dt*xTau**(-1))*myTwister.randNorm(mean,std)
#                
#                # Create Noise Signals
#                xNoise = xNoise - dt*xNoise*xTau**(-1) + noiseSigma*sqrt(2*dt*xTau**(-1))*myTwister.randNorm(mean,std)
#                
#                # Integrate Preferred Integrator based on chop
#                if fabs((xCurr+xNoise) + beta*yCurrP*K + 0) < chopHat*sqrt(xStd*xStd + noiseSigma*noiseSigma):
#                    yCurrP = yCurrP + 0
#                else:
#                    yCurrP = yCurrP + dt*yTau**(-1)*((xCurr+xNoise)*K**(-1) + beta*yCurrP + 0)
#
#                # Integrate Preferred Integrator based on chop                
#                if fabs(-(xCurr+xNoise) + beta*yCurrN*K + 0) < chopHat*sqrt(xStd*xStd + noiseSigma*noiseSigma):
#                    yCurrN = yCurrN + 0
#                else:
#                    yCurrN = yCurrN + dt*yTau**(-1)*(-(xCurr+xNoise)*K**(-1) + beta*yCurrN + 0)
#                
#                # Ensure both trains remain positive
#                if yCurrP < 0:
#                    yCurrP = 0
#                if yCurrN < 0:
#                    yCurrN = 0
#                
#                # Update Time Step
#                tCurr=tCurr+dt
#                
#                # Check if we are over time
#                if tCurr > tMax:
#                    overTime = 1
#                    break
#
#            crossTimes += tCurr
#            if FD:
#                if yCurrP > yCurrN:
#                    results += 1
#                elif yCurrP == yCurrN:
#                    tieBreak = myTwister.randNorm(mean,std)
#                    if tieBreak > 0:
#                        results += 1
#            else:
#                if not(overTime):
#                    if (yCurrP - yBegin >= theta) and (yCurrN - yBegin < theta):
#                        results += 1
#                    elif yCurrP == yCurrN:
#                        tieBreak = myTwister.randNorm(mean,std)
#                        if tieBreak > 0:
#                            results += 1
#                else:
#                    if yCurrP > yCurrN:
#                        results += 1
#                    elif yCurrP == yCurrN:
#                        tieBreak = myTwister.randNorm(mean,std)
#                        if tieBreak > 0:
#                            results += 1
                    
                    




    return (resultsArray, crossTimesArray)