#  DDMCube.pyx
#  Created by nicain on 11/1/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.

# Import floating point division from python 3.0
from __future__ import division

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
    import numpy as np
    
    # C initializations
    cdef float xCurr, tCurr, yCurrPSum, yCurrNSum, C, xStd, xTau, xNoise, tieBreak, tMax
    cdef float dt, theta, crossTimes, results, chop, beta, K, yTau,chopHat, noiseSigma
    cdef double mean = 0, std = 1
    cdef unsigned long mySeed[624]
    cdef c_MTRand myTwister
    cdef int i, overTime
    cdef float q = 1
    cdef float betaSigma
    cdef int numberOfNeurons
    cdef float kappa, p, a
    cdef float rMin = 0
    cdef float rMax = 50
    
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
    CNull = 0
    for currentSettings in settingsIterator:
        C, K, betaSigma, chopHat, dt, noiseSigma, numberOfNeurons, tMax, theta = currentSettings        # Must be alphabetized, with capitol letters coming first!
        
        yCurrP = np.zeros(numberOfNeurons)
        yCurrN = np.zeros(numberOfNeurons)
        linearRange = np.array(range(0,numberOfNeurons))
        
        xTau = 20
        yTau = 20
        kappa = K**(-1)
        a = kappa*numberOfNeurons*q
        
        xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
        
        p = chopHat*sqrt(xStd*xStd + noiseSigma*noiseSigma)*2*numberOfNeurons*kappa*(rMax-rMin)**(-1);
        
        if FD:
            theta = 1000000000
        crossTimes = 0
        results = 0
        for i in range(perLoc):
            
            yCurrP[:] = 0
            yCurrN[:] = 0
            yCurrPSum = yCurrP.sum()
            yCurrNSum = yCurrN.sum()
        
            beta = myTwister.randNorm(0,betaSigma);

            overTime = 0
            tCurr = 0
            xCurr = myTwister.randNorm(C*.6,xStd)
            xNoise = myTwister.randNorm(0,noiseSigma)

            while yCurrPSum/numberOfNeurons < theta and yCurrNSum/numberOfNeurons < theta:
                
                # Create Input Signal
                xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
                xCurr = xCurr+dt*(C*.6 - xCurr)*xTau**(-1) + xStd*sqrt(2*dt*xTau**(-1))*myTwister.randNorm(mean,std)
                
                # Create Noise Signals
                xNoise = xNoise - dt*xNoise*xTau**(-1) + noiseSigma*sqrt(2*dt*xTau**(-1))*myTwister.randNorm(mean,std)
                
#                # Integrate Preferred Integrator based on chop
#                if fabs((rMax + rMin)*(2*numberOfNeurons)**(-1)*beta + yCurrP*beta + kappa*(xCurr+xNoise)) > (rMax-rMin)*(p-q*beta)*(2*q*numberOfNeurons)**(-1):
#                    yCurrP = yCurrP + dt*yTau**(-1)*(yCurrP*beta + beta*(rMin+rMax)*(2*numberOfNeurons)**(-1) + kappa*(xCurr+xNoise))
#
#                # Integrate Preferred Integrator based on chop                
#                if fabs((rMax + rMin)*(2*numberOfNeurons)**(-1)*beta + yCurrN*beta + kappa*(-(xCurr+xNoise))) > (rMax-rMin)*(p-q*beta)*(2*q*numberOfNeurons)**(-1):
#                    yCurrN = yCurrN + dt*yTau**(-1)*(yCurrN*beta + beta*(rMin+rMax)*(2*numberOfNeurons)**(-1) + kappa*(-(xCurr+xNoise)))

                yCurrPSum = yCurrP.sum()
                yCurrNSum = yCurrN.sum()
                

                
                yCurrP = yCurrP + dt*yTau**(-1)*(-yCurrP + rMin + (rMax-rMin)*np.array(p/q*yCurrP+(1+beta)*(yCurrPSum-yCurrP)+a*(xCurr+xNoise) - ((p/q-1)*(rMax+rMin)/2 + numberOfNeurons*rMin + 
                                                 (rMax - rMin)*(linearRange) + (rMax-rMin)/2)>0,dtype=float))
                
                yCurrN = yCurrN + dt*yTau**(-1)*(-yCurrN + rMin + (rMax-rMin)*np.array(p/q*yCurrN+(1+beta)*(yCurrNSum-yCurrP)+a*(-(xCurr+xNoise)) - ((p/q-1)*(rMax+rMin)/2 + numberOfNeurons*rMin + 
                                                 (rMax - rMin)*(linearRange) + (rMax-rMin)/2)>0,dtype=float))
                                                                                                                                                            

                        
                        

                
                # Update Time Step
                tCurr=tCurr+dt
                
                # Check if we are over time
                if tCurr > tMax:
                    overTime = 1
                    break

            crossTimes += tCurr
            if FD:
                if yCurrPSum/numberOfNeurons > yCurrNSum:
                    results += 1
                elif yCurrPSum/numberOfNeurons == yCurrNSum/numberOfNeurons:
                    tieBreak = myTwister.randNorm(mean,std)
                    if tieBreak > 0:
                        results += 1
            else:
                if not(overTime):
                    if (yCurrPSum/numberOfNeurons >= theta) and (yCurrNSum/numberOfNeurons < theta):
                        results += 1
                    elif yCurrPSum/numberOfNeurons == yCurrNSum/numberOfNeurons:
                        tieBreak = myTwister.randNorm(mean,std)
                        if tieBreak > 0:
                            results += 1
                else:
                    if yCurrPSum/numberOfNeurons > yCurrNSum/numberOfNeurons:
                        results += 1
                    elif yCurrPSum/numberOfNeurons == yCurrNSum/numberOfNeurons:
                        tieBreak = myTwister.randNorm(mean,std)
                        if tieBreak > 0:
                            results += 1
                    
                    


        resultsArray[counter] = results
        crossTimesArray[counter] = crossTimes
        counter += 1

    return (resultsArray, crossTimesArray)