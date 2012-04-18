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
	cdef float xCurr, tCurr, yCurrP, yCurrN, C, xNoise, CPre, CPost, noiseSigma, xStd
	cdef float theta, chop, beta, K, A, B, tMax,chopHat, betaSigma, betaMu
	cdef float xTauInv, yTauInv, KInv
	cdef double mean = 0, std = 1
	cdef unsigned long mySeed[624]
	cdef c_MTRand myTwister
	cdef int i, overTime
    cdef float dt = .1
    cdef float tFrac = 1.0
    cdef float xTau = 20.0
    cdef float yBegin = 0.0
    cdef float yTau = 20.0
	
	# Convert settings dictionary to iterator:
	params = settings.keys()
	params.sort()
	settingsList = []
	totalLength = 1
	for parameter in params: 
		settingsList.append(settings[parameter])
		totalLength *= len(settings[parameter])
	settingsIterator = product.product(*settingsList)
	resultsArray = [0]*totalLength
	crossTimesArray = [0]*totalLength

	# Initialization of random number generator:
	myUUID = uuid.uuid4()
	random.seed(myUUID.int)
	for i in range(624): mySeed[i] = random.randint(0,int(exp(21)))
	myTwister.seed(mySeed)

	# Parameter space loop:
	counter = 0
	CPost = 0
	for currentSettings in settingsIterator:
		C, betaMu, betaSigma, chopHat, noiseSigma, tMax, theta = currentSettings		# Must be alphabetized, with capitol letters coming first!
		
		# Use reciprocal of taus, for speedup:
		xTauInv = 1./xTau
		yTauInv = 1./yTau
		KInv = 1./9.0

		crossTimesArray[counter] = zeros(perLoc)
		resultsArray[counter] = zeros(perLoc)

		counter2 = 0
		if FD:
			theta = 1000000000

		for i in range(perLoc):
			xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
			overTime = 0
			tCurr = 0
			xCurr = myTwister.randNorm(C*.6,xStd)
			xNoise = myTwister.randNorm(0,noiseSigma)
			beta = myTwister.randNorm(betaMu,betaSigma)
			yCurrP = yBegin
			yCurrN = yBegin
			while yCurrP - yBegin < theta and yCurrN - yBegin < theta:
				
				# Create Input Signal
				xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
				xCurr = xCurr+dt*(C*.6 - xCurr)*xTauInv + xStd*sqrt(2*dt*xTauInv)*myTwister.randNorm(mean,std)
				
				# Create Noise Signals
				xNoise = xNoise - dt*xNoise*xTauInv + noiseSigma*sqrt(2*dt*xTauInv)*myTwister.randNorm(mean,std)
				
				# Integrate Preferred Integrator based on chop
				if fabs((xCurr+xNoise) + beta*yCurrP*K) < chopHat*sqrt(xStd*xStd + noiseSigma*noiseSigma):
					yCurrP = yCurrP
				else:
					yCurrP = yCurrP + dt*yTauInv*((xCurr+xNoise)*KInv + beta*yCurrP)

				# Integrate Preferred Integrator based on chop				
				if fabs(-(xCurr+xNoise) + beta*yCurrN*K) < chopHat*sqrt(xStd*xStd + noiseSigma*noiseSigma):
					yCurrN = yCurrN
				else:
					yCurrN = yCurrN + dt*yTauInv*(-(xCurr+xNoise)*KInv + beta*yCurrN)
				
				# Ensure both trains remain positive
				if yCurrP < 0:
					yCurrP = 0
				if yCurrN < 0:
					yCurrN = 0
				
				# Update Time Step
				tCurr=tCurr+dt
				
				# Check if we are over time
				if tCurr > tMax:
					overTime = 1
					break

			crossTimesArray[counter][counter2] = tCurr
			if FD:
				if yCurrP > yCurrN:
					resultsArray[counter][counter2] = 1
			else:
				if not(overTime):
					if (yCurrP - yBegin >= theta) and (yCurrN - yBegin < theta):
						resultsArray[counter][counter2] = 1
				else:
					if yCurrP > yCurrN:
						resultsArray[counter][counter2] = 1
			
			counter2 += 1
		counter += 1

	return (resultsArray, crossTimesArray)