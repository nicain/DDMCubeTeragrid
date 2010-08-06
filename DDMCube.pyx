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
	float abs(float absMe)

################################################################################
######################## Main function, the workhorse:  ########################
################################################################################
def DDMOU(settings, int FD,int perLoc):

	# Import necessary python packages:
	import random, uuid, os, product
	from numpy import zeros
	
	# C initializations
	cdef float xCurr, tCurr, yCurrP, yCurrN, C, xStd, xTau, xNoise, CPre, CPost, tFrac
	cdef float dt, theta, crossTimes, results, chop, beta, K, yTau, A, B, yBegin, tMax,chopHat, noiseSigma
	cdef double mean = 0, std = 1
	cdef unsigned long mySeed[624]
	cdef c_MTRand myTwister
	cdef int i, overTime
	
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
	for i in range(624): mySeed[i] = random.randint(0,2**30)
	myTwister.seed(mySeed)

	# Parameter space loop:
	counter = 0
	CPost = 0
	for currentSettings in settingsIterator:
		A, B, CPre, K, beta, chopHat, dt, noiseSigma, tFrac, tMax, theta, xStd, xTau, yBegin, yTau = currentSettings		# Must be alphabetized, with capitol letters coming first!

		chop = sqrt(xStd*xStd + noiseSigma*noiseSigma)*chopHat

		if FD:
			theta = 1000000000
		crossTimes = np.zeros(perLoc)
		results = np.zeros(perLoc)
		counter2 = 0
		for i in range(perLoc):
			C = CPre
			xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
			overTime = 0
			tCurr = 0
			xCurr = myTwister.randNorm(C*.6,xStd)
			xNoise = myTwister.randNorm(0,noiseSigma)
			yCurrP = yBegin
			yCurrN = yBegin
			while yCurrP - yBegin < theta and yCurrN - yBegin < theta:
				
				# Create Input Signal
				if tCurr > tFrac*tMax:
					C = CPost
				xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
				xCurr = xCurr+dt*(C*.6 - xCurr)/xTau + xStd*sqrt(2*dt/xTau)*myTwister.randNorm(mean,std)
				
				# Create Noise Signals
				xNoise = xNoise - dt*xNoise/xTau + noiseSigma*sqrt(2*dt/xTau)*myTwister.randNorm(mean,std)
				
				# Integrate Preferred Integrator based on chop
				if abs((xCurr+xNoise) + beta*yCurrP*K + B) < chop:
					yCurrP = yCurrP + dt/yTau*A
				else:
					yCurrP = yCurrP + dt/yTau*((xCurr+xNoise)/K + beta*yCurrP + A)

				# Integrate Preferred Integrator based on chop				
				if abs(-(xCurr+xNoise) + beta*yCurrN*K + B) < chop:
					yCurrN = yCurrN + dt/yTau*A
				else:
					yCurrN = yCurrN + dt/yTau*(-(xCurr+xNoise)/K + beta*yCurrN + A)
				
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

			crossTimes[i] += tCurr
			if FD:
				if yCurrP > yCurrN:
					results[i] += 1
			else:
				if not(overTime):
					if (yCurrP - yBegin >= theta) and (yCurrN - yBegin < theta):
						results[i] += 1
				else:
					if yCurrP > yCurrN:
						results[i] += 1
					
					


		resultsArray[counter] = results
		crossTimesArray[counter] = crossTimes
		counter += 1

	return (resultsArray, crossTimesArray)