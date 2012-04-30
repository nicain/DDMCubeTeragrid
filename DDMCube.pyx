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
	cdef float xCurr, tCurr, yCurrP, yCurrN, C, xStd, xTau, xNoise, COn, CPost, tFrac, tieBreak,tBeginFrac, tBegin, CNull
	cdef float dt, theta, crossTimes, results, chop, beta, K, yTau, A, B, yBegin, tMax,chopHat, noiseSigma
	cdef double mean = 0, std = 1
	cdef unsigned long mySeed[624]
	cdef c_MTRand myTwister
	cdef int i, overTime
	cdef float chophahahaha
	
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
		A, B, COn, K, beta, chopHat, dt, noiseSigma, tBeginFrac, tFrac, tMax, theta, xStd, xTau, yBegin, yTau = currentSettings		# Must be alphabetized, with capitol letters coming first!


        # Initialize for current loop
		tBegin = tBeginFrac*tMax*(1-tFrac)
		if FD:
			theta = 1000000000
		crossTimes = 0
		results = 0

        # Loop across number of sims, at this point in parameter space
		for i in range(perLoc):
			if tBegin == 0:
				C = COn
			else:
				C = CNull
        
            # Set input standard deviation, pre additive noise
			xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
			overTime = 0
			tCurr = 0
			xCurr = myTwister.randNorm(C*.6,xStd)
			xNoise = myTwister.randNorm(0,noiseSigma)
			yCurrP = yBegin
			yCurrN = yBegin
                    
            # Loop until threshold crossing:
			while yCurrP - yBegin < theta and yCurrN - yBegin < theta:
				
				# Create Input Signal
				if (tCurr < tBegin) or (tCurr > tFrac*tMax):
					C = CNull
				else:
					C = COn
				xStd = sqrt(4.5*((20-.2*C) + (20+.4*C)))
				xCurr = xCurr+dt*(C*.6 - xCurr)*xTau**(-1) + xStd*sqrt(2*dt*xTau**(-1))*myTwister.randNorm(mean,std)
				
				# Create Noise Signals
				xNoise = xNoise - dt*xNoise*xTau**(-1) + noiseSigma*sqrt(2*dt*xTau**(-1))*myTwister.randNorm(mean,std)
				
				# Integrate Preferred Integrator based on chop
				if fabs((xCurr+xNoise) + beta*yCurrP*K + B) < chopHat*sqrt(xStd*xStd + noiseSigma*noiseSigma):
					yCurrP = yCurrP + dt*yTau**(-1)*A
				else:
					yCurrP = yCurrP + dt*yTau**(-1)*((xCurr+xNoise)*K**(-1) + beta*yCurrP + A)

				# Integrate Preferred Integrator based on chop				
				if fabs(-(xCurr+xNoise) + beta*yCurrN*K + B) < chopHat*sqrt(xStd*xStd + noiseSigma*noiseSigma):
					yCurrN = yCurrN + dt*(yTau)**(-1)*A
				else:
					yCurrN = yCurrN + dt*yTau**(-1)*(-(xCurr+xNoise)*K**(-1) + beta*yCurrN + A)
				
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

			crossTimes += tCurr
			if FD:
				if yCurrP > yCurrN:
					results += 1
				elif yCurrP == yCurrN:
					tieBreak = myTwister.randNorm(mean,std)
					if tieBreak > 0:
						results += 1
			else:
				if not(overTime):
					if (yCurrP - yBegin >= theta) and (yCurrN - yBegin < theta):
						results += 1
					elif yCurrP == yCurrN:
						tieBreak = myTwister.randNorm(mean,std)
						if tieBreak > 0:
							results += 1
				else:
                    
                    # A tie happend; break that tie
					if yCurrP > yCurrN:
						results += 1
					elif yCurrP == yCurrN:
						tieBreak = myTwister.randNorm(mean,std)
						if tieBreak > 0:
							results += 1
					
					

        # Record results:
		resultsArray[counter] = results
		crossTimesArray[counter] = crossTimes
		counter += 1

	return (resultsArray, crossTimesArray)